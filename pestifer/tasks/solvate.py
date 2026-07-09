# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`SolvateTask` class for solvating a molecular structure.
This class is a descendant of the :class:`BaseTask <pestifer.tasks.basetask.BaseTask>` class and is used to solvate a molecular structure
using the VMD solvate and autoionize packages.
It generates a solvated PDB and PSF file, optionally adding ions based on specified salt concentration and ion types.
The task can also handle cubic or rectangular boxes for solvation based on the provided specifications.
The task reads the input PDB and PSF files, calculates the bounding box for the solvation,
and generates the necessary Tcl commands to perform the solvation and ionization.
The resulting solvated structure is saved with the specified basename, and the state is updated with the new PSF, PDB, and XSC files.

Usage is described in the :ref:`subs_buildtasks_solvate` documentation.
"""
import logging
import os

import numpy as np

from pidibble.pdbparse import PDBParser

from .basetask import VMDTask
from ..core.artifacts import *
from ..core.errors import PestiferError
from ..psfutil.psfcontents import PSFContents
from ..util.util import cell_from_xsc, cell_to_xsc
from ..scripters import VMDScripter

logger = logging.getLogger(__name__)

class SolvateTask(VMDTask):
    """
    SolvateTask class for solvating a molecular structure.
    """
    _yaml_header = 'solvate'
    """
    YAML header for the SolvateTask, used to identify the task in configuration files as part of a ``tasks`` list.
    This header is used to declare SolvateTask objects in YAML task lists.
    """

    @classmethod
    def pipeline_contract(cls, specs):
        from .pipeline_contract import TaskContract, STATE, SOLVATED
        # adds bulk solvent; warns if the system is already solvated
        return TaskContract(requires=(STATE,), provides=(STATE, SOLVATED), warn_if_present=(SOLVATED,))

    # solvent names that mean "use VMD's built-in pre-equilibrated water box" (no custom box)
    _WATER_SOLVENTS = {'TIP3', 'TIP3P', 'WATER', 'WAT'}

    def _ionization_plan(self, solvent):
        """
        Return ``(wants_ions, is_water, neutralize)``.  Ions are added when the system is to
        be neutralized (``neutralize``, default True) or a salt concentration is requested.
        Water uses ``autoionize``; a non-water solvent is ionized by solvent replacement.
        """
        is_water = solvent is None or solvent.upper() in self._WATER_SOLVENTS
        neutralize = bool(self.specs.get('neutralize', True))
        salt = self.specs.get('salt_con') is not None
        return (neutralize or salt), is_water, neutralize

    def _solvent_box_entry(self, solvent: str):
        """
        Return the ``kind: box`` PDB-repository entry for a non-water ``solvent`` (checking it
        out), or ``None`` for water.  Raises if the solvent has no box entry or is only a
        single-molecule entry.
        """
        if solvent is None or solvent.upper() in self._WATER_SOLVENTS:
            return None
        CC = self.resource_manager.charmmff_content
        if CC.pdbrepository is None:
            CC.provision_pdbrepository()
        repo = CC.pdbrepository
        if solvent not in repo:
            self._generate_solvent_box(solvent, CC, repo)
        entry = repo.checkout(solvent)
        if not entry.is_box():
            raise PestiferError(
                f"solvent '{solvent}' is a single-molecule PDB entry, not a pre-equilibrated box; "
                f"solvate needs a kind:box entry (build one with 'make-pdb-collection solvent')")
        return entry

    def _generate_solvent_box(self, solvent: str, CC, repo):
        """Handle a solvent-box miss: unless generation is disabled, build the box on the fly,
        cache it under ``~/.pestifer/``, and register the cache collection so ``repo`` can find it.

        Raises :class:`PestiferError` (preserving the former behavior) when generation is disabled
        (``charmmff.generate_missing_coordinates: false``) or the solvent is not defined in the
        force field at all.
        """
        if not getattr(CC, 'generate_missing_coordinates', True):
            raise PestiferError(
                f"solvent '{solvent}' is not in the PDB repository; build a box for it with "
                f"'pestifer make-pdb-collection solvent --resname {solvent}' and install it into "
                f"the solvent collection, or set 'charmmff.generate_missing_coordinates: true' to "
                f"have pestifer build and cache one automatically")
        if solvent not in CC:
            raise PestiferError(
                f"solvent '{solvent}' is neither in the PDB repository nor defined in the CHARMM "
                f"force field; cannot auto-generate a box")
        from ..charmmff.autocache import ensure_solvent_box
        collection_dir = ensure_solvent_box(solvent, CC)
        repo.add_resource(str(collection_dir))
        if solvent not in repo:
            raise PestiferError(
                f"internal error: auto-generated solvent box for '{solvent}' but it did not "
                f"register in the PDB repository (cache at {collection_dir})")

    def _write_solvent_topology(self, solvent: str, outpath: str) -> str:
        """
        Write a small self-contained CHARMM topology for ``solvent`` -- the ``MASS`` records for
        the atom types its RESI uses (each atom knows its type/mass/element) followed by the RESI
        block -- so VMD solvate's isolated psfgen context can build the solvent's replica residues
        even when the solvent's stock topfile (e.g. a model-compound stream) omits those masses.
        Returns ``outpath``.
        """
        import io
        CC = self.resource_manager.charmmff_content
        CC.provision()
        topo = CC.get_resi(solvent)
        masses = {}
        for a in topo.atoms:
            masses.setdefault(a.type, (a.mass, a.element))
        resi = io.StringIO()
        topo.to_file(resi)
        with open(outpath, 'w') as f:
            f.write(f'* self-contained topology for {solvent} (pestifer-generated)\n*\n   36  1\n\n')
            for t, (m, el) in masses.items():
                f.write(f'MASS  -1  {t:<8s}{m:>10.5f} {el}\n')
            # CGenFF RESIs list only atoms and bonds; angles/dihedrals are auto-generated.  Without
            # this directive, solvate's replica `segment {residue <solvent>}` (which sets no `auto`
            # option) builds the molecules with no angle/dihedral terms and they distort under MD.
            f.write('\nAUTOGENERATE ANGLES DIHEDRALS\n\n')
            f.write(resi.getvalue())
            f.write('\nEND\n')
        return outpath

    def _solvent_box_args(self, entry, solvent: str) -> str:
        """
        Build the trailing ``solvate`` arguments (``-spsf/-spdb/-ws/-ks``) that tile the given
        box ``entry``, or ``''`` for water (``entry is None``).
        """
        if entry is None:
            return ''
        spsf, spdb, ws, ks = (entry.get_box_psf(), entry.get_box_pdb(),
                              entry.get_box_edge(), entry.get_key_atom())
        logger.info(f"solvating with pre-equilibrated {solvent} box (edge {ws} A, key atom {ks})")
        args = f' -spsf {spsf} -spdb {spdb} -ws {ws}'
        if ks:
            args += f' -ks "name {ks}"'
        return args

    # ion valences for neutralization/salt counts (CHARMM ion RESIs are monatomic; charge is
    # taken from the PDB repository when available, else this fallback map)
    _ION_CHARGE = {'SOD': 1, 'POT': 1, 'LIT': 1, 'CES': 1, 'RUB': 1,
                   'CLA': -1, 'CAL': 2, 'MG': 2, 'ZN2': 2, 'BAR': 2, 'CD2': 2}

    def _ion_charge(self, resname):
        repo = getattr(self.resource_manager.charmmff_content, 'pdbrepository', None)
        if repo is not None and resname in repo:
            q = repo.checkout(resname).get_charge()
            if q:
                return q
        return self._ION_CHARGE.get(resname.upper())

    def _ion_counts(self, net_charge, box_volume_A3, cation, anion, salt_con, neutralize=True):
        """
        Return ``{resname: count}`` of ions to add: ``salt_con`` (mol/L) worth of
        cation/anion pairs sized from the box volume, plus (when ``neutralize``) enough
        counter-ions to cancel ``net_charge``.  Ion valences come from :meth:`_ion_charge`.
        """
        qcat = abs(self._ion_charge(cation) or 1)
        qani = abs(self._ion_charge(anion) or 1)
        counts = {}
        if salt_con:
            n_pairs = int(round(salt_con * box_volume_A3 * 1e-27 * 6.02214076e23))
            counts[cation] = counts.get(cation, 0) + n_pairs
            counts[anion] = counts.get(anion, 0) + n_pairs
        if neutralize:
            net = int(round(net_charge))
            if net > 0:
                counts[anion] = counts.get(anion, 0) + int(round(net / qani))
            elif net < 0:
                counts[cation] = counts.get(cation, 0) + int(round(-net / qcat))
        return {k: v for k, v in counts.items() if v > 0}

    def do(self):
        """
        Execute the solvate task.
        This method initializes the task, inherits the state from prior tasks, and performs the solvation and ionization.
        """
        self.next_basename()
        state: StateArtifacts = self.get_current_artifact('state')
        # psf: Path = self.get_current_artifact_path('psf')
        # pdb: Path = self.get_current_artifact_path('pdb')
        # xsc: Path = self.get_current_artifact_path('xsc')
        # self.stash_current_artifact('vel')
        use_minmax = True
        if state.xsc is not None:
            box, origin = cell_from_xsc(state.xsc.name)
            if box is not None and origin is not None:
                use_minmax = False
                basisvec = box.diagonal()
                LL = origin - 0.5 * basisvec
                UR = origin + 0.5 * basisvec
            xsc = state.xsc
        if use_minmax:
            p = PDBParser(filepath=state.pdb.name).parse()
            x = np.array([a.x for a in p.parsed['ATOM']])
            y = np.array([a.y for a in p.parsed['ATOM']])
            z = np.array([a.z for a in p.parsed['ATOM']])
            minmax = np.array([[x.min(), y.min(), z.min()], [x.max(), y.max(), z.max()]])
            spans = minmax[1] - minmax[0]
            maxspan = spans.max()
            cubic = self.specs.get('cubic', False)
            sympad = np.array([0., 0., 0.])
            if cubic:
                sympad = 0.5 * (maxspan - spans)
            pad = self.specs.get('pad', 10.0)
            LL = minmax[0] - pad - sympad
            UR = minmax[1] + pad + sympad
            basisvec = UR - LL
            origin = 0.5 * (UR + LL)
            cell_to_xsc(np.diag(basisvec), origin, f'{self.basename}.xsc')
            xsc = NAMDXscFileArtifact(self.basename)

        ll_tcl = r'{ ' + ' '.join([str(_) for _ in LL.tolist()]) + r' }'
        ur_tcl = r'{ ' + ' '.join([str(_) for _ in UR.tolist()]) + r' }'
        box_tcl = r'{ ' + ll_tcl + ' ' + ur_tcl + r' }'

        # solvent species: TIP3/water -> VMD's built-in water box; anything else -> a
        # pre-equilibrated kind:box entry from the 'solvent' PDB collection (supplied to
        # solvate as -spsf/-spdb/-ws/-ks)
        solvent = self.specs.get('solvent', 'TIP3')
        box_entry = self._solvent_box_entry(solvent)
        box_args = self._solvent_box_args(box_entry, solvent)

        # For a non-water solvent, solvate tiles the box by building `segment {residue <solvent>}`
        # replicas in its own (isolated) psfgen context, which needs the solvent's topology --
        # the RESI block *and* the atom-type MASS records it references.  A solvent's stock topfile
        # (e.g. a model-compound stream) may not carry those masses, so hand solvate a small
        # self-contained topology we generate for just this solvent (-stop).
        if box_entry is not None:
            stop = self._write_solvent_topology(solvent, f'{self.basename}-{solvent}.rtf')
            box_args += f' -stop {stop}'

        # Ionization.  For water, VMD autoionize replaces water molecules with ions.  For a
        # non-water solvent, autoionize cannot work (it needs water to make room), so we ionize
        # only when the user explicitly asked -- by replacing solvent molecules with ions
        # (PestiferIonize).  With no ions requested for a non-water solvent, the system keeps
        # its net charge (NAMD neutralizes with a uniform background under PME).
        sc = self.specs.get('salt_con', None)
        cation = self.specs.get('cation', None) or 'SOD'
        anion = self.specs.get('anion', None) or 'CLA'
        wants_ions, is_water, neutralize = self._ionization_plan(solvent)
        replace_ionize = wants_ions and not is_water

        ion_counts, ion_topfile = {}, None
        if replace_ionize:
            CC = self.resource_manager.charmmff_content
            net_charge = sum(a.charge for a in PSFContents(state.psf.name).atoms)
            box_volume = float(abs(np.prod(basisvec)))
            ion_counts = self._ion_counts(net_charge, box_volume, cation, anion, sc, neutralize)
            if ion_counts:
                topbn = (CC.get_topfile_of_resname(anion) or CC.get_topfile_of_resname(cation)
                         or 'toppar_water_ions.str')
                ion_topfile = CC.copy_charmmfile_local(topbn)

        packages = ['PestiferIonize'] if (replace_ionize and ion_counts) else []
        vt: VMDScripter = self.scripters['vmd']
        vt.newscript(self.basename, packages=packages)
        vt.addline( 'package require solvate')
        vt.addline( 'package require autoionize')
        vt.addline( 'psfcontext mixedcase')
        vt.addline(f'mol new {state.psf.name}')
        vt.addline(f'mol addfile {state.pdb.name} waitfor all')
        vt.addline(f'solvate {state.psf.name} {state.pdb.name} -minmax {box_tcl}{box_args} -o {self.basename}_solv')

        if is_water and wants_ions:
            # water: VMD autoionize (-sc also neutralizes; else -neutralize)
            ai_args = [f'-sc {sc}'] if sc is not None else ['-neutralize']
            if self.specs.get('cation'):
                ai_args.append(f"-cation {self.specs['cation']}")
            if self.specs.get('anion'):
                ai_args.append(f"-anion {self.specs['anion']}")
            vt.addline(f'autoionize -psf {self.basename}_solv.psf -pdb {self.basename}_solv.pdb {" ".join(ai_args)} -o {self.basename}')
        elif replace_ionize and ion_counts:
            # non-water: replace solvent molecules with ions (autoionize can't ionize a
            # non-aqueous box; see PestiferIonize)
            ks = box_entry.get_key_atom()
            ionlist = ' '.join(f'{rn} {cnt}' for rn, cnt in ion_counts.items())
            logger.info(f'ionizing {solvent} box by solvent replacement: {ion_counts}')
            vt.addline(f'PestiferIonize::ionize_by_replacement {self.basename}_solv.psf '
                       f'{self.basename}_solv.pdb -o {self.basename} -solvent {solvent} '
                       f'-keyatom {ks} -ions {{{ionlist}}} -topology {ion_topfile}')
        else:
            if not wants_ions:
                logger.info('ionization disabled (neutralize=false, no salt); the system keeps '
                            'its net charge')
            else:
                logger.info(f'{solvent} system is already net-neutral; no ions needed')
            # promote the solvate output to the task basename so downstream naming is uniform
            vt.addline(f'file copy -force {self.basename}_solv.psf {self.basename}.psf')
            vt.addline(f'file copy -force {self.basename}_solv.pdb {self.basename}.pdb')

        vt.writescript()
        self.register(self.basename, key='tcl', artifact_type=VMDScriptArtifact)
        self.result = vt.runscript(progress_title='solvate', progress_color='fluorescent_blue')
        self.register(self.basename, key='log', artifact_type=VMDLogFileArtifact)
        self.register(self.basename+'_solv', key='log', artifact_type=VMDLogFileArtifact)
        if self.result != 0:
            return self.result
        # VMD exits 0 even when solvate/autoionize catch their own errors, so verify the
        # expected outputs were actually produced before continuing
        for ext in ('psf', 'pdb'):
            if not os.path.isfile(f'{self.basename}.{ext}'):
                if is_water and wants_ions:
                    hint = 'autoionize failed to place ions; check the VMD log'
                elif replace_ionize:
                    hint = ('solvent-replacement ionization failed to place all ions; lower the '
                            'ion count/concentration or use a larger box')
                else:
                    hint = 'solvate produced no output; check the VMD log'
                raise PestiferError(
                    f'solvate task {self.taskname}: no {self.basename}.{ext} was produced '
                    f'(VMD returned {self.result}) -- {hint}. See {self.basename}.log')
        psf1 = PSFFileArtifact(f'{self.basename}_solv.psf')
        pdb1 = PDBFileArtifact(f'{self.basename}_solv.pdb')
        self.register(dict(psf=psf1, pdb=pdb1, xsc=xsc), key='state', artifact_type=StateArtifacts)
        self.pdb_to_coor(f'{self.basename}.pdb')
        self.register(dict(pdb=PDBFileArtifact(f'{self.basename}', pytestable=True), coor=NAMDCoorFileArtifact(f'{self.basename}'), psf=PSFFileArtifact(f'{self.basename}', pytestable=True), xsc=xsc), key='state', artifact_type=StateArtifacts)
        return self.result