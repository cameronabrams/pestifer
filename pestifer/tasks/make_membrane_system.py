# Author: Cameron F. Abrams, <cfa2@drexel.edu>
"""
Definition of the :class:`MakeMembraneSystemTask` class for handling embedding proteins into bilayers.

Usage is described in the :ref:`subs_buildtasks_make_membrane_system` documentation.
"""

import copy
import glob
import logging
import numpy as np
from pathlib import Path

from .basetask import BaseTask
# from .terminate import TerminateTask

from ..charmmff.charmmffcontent import CHARMMFFContent

from ..core.artifacts import *
from ..core.errors import PestiferBuildError
from ..util.stringthings import __pestifer_version__

from ..molecule.bilayer import Bilayer, BilayerSpecString, specstrings_builddict

from ..psfutil.psfcontents import get_toppar_from_psf

from ..scripters import PsfgenScripter, PackmolScripter, VMDScripter

from ..util.util import cell_to_xsc,cell_from_xsc, protect_str_arg
from ..util.units import _UNITS_

sA2_ = _UNITS_['SQUARE-ANGSTROMS']

# Fractional lateral-area drift over the final half of a bilayer relaxation above which
# the cell is flagged as not yet equilibrated.  An unconverged calibration cell yields an
# untrustworthy preferred area-per-lipid, which in turn corrupts the stress-free leaflet
# counts of an asymmetric build (the failure mode where an under-relaxed cholesterol
# leaflet fails to condense and its APL comes out spuriously close to a pure-PC leaflet).
_AREA_CONVERGENCE_TOL = 0.005

# Assumed per-leaflet area-compressibility modulus K_A (mN/m), used only to translate a
# measured residual differential stress into a suggested leaflet-count adjustment.  A
# representative value for fluid lipid monolayers; the suggestion is a guide, not a fit.
_DEFAULT_AREA_MODULUS = 250.0

# Minimum half-midplane gap (A) for deterministic grid placement.  The grid packer puts
# opposing-leaflet tail atoms at ``midplane_z +/- half_mid_zgap``; too small a gap makes
# them coincide, and the relaxation minimize then diverges (NaN forces -> segfault).  packmol
# resolves such overlaps with its packing tolerance, so this clamp applies only to grids.
_MIN_GRID_HALF_MID_ZGAP = 1.0

logger = logging.getLogger(__name__)


def _per_leaflet_tension(pp_df, c_z, midplane_slab=None):
    """Per-leaflet surface tension from a NAMD total pressure-profile trajectory.

    ``pp_df`` is the ``pressureprofile`` dataframe produced by the NAMD log parser: a row
    per output step with columns ``x_<i>``, ``y_<i>``, ``z_<i>`` giving the lateral/normal
    pressure components (bar) of z-slab ``i`` (slab 0 lowest).  The kinetic contribution is
    isotropic on average and cancels in ``Pzz - (Pxx+Pyy)/2``, so the total profile suffices.

    The trajectory's second half is averaged (discarding equilibration).  The bilayer is
    assumed centred in the cell, so the midplane is the central slab unless ``midplane_slab``
    is given.  Per-slab tension ``gamma = 0.01 * [Pzz - (Pxx+Pyy)/2] * dz`` (mN/m, with
    ``dz`` in A); leaflet tensions are the slab sums below/above the midplane.

    Returns ``dict(total, lower, upper, dgamma, nslabs, nframes)`` (tensions in mN/m,
    ``dgamma = upper - lower``), or ``None`` if the profile is unavailable/degenerate."""
    if pp_df is None or len(pp_df) == 0:
        return None
    nslabs = sum(1 for col in pp_df.columns if str(col).startswith('z_'))
    if nslabs == 0:
        return None
    frames = pp_df
    if 'TS' in pp_df.columns:                       # drop the step-0 (pre-dynamics) frame
        frames = pp_df[pp_df['TS'] > 0]
    if len(frames) == 0:
        return None
    frames = frames.iloc[len(frames) // 2:]         # average the equilibrated second half
    if len(frames) == 0:
        return None
    pxx = np.array([frames[f'x_{i}'].mean() for i in range(nslabs)])
    pyy = np.array([frames[f'y_{i}'].mean() for i in range(nslabs)])
    pzz = np.array([frames[f'z_{i}'].mean() for i in range(nslabs)])
    dz = float(c_z) / nslabs
    gamma_slab = 0.01 * (pzz - 0.5 * (pxx + pyy)) * dz   # mN/m per slab
    if midplane_slab is None:
        midplane_slab = nslabs // 2
    lower = float(gamma_slab[:midplane_slab].sum())
    upper = float(gamma_slab[midplane_slab:].sum())
    return dict(total=float(gamma_slab.sum()), lower=lower, upper=upper,
                dgamma=upper - lower, nslabs=nslabs, nframes=int(len(frames)))

class MakeMembraneSystemTask(BaseTask):
    """ 
    A class for handling embedding proteins into bilayers
    """
    _yaml_header = 'make_membrane_system'
    """
    YAML header for the MakeMembraneSystemTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """

    @classmethod
    def pipeline_contract(cls, specs):
        from .pipeline_contract import TaskContract, STATE, SOLVATED
        embedding = bool(specs.get('embed'))
        # embedding requires a protein state to insert; a bare membrane build is an origin.
        # The builder solvates the assembled system itself, so a prior solvation is suspect.
        return TaskContract(requires=(STATE,) if embedding else (),
                            provides=(STATE, SOLVATED),
                            discards_state=not embedding,
                            warn_if_present=(SOLVATED,))

    def provision(self, packet: dict):
        logger.debug(f'Provisioning MakeMembraneSystemTask with packet:')
        my_logger(packet, logger.debug)
        super().provision(packet)
        self.patchA: Bilayer = None
        self.patchB: Bilayer = None
        self.patch: Bilayer = None
        self.bilayer_specs: dict = self.specs.get('bilayer', {})
        self.embed_specs: dict = self.specs.get('embed', {})
        self.using_prebuilt_bilayer: bool = False
        self.progress = self.provisions.get('progress-flag', True)
        self.charmmff_content: CHARMMFFContent = self.resource_manager.charmmff_content
    
        if 'prebuilt' in self.bilayer_specs and 'pdb' in self.bilayer_specs['prebuilt']:
            logger.debug('Using prebuilt bilayer')
            self.using_prebuilt_bilayer = True
            self.quilt: Bilayer = Bilayer()
            quilt_state = dict(pdb=self.bilayer_specs['prebuilt']['pdb'], psf=self.bilayer_specs['prebuilt']['psf'], xsc=self.bilayer_specs['prebuilt']['xsc'])
            self.register(quilt_state, key='quilt_state', artifact_type=StateArtifacts)
            self.quilt.box, self.quilt.origin = cell_from_xsc(quilt_state.xsc.path)
            self.quilt.area = self.quilt.box[0][0] * self.quilt.box[1][1]
            additional_topologies = get_toppar_from_psf(quilt_state.psf.name)
            # these will be registered as artifacts when psfgen executes
            self.quilt.addl_streamfiles = additional_topologies
        else:
            self.initialize()

    def initialize(self):
        """
        Initialize the MakeMembraneSystemTask by building the bilayer patch and quilt.
        This method sets up the bilayer patch based on the specifications provided in the configuration.
        If the bilayer is asymmetric, it builds two symmetric patches.
        It also sets up the quilt from the bilayer patch, which will be used for embedding proteins. 
        If a prebuilt bilayer is specified, it uses that instead of building a new one.
        """
        lipid_specstring = self.bilayer_specs.get('lipids', '')
        ratio_specstring = self.bilayer_specs.get('mole_fractions', '')
        conformers_specstring = self.bilayer_specs.get('conformers', '')
        solvent_specstring = self.bilayer_specs.get('solvents', 'TIP3')
        solvent_ratio_specstring = self.bilayer_specs.get('solvent_mole_fractions', '1.0')
        solvent_to_lipid_ratio = self.bilayer_specs.get('solvent_to_lipid_ratio', 32.0)
        patch_nlipids = self.bilayer_specs.get('patch_nlipids', dict(upper=100, lower=100))
        cation_name = self.bilayer_specs.get('cation', 'POT')
        anion_name = self.bilayer_specs.get('anion', 'CLA')
        neutralizing_salt = [cation_name, anion_name]
        salt_con = self.bilayer_specs.get('salt_con', 0.0)  # Molar concentration
        composition_dict = self.bilayer_specs.get('composition', {})

        if not composition_dict['upper_leaflet'] or not composition_dict['lower_leaflet']:
            logger.debug('No upper or lower leaflet specified in composition; building from memgen-format specstrings')
            composition_dict = specstrings_builddict(lipid_specstring,
                                                     ratio_specstring,
                                                     conformers_specstring,
                                                     solvent_specstring,
                                                     solvent_ratio_specstring)
        # Ensure the solvent chambers are present: an explicit `composition:` block specifies
        # only the leaflets, and the asymmetric branch below reads/swaps `upper_chamber`/
        # `lower_chamber` in place.  (Bilayer's constructor fills these in too, but now on its
        # own copy, so it no longer back-populates this shared dict.)
        if 'upper_chamber' not in composition_dict or 'lower_chamber' not in composition_dict:
            C = BilayerSpecString(specstring=solvent_specstring, fracstring=solvent_ratio_specstring)
            composition_dict.setdefault('upper_chamber', C.left)
            composition_dict.setdefault('lower_chamber', C.right)
        logger.debug(f'Naive main composition dict:')
        my_logger(composition_dict, logger.debug)
        self.patch = Bilayer(composition_dict,
                             neutralizing_salt = neutralizing_salt,
                             salt_concentration = salt_con,
                             solvent_specstring = solvent_specstring,
                             solvent_ratio_specstring = solvent_ratio_specstring,
                             solvent_to_key_lipid_ratio = solvent_to_lipid_ratio,
                             leaflet_nlipids = patch_nlipids,
                             charmmffcontent = self.charmmff_content)
        # stash the construction kwargs so the grid packer can rebuild a Bilayer at the
        # full target-box size (direct grid, no patch->quilt tiling). deepcopy the
        # composition because the asymmetric branch below mutates it in place.
        self._bilayer_kwargs = dict(
            composition_dict=copy.deepcopy(composition_dict),
            neutralizing_salt=neutralizing_salt,
            salt_concentration=salt_con,
            solvent_specstring=solvent_specstring,
            solvent_ratio_specstring=solvent_ratio_specstring,
            solvent_to_key_lipid_ratio=solvent_to_lipid_ratio,
        )
        for spdb in self.patch.register_species_pdbs:
            self.register(spdb, key='species_pdb_for_packmol', artifact_type=PDBFileArtifact)
        logger.debug(f'Mature main composition dict:')
        my_logger(composition_dict, logger.debug)
        if self.patch.asymmetric:
            logger.debug(f'Requested patch is asymmetric; generating two symmetric patches')
            logger.debug(f'Symmetrizing bilayer to upper leaflet')
            composition_dict['lower_leaflet_saved'] = composition_dict['lower_leaflet']
            composition_dict['lower_chamber_saved'] = composition_dict['lower_chamber']
            composition_dict['lower_leaflet'] = composition_dict['upper_leaflet']
            composition_dict['lower_chamber'] = composition_dict['upper_chamber']
            self.patchA = Bilayer(composition_dict,
                                  neutralizing_salt = neutralizing_salt,
                                  salt_concentration = salt_con,
                                  solvent_specstring = solvent_specstring,
                                  solvent_ratio_specstring = solvent_ratio_specstring,
                                  solvent_to_key_lipid_ratio = solvent_to_lipid_ratio,
                                  leaflet_nlipids = patch_nlipids,
                                  charmmffcontent = self.charmmff_content)
            logger.debug(f'Symmetrizing bilayer to lower leaflet')
            composition_dict['upper_leaflet_saved'] = composition_dict['upper_leaflet']
            composition_dict['upper_chamber_saved'] = composition_dict['upper_chamber']
            composition_dict['lower_leaflet'] = composition_dict['lower_leaflet_saved']
            composition_dict['lower_chamber'] = composition_dict['lower_chamber_saved']
            composition_dict['upper_leaflet'] = composition_dict['lower_leaflet']
            composition_dict['upper_chamber'] = composition_dict['lower_chamber']
            self.patchB = Bilayer(composition_dict,
                                  neutralizing_salt = neutralizing_salt,
                                  solvent_specstring = solvent_specstring,
                                  solvent_ratio_specstring = solvent_ratio_specstring,
                                  solvent_to_key_lipid_ratio = solvent_to_lipid_ratio,
                                  leaflet_nlipids = patch_nlipids,
                                  charmmffcontent = self.charmmff_content)
            composition_dict['upper_leaflet'] = composition_dict['upper_leaflet_saved']
            composition_dict['upper_chamber'] = composition_dict['upper_chamber_saved']
            self.patch = None

    def do(self) -> int:
        """
        Execute the MakeMembraneSystemTask.
        """
        # as part of a list of tasks, this task expects to be fed a protein system to embed
        protein_state: StateArtifacts = self.get_current_artifact('state')
        if protein_state is None or protein_state.psf is None or protein_state.pdb is None:
            self.embedding = False
        else:
            self.embedding = True
            logger.debug(f'Using psf {protein_state.psf.name} and pdb {protein_state.pdb.name} as inputs')

        if self.embedding:
            self.orient_protein()
        if not self.using_prebuilt_bilayer:
            packer = self.bilayer_specs.get('packer', 'packmol')
            # grid placement fills the whole box directly, retiring the patch->quilt tiling
            # for every case: symmetric builds grid the box directly; asymmetric builds
            # calibrate two symmetric patches then grid at stress-free per-leaflet counts.
            # Either way the box is sized to the protein footprint when embedding. The
            # patch->quilt path now serves only the packmol packer.
            if packer == 'grid':
                if self.patch is not None:
                    self.build_grid_membrane()
                else:
                    self.build_patch()
                    self.build_grid_membrane_asymmetric()
            else:
                self.build_patch()
                self.make_quilt_from_patch()
        if self.embedding:
            self.embed_protein()
        else:
            self.pipeline.rekey('quilt_state', 'state')
        return 0

    def orient_protein(self):
        """
        Orient the protein so that its principal axis is aligned with the z-axis (membrane normal).
        This method uses the Orient package in VMD to perform the orientation based on specified atom selections.
        The oriented protein structure is then saved as a new PDB file and registered as an artifact.
        """
        if self.embed_specs.get('no_orient', True):
            logger.debug('Skipping orientation of protein as per embed specs')
            return
        
        z_head_group = self.embed_specs.get('z_head_group', None)
        z_tail_group = self.embed_specs.get('z_tail_group', None)
        if z_head_group is None or z_tail_group is None:
            logger.debug('No z_head_group or z_tail_group specified; skipping orientation of protein')
            return
        
        self.next_basename('orient')
        logger.debug(f'Orienting protein')
        vm: VMDScripter = self.scripters['vmd']
        vm.newscript(self.basename)
        vm.usescript('bilayer_orient')
        vm.writescript(self.basename)
        self.register(self.basename, key='tcl', artifact_type=VMDScriptArtifact)
        protein_state: StateArtifacts = self.get_current_artifact('state')
        protein_psf: str = protein_state.psf.name
        protein_pdb: str = protein_state.pdb.name
        result = vm.runscript(psf=protein_psf,
                              pdb=protein_pdb,
                              z_head_group=protect_str_arg(z_head_group),
                              z_tail_group=protect_str_arg(z_tail_group),
                              o=self.basename)
        if result != 0:
            raise PestiferBuildError(f'vmd failed with result {result} for {self.basename}')
        self.register(self.basename, key='log', artifact_type=VMDLogFileArtifact)
        self.register(dict(
            psf=protein_psf, 
            pdb=PDBFileArtifact(self.basename, pytestable=True)), 
            key='state', artifact_type=StateArtifacts)

    def build_patch(self):
        """
        Build the bilayer patch or patches based on the specifications provided in the configuration.
        This method retrieves the bilayer specifications, including solution conditions, rotation parameters,
        and other relevant settings.
        It then constructs the patch or patches, packs them using Packmol, generates the PSF file, and then does a short series of equilibration MD simulations.
        """
        logger.debug(f'Bilayer specs:')
        my_logger(self.bilayer_specs, logger.debug)
        solution_gcc: float = self.bilayer_specs.get('solution_gcc',1.0)
        rotation_pm: float = self.bilayer_specs.get('rotation_pm',10.)
        half_mid_zgap: float = self.bilayer_specs.get('half_mid_zgap',1.0)
        SAPL: float = self.bilayer_specs.get('SAPL',75.0)
        seed: int = self.bilayer_specs.get('seed',27021972)
        tolerance: float = self.bilayer_specs.get('tolerance',2.0)
        xy_aspect_ratio: float = self.bilayer_specs.get('xy_aspect_ratio',1.0)
        nloop: int = self.bilayer_specs.get('nloop',100)
        nloop_all: int = self.bilayer_specs.get('nloop_all',100)
        relaxation_protocols: dict = self.bilayer_specs.get('relaxation_protocols',{})
        relaxation_protocol: dict = relaxation_protocols.get('patch',{})
        logger.debug('Relaxation protocols:')
        my_logger(relaxation_protocols, logger.debug)
        packer = self.bilayer_specs.get('packer', 'packmol')
        if packer == 'grid':
            half_mid_zgap = self._grid_half_mid_zgap(half_mid_zgap)
        # we now build the patch, or if asymmetric, two patches
        for patch, specbyte in zip([self.patch, self.patchA, self.patchB], ['', 'A', 'B']):
            if patch is None:
                continue
            patch.spec_out(SAPL=SAPL, xy_aspect_ratio=xy_aspect_ratio,
                            rotation_pm=rotation_pm, solution_gcc=solution_gcc,
                            half_mid_zgap=half_mid_zgap)
            if packer == 'grid':
                self.grid_patch(patch, patch_name=f'patch{specbyte}', seed=seed,
                                half_mid_zgap=half_mid_zgap)
            else:
                self.pack_patch(patch, patch_name=f'patch{specbyte}', seed=seed,
                                tolerance=tolerance,
                                nloop_all=nloop_all,
                                half_mid_zgap=half_mid_zgap,
                                rotation_pm=rotation_pm,
                                nloop=nloop)
            self.do_psfgen(patch, bilayer_name=f'patch{specbyte}')
            self.equilibrate_bilayer(patch, bilayer_name=f'patch{specbyte}', relaxation_protocol=relaxation_protocol)

    def _guard_pierced_lipids(self, protocol):
        """Insert a lipid ``ring_check`` (and a follow-up minimize) before the first dynamics
        stage of a *grid*-packed membrane's relaxation.

        Grid placement itself leaves no pierced rings, but the condensing minimize that opens
        the relaxation can pull a lipid tail through a sterol ring.  Such a threading survives
        further minimization (it is topological, not just strained) and then RATTLE-fails the
        first dynamics step.  ``ring_check`` deletes the pierced lipid on the minimized
        coordinates; the deletion frees the *other* (threading) lipid into a strained pose, so
        a short minimize follows to relax it before any dynamics.  The minimize runs
        regardless, but does real work only when a lipid was removed (otherwise it re-relaxes
        an already-minimized structure).  Only grid builds need this; packmol places lipids
        with a hard overlap tolerance."""
        if not isinstance(protocol, list) or self.bilayer_specs.get('packer') != 'grid':
            return protocol
        guarded, inserted = [], False
        for stage in protocol:
            ensemble = str(stage.get('md', {}).get('ensemble', '') or '').casefold()
            if not inserted and ensemble and ensemble != 'minimize':
                # cutoff 4.0 A: the piercing test itself accepts a bond only within ~3.5 A of
                # the ring centre, so a slightly larger link-cell search radius finds every
                # pierced lipid without the ~10x cost of the schema default (10.0 A), which on
                # a ~350k-atom membrane only adds duplicate reports of the same lipid.
                guarded.append({'ring_check': dict(segtypes=['lipid'], delete='piercee', cutoff=4.0)})
                # relax the lipid freed by the deletion (and the void) before dynamics
                guarded.append({'md': dict(ensemble='minimize', minimize=1000)})
                inserted = True
            guarded.append(stage)
        return guarded

    @staticmethod
    def _autostage_protocol(protocol, chunk_min=500, cap=2000):
        """Split a gridded membrane's long barostatted stages into a ramp of restarts.

        A gridded membrane is placed loose and -- because each leaflet is built at the
        extended lipid length -- much thicker in z than its relaxed state, so it condenses
        substantially (tens of percent) once a barostat is switched on.  NAMD fixes its
        spatial-decomposition patch grid at startup, so a single long run can shrink the
        cell past that grid and abort (``Periodic cell has become too small for original
        patch grid``).  The NAMD-recommended remedy is to **restart** once the cell has
        shrunk -- a fresh run rebuilds the patch grid for the smaller cell -- so we split
        each long NPT/NPgT stage into a ramp of sub-runs (``chunk_min``, ``2*chunk_min``,
        ... up to ``cap`` steps).  Each short sub-run condenses only a little before the
        next restart re-grids, so no single run outgrows its startup grid.  (This is the
        proper fix; raising NAMD's ``margin`` is only an easier-but-slower workaround and is
        not used.)

        The first sub-run is kept at ``chunk_min`` (>> the ~200-step Langevin-piston period)
        so the barostat can actually move the cell: chunks shorter than the piston period
        condense nothing and silently defer all the shrinkage to a later, longer -- and
        overflow-prone -- run.  Non-barostatted stages pass through unchanged.  Returns the
        input untouched when it is empty or not a list."""
        if not protocol or not isinstance(protocol, list):
            return protocol
        staged = []
        for stage in protocol:
            md = stage.get('md', {}) if isinstance(stage, dict) else {}
            ens = str(md.get('ensemble', '')).casefold()
            if ens not in ('npt', 'npat', 'npgt'):   # barostatted stages only
                staged.append(stage)
                continue

            def _emit(steps):
                sub = copy.deepcopy(stage)
                sub['md']['nsteps'] = steps
                staged.append(sub)

            nsteps = md.get('nsteps', 0)
            if not isinstance(nsteps, int) or nsteps <= chunk_min:
                _emit(nsteps)
                continue
            remaining, chunk = nsteps, chunk_min
            while remaining > 0:
                this = min(chunk, remaining)
                _emit(this)
                remaining -= this
                chunk = min(chunk * 2, cap)
        return staged

    def _grid_membrane(self, leaflet_nlipids, box_SAPL, aspect=None):
        """Grid, build, and relax a full-size membrane with the given per-leaflet counts.

        ``box_SAPL`` sizes the box: :meth:`Bilayer.spec_out` sets the lateral area to
        ``box_SAPL * leaflet_nlipids['upper']``.  Pass the upper leaflet's preferred APL as
        ``box_SAPL`` so that area is the stress-free box area for the requested upper count;
        the lower leaflet's (separately chosen) count then sits at its own preferred APL.
        ``aspect`` (Ly/Lx) defaults to the configured ``xy_aspect_ratio``; pass an explicit
        value to match a non-square target box (e.g. a protein footprint).  The result is
        registered as ``quilt_state`` and relaxed with the quilt protocol.
        """
        bs = self.bilayer_specs
        if aspect is None:
            aspect = bs.get('xy_aspect_ratio', 1.0)
        half_mid_zgap = self._grid_half_mid_zgap(bs.get('half_mid_zgap', 1.0))
        membrane = Bilayer(leaflet_nlipids=leaflet_nlipids,
                           charmmffcontent=self.charmmff_content, **self._bilayer_kwargs)
        membrane.spec_out(SAPL=box_SAPL, xy_aspect_ratio=aspect,
                          rotation_pm=bs.get('rotation_pm', 10.0),
                          solution_gcc=bs.get('solution_gcc', 1.0),
                          half_mid_zgap=half_mid_zgap)
        self.quilt = membrane
        for spdb in membrane.register_species_pdbs:
            self.register(spdb, key='species_pdb_for_packmol', artifact_type=PDBFileArtifact)
        self.grid_patch(membrane, patch_name='quilt', seed=bs.get('seed', 27021972),
                        half_mid_zgap=half_mid_zgap)
        self.do_psfgen(membrane, bilayer_name='quilt')
        # the gridded membrane starts loose and condenses; stage the relaxation so the cell
        # never shrinks past NAMD's startup patch grid in a single run
        quilt_protocol = self._autostage_protocol(bs.get('relaxation_protocols', {}).get('quilt', {}))
        self.equilibrate_bilayer(membrane, bilayer_name='quilt', relaxation_protocol=quilt_protocol)

    def _protein_box_dims(self):
        """Lateral (Lx, Ly) the gridded membrane must span to embed the oriented protein:
        the protein's xy extent plus ``embed.xydist`` margin on each side (the same sizing
        the quilt uses, but applied directly so no patch tiling is needed)."""
        margin = self.embed_specs.get('xydist', 10.0)
        state: StateArtifacts = self.get_current_artifact('state')
        xs, ys = [], []
        with open(state.pdb.name) as fh:
            for ln in fh:
                if ln.startswith(('ATOM', 'HETATM')):
                    xs.append(float(ln[30:38]))
                    ys.append(float(ln[38:46]))
        pro_Lx, pro_Ly = max(xs) - min(xs), max(ys) - min(ys)
        Lx, Ly = pro_Lx + 2 * margin, pro_Ly + 2 * margin
        logger.info(f'Sizing gridded membrane to protein {pro_Lx:.1f} x {pro_Ly:.1f} A '
                    f'+ {margin} A margin -> box {Lx:.1f} x {Ly:.1f} A')
        return Lx, Ly

    def _npatch(self):
        """Return the ``[nx, ny]`` patch tiling used to size a non-embedding build,
        defaulting to ``[1, 1]`` when ``npatch`` is unset, empty (the schema default
        for a list type), or malformed."""
        npatch = self.bilayer_specs.get('npatch') or [1, 1]
        if len(npatch) != 2:
            npatch = [1, 1]
        return npatch

    def build_grid_membrane(self):
        """Build a symmetric full-size membrane directly by grid placement, bypassing the
        patch -> quilt tiling.  When embedding a protein the box is sized to the protein
        footprint + margin; otherwise it is sized from ``npatch`` x the patch size.  The
        corresponding number of lipids per leaflet is gridded directly."""
        bs = self.bilayer_specs
        SAPL = bs.get('SAPL', 75.0)
        if self.embedding:
            Lx, Ly = self._protein_box_dims()
            n_target = int(round(Lx * Ly / SAPL))
            aspect = Ly / Lx
        else:
            patch_nlipids = bs.get('patch_nlipids', dict(upper=100, lower=100))
            npatch = self._npatch()
            n_target = int(round(patch_nlipids['upper'] * npatch[0] * npatch[1]))
            aspect = bs.get('xy_aspect_ratio', 1.0)
        logger.debug(f'build_grid_membrane: {n_target} lipids per leaflet (aspect {aspect:.2f})')
        self._grid_membrane(dict(upper=n_target, lower=n_target), SAPL, aspect)

    def build_grid_membrane_asymmetric(self):
        """Build a stress-free asymmetric membrane directly by grid placement.

        ``build_patch`` has already relaxed two symmetric calibration patches -- ``patchA``
        (upper composition) and ``patchB`` (lower composition) -- and recorded their
        equilibrium areas.  Each symmetric patch is tensionless by symmetry, so its
        equilibrium area / lipid count gives the leaflet's preferred area-per-lipid (APL).
        At a common box area both leaflets are tensionless when each is at its own preferred
        APL, i.e. ``n_leaflet = box_area / APL_leaflet``; placing ``n_upper`` and ``n_lower``
        in that ratio gives zero differential stress by construction (to first order; a
        pressure-profile diagnostic can follow).

        When embedding a protein the common box area is sized to the protein footprint +
        margin (so the relaxed membrane is just large enough to embed into); otherwise it is
        sized from the requested upper count (``patch_nlipids`` x ``npatch``).
        """
        bs = self.bilayer_specs
        patch_nlipids = bs.get('patch_nlipids', dict(upper=100, lower=100))
        npatch = self._npatch()
        n_cal = patch_nlipids['upper']                        # per-leaflet count in patchA/patchB
        apl_upper = self.patchA.area / n_cal
        apl_lower = self.patchB.area / n_cal
        # the stress-free counts are only as trustworthy as the calibrated APLs: if either
        # calibration cell had not equilibrated, warn that the resulting membrane may carry
        # residual differential stress (the classic symptom is two leaflets of different
        # composition reporting near-equal APLs because one under-condensed)
        for patch, lname in ((self.patchA, 'upper'), (self.patchB, 'lower')):
            drift = getattr(patch, 'area_drift', None)
            if drift is not None and abs(drift) > _AREA_CONVERGENCE_TOL:
                logger.warning(
                    f'Calibration cell for the {lname} leaflet had not equilibrated '
                    f'(area drift {100 * drift:+.2f}%); its preferred APL '
                    f'({patch.area / n_cal:.2f} {sA2_}) and the resulting stress-free leaflet '
                    f'counts may be unreliable. Lengthen bilayer.relaxation_protocols.patch.')
        if self.embedding:
            # size the common box to the protein footprint and place each leaflet at its own
            # preferred APL there; this keeps the stress-free ratio n_upper/n_lower = apl_lo/apl_up
            Lx, Ly = self._protein_box_dims()
            target_area = Lx * Ly
            aspect = Ly / Lx
            n_upper = int(round(target_area / apl_upper))
            n_lower = int(round(target_area / apl_lower))
        else:
            aspect = bs.get('xy_aspect_ratio', 1.0)
            n_upper = int(round(n_cal * npatch[0] * npatch[1]))   # requested upper count
            n_lower = int(round(n_upper * apl_upper / apl_lower)) # stress-free count ratio
            target_area = n_upper * apl_upper                     # target equilibrium box area
        # Grid ~15% looser than the calibrated preferred APL so the placement has no severe
        # clashes; the NPT relaxation then condenses it to target_area, where each leaflet
        # reaches its own preferred APL. Base the loose target on the *calibrated* APL, not
        # the nominal SAPL: SAPL can badly overestimate a condensed (e.g. cholesterol-rich)
        # leaflet, which would grid the box so loose it has to shrink past NAMD's patch grid.
        box_SAPL = 1.15 * max(apl_upper, apl_lower)
        logger.info(f'Calibrated preferred APL: upper {apl_upper:.2f}, lower {apl_lower:.2f} A^2 '
                    f'(from relaxed symmetric patches)')
        logger.info(f'Stress-free leaflet counts: upper {n_upper}, lower {n_lower} '
                    f'(target equilibrium box area {target_area:.0f} A^2); gridding at SAPL {box_SAPL:.1f}')
        self._grid_membrane(dict(upper=n_upper, lower=n_lower), box_SAPL, aspect)
        diag_specs = self.specs.get('diagnose_differential_stress', {})
        if diag_specs.get('enabled', False):
            self.diagnose_differential_stress(self.quilt, n_upper, n_lower, diag_specs)

    def diagnose_differential_stress(self, membrane: Bilayer, n_upper: int, n_lower: int, diag_specs: dict):
        """Measure and report the residual differential stress of an assembled asymmetric
        membrane.  Diagnostic only -- the membrane is *not* rebuilt.

        Runs a short tensionless (NPgammaT) pressure-profile pass continuing from the
        relaxed membrane, splits the per-slab lateral pressure at the midplane into
        per-leaflet surface tensions (gamma_upper, gamma_lower), and reports their
        difference Dgamma along with the number of lipids that would have to move between
        leaflets to null it (estimated from an assumed area-compressibility modulus).
        Requires a CPU NAMD build (pressureProfile is unsupported on CUDA NAMD)."""
        if self.provisions['processor-type'] == 'gpu':
            logger.warning('diagnose_differential_stress is enabled but processor-type is '
                           '"gpu"; CUDA-enabled NAMD does not support pressureProfile. '
                           'Skipping the differential-stress diagnostic.')
            return
        nsteps = int(diag_specs.get('nsteps', 100000))
        slabs = int(diag_specs.get('slabs', 40))
        freq = int(diag_specs.get('freq', 500))
        K_A = float(diag_specs.get('area_modulus', _DEFAULT_AREA_MODULUS))
        logger.info(f'Differential-stress diagnostic: tensionless pressure-profile pass '
                    f'({nsteps} steps, {slabs} slabs)')
        # The diagnostic continues from the relaxed membrane's PDB (3-decimal precision)
        # with fresh velocities, so a short minimize first reconditions the coordinates;
        # going straight to NPT dynamics off the PDB triggers RATTLE blow-ups.  Pin xstfreq
        # to the profile frequency so the cell-vector and pressure-profile samples share
        # timesteps (mdplot merges them by matching TS).
        pp_protocol = [
            {'md': dict(ensemble='minimize', nsteps=500)},
            {'md': dict(ensemble='NPT', nsteps=nsteps, xstfreq=freq,
                        other_parameters=dict(pressureProfile=True,
                                              pressureProfileSlabs=slabs,
                                              pressureProfileFreq=freq,
                                              surfaceTensionTarget=0.0))}]
        # stage the long pp pass like any other gridded-membrane relaxation: if the membrane
        # is still condensing, the unstaged single run would overflow the patch grid
        self.equilibrate_bilayer(membrane, bilayer_name='quilt',
                                 relaxation_protocol=self._autostage_protocol(pp_protocol))
        mdplot_task = self.subcontroller.tasks[-1]
        pp_df = getattr(mdplot_task, 'dataframes', {}).get('pressureprofile')
        # final cell depth is a good approximation of dz*nslabs for the averaged tail
        c_z = membrane.box[2][2]
        res = _per_leaflet_tension(pp_df, c_z)
        if res is None:
            logger.warning('Differential-stress diagnostic: no usable pressure profile was '
                           'produced; cannot report per-leaflet tensions.')
            return
        membrane.dgamma = res['dgamma']
        logger.info(f'Differential stress (averaged over {res["nframes"]} frames, '
                    f'{res["nslabs"]} slabs):')
        logger.info(f'  gamma_upper          = {res["upper"]:+.3f} mN/m')
        logger.info(f'  gamma_lower          = {res["lower"]:+.3f} mN/m')
        logger.info(f'  gamma_total          = {res["total"]:+.3f} mN/m '
                    f'(sanity: ~0 for a tensionless run)')
        logger.info(f'  Dgamma (upper-lower) = {res["dgamma"]:+.3f} mN/m '
                    f'(target ~0 for a stress-free membrane)')
        # gamma_leaflet ~ K_A*(A/A0 - 1); moving dm lipids from upper to lower changes
        # Dgamma by K_A*dm*(1/n_upper + 1/n_lower), so dm = -Dgamma / [K_A*(...)]
        denom = K_A * (1.0 / n_upper + 1.0 / n_lower)
        dm = -res['dgamma'] / denom if denom else 0.0
        if abs(dm) < 0.5:
            logger.info(f'  Leaflet counts look balanced (suggested move < 1 lipid, '
                        f'assuming K_A = {K_A:.0f} mN/m).')
        else:
            src, dst = ('upper', 'lower') if dm > 0 else ('lower', 'upper')
            logger.info(f'  To reduce Dgamma, move ~{abs(dm):.1f} lipid(s) from the {src} '
                        f'to the {dst} leaflet (assuming K_A = {K_A:.0f} mN/m): rebuild with '
                        f'n_{src} lowered and n_{dst} raised by that count.')

    def register_tops_streams_from_psfgen(self, filelist):
        self.register([CharmmffTopFileArtifact(x) for x in filelist if x.endswith('rtf')], key='charmmff_topfiles', artifact_type=CharmmffTopFileArtifacts)
        self.register([CharmmffStreamFileArtifact(x) for x in filelist if x.endswith('str')], key='charmmff_streamfiles', artifact_type=CharmmffStreamFileArtifacts)

    def register_tmpfiles_from_tcl(self, logname: str):
        tmp_files = {}
        with open(logname, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'register_as_artifact' in line:
                    artifactname = line.split()[-1]
                    ext = os.path.splitext(artifactname)[1].replace('.','')
                    if not ext in tmp_files:
                        tmp_files[ext] = []
                    tmp_files[ext].append(artifactname)
        logger.debug(f'tmp_files: {tmp_files}')
        if 'pdb' in tmp_files:
            self.register([PDBFileArtifact(x) for x in tmp_files['pdb']], key='tmp_patch_pdbs', artifact_type=PDBFileArtifactList)
        if 'psf' in tmp_files:
            self.register([PSFFileArtifact(x) for x in tmp_files['psf']], key='tmp_patch_psfs', artifact_type=PSFFileArtifactList)
        if 'coor' in tmp_files:
            self.register([NAMDCoorFileArtifact(x) for x in tmp_files['coor']], key='tmp_patch_coors', artifact_type=NAMDCoorFileArtifactList)
        if 'log' in tmp_files:
            self.register([NAMDLogFileArtifact(x) for x in tmp_files['log']], key='tmp_patch_logs', artifact_type=NAMDLogFileArtifactList)

    def do_psfgen(self, patch: Bilayer, bilayer_name: str):
        """
        Perform the psfgen operation to generate the PSF and PDB files for the bilayer patch from the packmol output.
        """
        self.next_basename(f'psfgen-{bilayer_name}')
        pg: PsfgenScripter = self.get_scripter('psfgen')
        pg.newscript(self.basename, additional_topologies=patch.addl_streamfiles)
        self.register_tops_streams_from_psfgen(pg.topologies)
        pg.usescript('bilayer_patch')
        pg.writescript(self.basename, guesscoord=False, regenerate=True, force_exit=True)
        state: StateArtifacts = self.get_current_artifact(f'{bilayer_name}_state')
        pdb: Path = state.pdb.path
        result = pg.runscript(pdb=pdb.name, o=self.basename)
        cell_to_xsc(patch.box, patch.origin, f'{self.basename}.xsc')
        patch.area = patch.box[0][0] * patch.box[1][1]
        self.register(dict(
                        psf=PSFFileArtifact(self.basename, pytestable=True), 
                        pdb=PDBFileArtifact(self.basename, pytestable=True), 
                        xsc=NAMDXscFileArtifact(self.basename)), 
                        key=f'{bilayer_name}_state', artifact_type=StateArtifacts)
        self.register(self.basename, key='tcl', artifact_type=PsfgenInputScriptArtifact)
        self.register(self.basename, key='log', artifact_type=PsfgenLogFileArtifact)
        self.register_tmpfiles_from_tcl(f'{self.basename}.log')

    def pack_patch(self, patch: Bilayer, patch_name: str = None, seed=None, tolerance=None, nloop_all=200, nloop=200, half_mid_zgap=1.0, rotation_pm=20):
        """
        Packs the bilayer patch using Packmol.
        
        Parameters
        ----------
        seed : int, optional
            The random seed for the packing process. Default is None.
        tolerance : float, optional
            The tolerance for the packing process. Default is None.
        nloop_all : int, optional
            The total number of loops for the packing process. Default is 200.
        nloop : int, optional
            The number of loops for each individual structure in the packing process. Default is 200.
        half_mid_zgap : float, optional
            The half mid-plane gap in Å. Default is 1.0 Å.
        rotation_pm : float, optional
            The rotation angle in degrees for the patch. Default is 20.0 degrees.
        """
        self.next_basename(f'packmol-{patch_name}')
        pm: PackmolScripter = self.get_scripter('packmol')
        pm.newscript(self.basename)
        packmol_output_pdb = f'{self.basename}.pdb'
        pm.comment(f'packmol input automatically generated by pestifer {__pestifer_version__}')
        pm.addline(f'output {packmol_output_pdb}')
        pm.addline(f'filetype pdb')
        if seed is not None:
            pm.addline(f'seed {seed}')
        pm.addline(f'tolerance {tolerance}')
        pm.addline(f'nloop {nloop_all}')
        patch.write_packmol(pm, half_mid_zgap=half_mid_zgap, rotation_pm=rotation_pm, nloop=nloop)
        pm.writefile()
        result = pm.runscript()
        logger.debug(f'{self.basename} packmol result {result}')
        if result != 0:
            raise PestiferBuildError(f'Packmol failed with result {result}')
        self.register(dict(pdb=PDBFileArtifact(self.basename, pytestable=True)), key=f'{patch_name}_state', artifact_type=StateArtifacts)
        if os.path.exists(f'{self.basename}.pdb_FORCED'):
            self.register(f'{self.basename}.pdb_FORCED', key=f'{patch_name}_packmol_forced', artifact_type=PackMolPDBForcedFileArtifact)
        self.register(self.basename, key='packmol_tcl', artifact_type=PackmolInputScriptArtifact)
        self.register(self.basename, key='packmol_log', artifact_type=PackmolLogFileArtifact)
        if os.path.exists(f'{self.basename}_packmol-results.yaml'):
            self.register(f'{self.basename}_packmol-results.yaml', key=f'{patch_name}_packmol_results', artifact_type=YAMLFileArtifact)
        csvs = glob.glob(f'{self.basename}*_packmol.csv')
        if len(csvs) > 0:
            logger.debug(f'CSVs:')
            my_logger(csvs, logger.debug)
            self.register([CSVDataFileArtifact(x) for x in csvs], key=f'{patch_name}_packmol_csvs', artifact_type=CSVDataFileArtifactList)
        pngs = glob.glob(f'{self.basename}*_packmol.png')
        if len(pngs) > 0:
            logger.debug(f'PNGs:')
            my_logger(pngs, logger.debug)
            self.register([PNGImageFileArtifact(x) for x in pngs], key=f'{patch_name}_packmol_pngs', artifact_type=PNGImageFileArtifactList)

    def _grid_half_mid_zgap(self, configured):
        """Clamp ``half_mid_zgap`` to a minimum safe for deterministic grid placement.

        With too small a gap, opposing-leaflet tail atoms coincide and the relaxation
        minimize diverges (a ``half_mid_zgap: 0`` config -- harmless for packmol -- segfaults
        the grid path).  Warns and uses :data:`_MIN_GRID_HALF_MID_ZGAP` when below it."""
        if configured < _MIN_GRID_HALF_MID_ZGAP:
            logger.warning(
                f'half_mid_zgap={configured} A is too small for grid placement '
                f'(opposing-leaflet tails would overlap); using {_MIN_GRID_HALF_MID_ZGAP} A. '
                f'Set bilayer.half_mid_zgap >= {_MIN_GRID_HALF_MID_ZGAP} to silence this.')
            return _MIN_GRID_HALF_MID_ZGAP
        return configured

    def grid_patch(self, patch: Bilayer, patch_name: str = None, seed=None, half_mid_zgap=1.0):
        """Build a bilayer patch by deterministic grid placement -- a fast, packmol-free
        alternative to :meth:`pack_patch`.

        Produces a coordinate PDB equivalent to packmol's output (lipids on a per-leaflet
        2D lattice, solvent on a 3D lattice) and registers it as the patch state for the
        same downstream psfgen patch build.  Initial overlaps are left for the relaxation
        MD, which is far cheaper than packmol's constrained packing.
        """
        self.next_basename(f'grid-{patch_name}')
        output_pdb = f'{self.basename}.pdb'
        patch.write_grid_pdb(output_pdb, half_mid_zgap=half_mid_zgap, seed=seed)
        self.register(dict(pdb=PDBFileArtifact(self.basename, pytestable=True)),
                      key=f'{patch_name}_state', artifact_type=StateArtifacts)
        logger.debug(f'grid_patch: built {output_pdb} for {patch_name}')
        return 0

    def _area_convergence(self, mdplot_task, name: str):
        """Estimate whether a just-relaxed bilayer's lateral area has plateaued.

        Reads the cell-vector trajectory (``a_x``, ``b_y`` over all relaxation stages)
        that the trailing mdplot task has already parsed and concatenated, and compares
        the mean area over the final quarter of the run to the mean over the preceding
        quarter.  A large systematic change means the area -- and any preferred APL
        calibrated from it -- has not yet converged.  Returns the signed fractional drift
        (second-half minus first-half of the trajectory tail, divided by the mean), or
        ``None`` if the trajectory is unavailable or too short to judge."""
        xst = getattr(mdplot_task, 'dataframes', {}).get('xst')
        if xst is None or not {'a_x', 'b_y'}.issubset(xst.columns):
            logger.debug(f'{name}: no cell-vector trajectory; skipping area-convergence check')
            return None
        area = xst['a_x'].to_numpy() * xst['b_y'].to_numpy()
        # restrict to the final half: skips the initial condensation transient and the
        # short early stages, leaving the long equilibrium-sampling tail
        tail = area[len(area) // 2:]
        if tail.size < 8:
            logger.debug(f'{name}: only {tail.size} cell-vector frames in tail; '
                         'skipping area-convergence check')
            return None
        h = tail.size // 2
        m1, m2 = tail[:h].mean(), tail[h:].mean()
        return float((m2 - m1) / m1)

    def equilibrate_bilayer(self, bilayer: Bilayer, bilayer_name: str, relaxation_protocol: list[dict] = None):
        """
        Equilibrates the bilayer patch using the specified user dictionary and relaxation protocol.
        
        Parameters
        ----------
        bilayer : Bilayer
            The bilayer object to be equilibrated.
        bilayer_name : str
            A string identifier of this bilayer used for keeping track of artifacts.
        relaxation_protocol : list, optional
            A list of dictionaries specifying the stages of the relaxation protocol.
            If not provided, a hard-coded relaxation protocol will be used. 
        """
        self.next_basename(f'equilibration-{bilayer_name}')
        # user_dict=deepcopy(self.config['user'])
        state: StateArtifacts = self.get_current_artifact(f'{bilayer_name}_state')
        logger.debug(f'Bilayer area before equilibration: {bilayer.area:.3f} {sA2_}')
        if not relaxation_protocol:
            logger.debug(f'Using hard-coded relaxation protocol for {self.basename}!!')
            relaxation_protocol=[
                {'md': dict(ensemble='minimize', minimize=1000)},
                {'md': dict(ensemble='NVT', nsteps=1000)},
                {'md': dict(ensemble='NPT', nsteps=200)},
                {'md': dict(ensemble='NPT', nsteps=400)},
                {'md': dict(ensemble='NPT', nsteps=800)},
                {'md': dict(ensemble='NPAT', nsteps=1600)},
                {'md': dict(ensemble='NPAT', nsteps=3200)},
                {'md': dict(ensemble='NPAT', nsteps=6400)},
                {'md': dict(ensemble='NPAT', nsteps=12800)},
                {'md': dict(ensemble='NPAT', nsteps=25600)}]
        else:
            logger.debug(f'Using user-specified relaxation protocol:')
            my_logger(relaxation_protocol, logger.debug)
        pp_specs = self.specs.get('compute_pressure_profile', {})
        # pressure profiling is requested either globally (compute_pressure_profile) or
        # for an individual stage (e.g. the differential-stress diagnostic pass, which sets
        # pressureProfile directly in that stage's other_parameters)
        pp_in_protocol = any(stage.get('md', {}).get('other_parameters', {}).get('pressureProfile')
                             for stage in relaxation_protocol)
        pp_requested = bool(pp_specs.get('enabled', False)) or pp_in_protocol
        is_gpu = self.provisions['processor-type'] == 'gpu'
        if pp_requested and is_gpu:
            raise PestiferBuildError(
                'compute_pressure_profile.enabled is true, but processor-type is "gpu"; '
                'CUDA-enabled NAMD does not support pressureProfile. '
                'Either disable compute_pressure_profile or use a CPU NAMD build.')
        # a grid-packed membrane may get a ring_check (+ minimize) spliced in before its
        # first dynamics stage; build that list now so the inserted minimize is configured
        # (parameter files, flexible cell) alongside the user's own stages
        guarded_protocol = self._guard_pierced_lipids(relaxation_protocol)
        for stage in guarded_protocol:
            if 'md' not in stage:
                continue  # e.g. the spliced-in ring_check
            specs = stage['md']
            specs['addl_paramfiles'] = bilayer.addl_streamfiles
            other = specs.setdefault('other_parameters', {})
            if is_gpu and other.get('pressureProfile'):
                raise PestiferBuildError(
                    f'Stage {specs.get("ensemble", "?")} has pressureProfile set in other_parameters, '
                    f'but processor-type is "gpu"; CUDA-enabled NAMD does not support pressureProfile.')
            if str(specs.get('ensemble', '') or '').casefold() in ('npt', 'npat', 'npgt'):
                other.setdefault('useflexiblecell', True)
                other.setdefault('useconstantratio', True)
                if pp_requested:
                    other.setdefault('pressureProfile', True)
                    other.setdefault('pressureProfileSlabs', pp_specs.get('slabs', 30))
                    other.setdefault('pressureProfileFreq', pp_specs.get('freq', 100))
        timeseries = ['density', ['a_x', 'b_y', 'c_z']]
        profiles = []
        if pp_requested:
            profiles.append('pressure')
            timeseries.append('pressure')  # To do: change this to pressureProfile plotting
        tasklist_user = [{'continuation': dict(psf=state.psf.name, pdb=state.pdb.name, xsc=state.xsc.name)}]
        tasklist_user.extend(guarded_protocol)
        tasklist_user.extend([
            {'mdplot': dict(timeseries=timeseries, profiles=profiles, legend=True, grid=True, basename=self.basename)},
            # {'terminate': dict(basename=self.basename, cleanup=False, chainmapfile=f'{self.basename}-chainmap.yaml', statefile=f'{self.basename}-state.yaml')}
        ])
        subcontroller = self.subcontroller
        subcontroller.config['user']['title'] = f'Bilayer equilibration from {self.basename}'
        subcontroller.reconfigure_tasks(tasklist_user)
        for task in subcontroller.tasks:
            save_task_name = task.taskname
            task_name = f'{self.taskname}-{task.taskname}-{bilayer_name}'
            task.override_taskname(task_name)
            # logger.debug(f'Subcontroller overrides task name {save_task_name} with {task.taskname}')
        subcontroller.do_tasks()
        last_task = subcontroller.tasks[-1]
        bilayer_state: StateArtifacts = last_task.get_current_artifact('state')
        assert bilayer_state is not None
        assert bilayer_state.psf.exists()
        assert bilayer_state.pdb.exists()
        assert bilayer_state.vel.exists()
        assert bilayer_state.xsc.exists()
        assert bilayer_state.coor.exists()
        self.subcontroller.pipeline.rekey('state', f'{bilayer_name}_state')
        self.import_artifacts(subcontroller.pipeline)
        bilayer.box, bilayer.origin = cell_from_xsc(bilayer_state.xsc.name)
        bilayer.area = bilayer.box[0][0] * bilayer.box[1][1]
        logger.debug(f'{self.basename} area after equilibration: {bilayer.area:.3f} {sA2_}')
        # flag (and record on the bilayer) whether its area actually equilibrated, so an
        # asymmetric build can tell whether the preferred APLs it calibrates are trustworthy
        bilayer.area_drift = self._area_convergence(last_task, bilayer_name)
        if bilayer.area_drift is not None:
            pct = 100.0 * bilayer.area_drift
            if abs(bilayer.area_drift) > _AREA_CONVERGENCE_TOL:
                logger.warning(
                    f'{bilayer_name} lateral area still drifting {pct:+.2f}% over the final '
                    f'half of relaxation (tolerance {100 * _AREA_CONVERGENCE_TOL:.1f}%); its '
                    f'equilibrium area may be unconverged. Consider lengthening the '
                    f'relaxation protocol.')
            else:
                logger.debug(f'{bilayer_name} area converged (final-half drift {pct:+.2f}%)')

    def make_quilt_from_patch(self):
        """
        Create a quilt from the bilayer patch or patches.
        This method generates a quilt that combines the bilayer patches into a single structure.
        It uses the psfgen scripter to create a script that builds the quilt based on the provided patches.
        The quilt is then equilibrated and saved with the appropriate state variables.
        """
        
        logger.debug(f'Creating quilt from patch')
        self.next_basename('quilt')
        additional_topologies=[]
        if self.patch is not None:
            patch_state: StateArtifacts = self.get_current_artifact('patch_state')
            pdb = patch_state.pdb.name
            xsc = patch_state.xsc.name
            psf = patch_state.psf.name
            pdbA = pdbB = pdb
            psfA = psfB = psf
            xscA = xscB = xsc
        elif self.patchA is not None and self.patchB is not None:
            patchA_state: StateArtifacts = self.get_current_artifact('patchA_state')
            psfA = patchA_state.psf.name
            pdbA = patchA_state.pdb.name
            xscA = patchA_state.xsc.name
            patchB_state: StateArtifacts = self.get_current_artifact('patchB_state')
            psfB = patchB_state.psf.name
            pdbB = patchB_state.pdb.name
            xscB = patchB_state.xsc.name
        else:
            raise PestiferBuildError("No valid patch state found.")
        for patch in [self.patch, self.patchA, self.patchB]:
            if patch is None:
                continue
            additional_topologies += patch.addl_streamfiles
        additional_topologies = list(set(additional_topologies))
        pg: PsfgenScripter = self.get_scripter('psfgen')
        pg.newscript(self.basename, additional_topologies=additional_topologies)
        self.register_tops_streams_from_psfgen(pg.topologies)
        pg.usescript('bilayer_quilt')
        pg.writescript(self.basename, guesscoord=False, regenerate=False, force_exit=True, writepsf=False, writepdb=False)
        self.register(self.basename, key='tcl', artifact_type=PsfgenInputScriptArtifact)
        margin = self.embed_specs.get('xydist', 10.0)
        protein_state: StateArtifacts = self.get_current_artifact('state')
        if protein_state is not None and protein_state.pdb.exists():
            # we will eventually embed a protein in here, so send its pdb along to help size the bilayer
            result = pg.runscript(propdb=protein_state.pdb.name, margin=margin, psfA=psfA, pdbA=pdbA, psfB=psfB, pdbB=pdbB, xscA=xscA, xscB=xscB, o=self.basename)
        else:
            dims = self.bilayer_specs.get('dims', [0,0])
            if dims is None or len(dims) != 2:
                dims = [0,0]
            dimx, dimy = dims
            npatch = self.bilayer_specs.get('npatch', [0,0])
            if npatch is None or len(npatch) != 2:
                npatch = [0,0]
            npatchx, npatchy = npatch
            if npatchx != 0 and npatchy != 0:
                result = pg.runscript(nx=npatchx, ny=npatchy, psfA=psfA, pdbA=pdbA,
                                      psfB=psfB, pdbB=pdbB, xscA=xscA, xscB=xscB, o=self.basename)
            elif dimx != 0 and dimy != 0:
                result = pg.runscript(dimx=dimx, dimy=dimy, psfA=psfA, pdbA=pdbA,
                                      psfB=psfB, pdbB=pdbB, xscA=xscA, xscB=xscB, o=self.basename)
        if result != 0:
            raise PestiferBuildError(f'psfgen failed with result {result} for {self.basename}')
        self.register(self.basename, key='log', artifact_type=PsfgenLogFileArtifact)
        self.register_tmpfiles_from_tcl(f'{self.basename}.log')
        self.register(dict(
            pdb=PDBFileArtifact(self.basename, pytestable=True), 
            psf=PSFFileArtifact(self.basename, pytestable=True), 
            xsc=NAMDXscFileArtifact(self.basename)), key='quilt_state', artifact_type=StateArtifacts)
        self.quilt = Bilayer()
        self.quilt.addl_streamfiles = additional_topologies
        self.quilt.box, self.quilt.origin = cell_from_xsc(f'{self.basename}.xsc')
        self.quilt.area = self.quilt.box[0][0] * self.quilt.box[1][1]
        relaxation_protocol = self.bilayer_specs.get('relaxation_protocols', {}).get('quilt', {})
        self.equilibrate_bilayer(self.quilt,
                                 bilayer_name='quilt',
                                 relaxation_protocol=relaxation_protocol)

    def embed_protein(self):
        """
        Embed the protein into the bilayer patch or quilt.
        This method uses the psfgen scripter to create a script that embeds the protein
        into the bilayer based on the provided specifications.
        It retrieves the embedding specifications, such as head and tail groups, reference groups,
        and orientation settings, and runs the script to perform the embedding.
        The resulting state variables are updated with the new PSF, PDB, and coordinate files.
        """
        if not self.embed_specs:
            logger.debug('No embed specs.')
            return
        zvals = np.zeros(2)
        dum_pdb: PDBFileArtifact = self.get_current_artifact('base_coordinates_dum')
        if dum_pdb and dum_pdb.exists():
            logger.debug(f'Using DUM pdb {dum_pdb.name} for embedding coordinates')
            z_ref_group = ""
            dum_struct = PDBParser(filepath=dum_pdb.path).parse().parsed
            zvals = np.array(list(set(a.z for a in dum_struct['HETATM'])))
            logger.debug(f'DUM z-values: {zvals}')
            assert len(zvals) == 2, "DUM pdb must have exactly two distinct z-values for head and tail groups"
            z_value = zvals.mean()
        else:
            z_ref_group = self.embed_specs.get('z_ref_group', {}).get('text', None)
            z_value = self.embed_specs.get('z_ref_group', {}).get('z_value', 0.0)
        self.next_basename('embed')
        pg: PsfgenScripter = self.scripters['psfgen']
        pg.newscript(self.basename, additional_topologies=self.quilt.addl_streamfiles)
        pg.usescript('bilayer_embed')
        # regenerate=False: the bilayer_embed script does its own regenerate+writepsf for
        # every structure it writes (prefill, filled, and the autoionized final psf). A
        # trailing 'regenerate angles dihedrals' appended here would run after the final
        # files are already written, with no consumer -- pure wasted work on a ~2M-atom
        # psfgen context.
        pg.writescript(self.basename, guesscoord=False, regenerate=False, force_exit=True, writepsf=False, writepdb=False)
        self.register(self.basename, key='tcl', artifact_type=PsfgenInputScriptArtifact)
        quilt_state: StateArtifacts = self.get_current_artifact('quilt_state')
        bilayer_psf: str = quilt_state.psf.name
        bilayer_pdb: str = quilt_state.pdb.name
        bilayer_xsc: str = quilt_state.xsc.name
        protein_state: StateArtifacts = self.get_current_artifact('state')
        protein_psf: str = protein_state.psf.name
        protein_pdb: str = protein_state.pdb.name
        logger.debug(f'Embedding {protein_pdb} with z_ref_group {z_ref_group} at z={z_value:.3f} into bilayer {bilayer_pdb}')
        result = pg.runscript(psf=protein_psf,
                              pdb=protein_pdb,
                              bilayer_psf=bilayer_psf,
                              bilayer_pdb=bilayer_pdb,
                              bilayer_xsc=bilayer_xsc,
                              z_ref_group=protect_str_arg(z_ref_group),
                              z_lo_dum=zvals[0],
                              z_hi_dum=zvals[1],
                              z_value=z_value,
                              o=self.basename)
        if result != 0:
            raise PestiferBuildError(f'psfgen failed with result {result} for {self.basename}')
        self.register(self.basename, key='log', artifact_type=PsfgenLogFileArtifact)
        self.register(dict(
            psf=PSFFileArtifact(self.basename, pytestable=True), 
            pdb=PDBFileArtifact(self.basename, pytestable=True), 
            xsc=NAMDXscFileArtifact(self.basename)), key='state', artifact_type=StateArtifacts)
        self.register_tmpfiles_from_tcl(f'{self.basename}.log')
        logger.debug(f'Embedding completed with result {result}')
        return result

    
