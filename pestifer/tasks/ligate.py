# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`LigateTask` class for ligating loops in molecular dynamics simulations.
This class is a descendant of the :class:`MDTask <pestifer.tasks.mdtask.MDTask>` class and is used to ligate loops
in a molecular structure using the NAMD molecular dynamics engine.
It measures the distances between loop termini, steers them toward each other, and connects them using a specified patch.
The resulting structure is then saved as a PSF/PDB files.

Usage is described in the :ref:`subs_buildtasks_ligate` documentation.
"""
import logging

from pathlib import Path

from .mdtask import MDTask

from ..charmmff.charmmffcontent import CHARMMFFContent
from ..core.artifacts import *
from ..core.resourcemanager import ResourceManager
from ..molecule.molecule import Molecule
from ..scripters import GenericScripter, VMDScripter, PsfgenScripter, NAMDColvarInputScripter

logger = logging.getLogger(__name__)

class LigateTask(MDTask):
    """
    LigateTask class for ligating loops in molecular dynamics simulations.
    """
    _yaml_header = 'ligate'
    """
    YAML header for the LigateTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """

    @classmethod
    def pipeline_contract(cls, specs):
        from .pipeline_contract import TaskContract, STATE, MOLECULE
        # closes protein loops on the in-memory molecule built by psfgen
        return TaskContract(requires=(MOLECULE,), provides=(STATE, MOLECULE))

    def do(self) -> int:
        """
        Execute the ligate task. This method checks if the base molecule has loops,
        measures the distances between loop termini, steers them toward each other,
        and connects them using a specified patch. The resulting structure is saved as a PSF/PDB fileset.
        It also writes the gaps to a data file and measures the distances between the termini of the loops.
        If the base molecule does not have loops, the task is bypassed.
        The method returns the result of the NAMD run, which is 0 on success or a non-zero error code on failure.
        If the task is successful, it saves the state of the simulation with the specified extensions.
        If the task is bypassed, it logs a message and returns without performing any operations.
        """
        self.base_molecule: Molecule = self.get_current_artifact_data('base_molecule')
        if not self.base_molecule.has_protein_loops:
            self.log_message('bypassed')
            return
        method = str(self.specs.get('method', 'ccd')).casefold()
        if method == 'steer':
            logger.debug('Storing sequence gaps.')
            self.write_gaps()
            logger.debug('Measuring gap distances.')
            steering_specs = self.specs.get('steer', {})
            if not steering_specs:
                logger.debug(f'No steering specifications for ligate task; this is a bug; bypassing')
                return
            self.measure_distances(steering_specs)
            logger.debug('Steering loop C-termini toward their partner N-termini')
            self.result = self.do_steered_md(self.specs['steer'])
        else:
            logger.debug('Closing loops onto their downstream anchors by cyclic coordinate descent')
            self.result = self.close_loops_ccd(self.specs.get('ccd', {}))
        if self.result != 0:
            return self.result
        logger.debug('Connecting loop C-termini to their partner N-termini')
        connect_specs = self.specs.get('connect', {})
        self.result = self.connect(connect_specs)
        # Relaxation is intentionally left to an explicit downstream `minimize`/`md` task
        # (the standard pattern is to follow `ligate` with one). The clash-filtered ensemble
        # already hands off a closed, clash-free structure ready for that minimization.
        return self.result

    def close_loops_ccd(self, ccd_specs: dict) -> int:
        """
        Close each modeled interior loop onto its downstream anchor by cyclic coordinate
        descent (offline, deterministic). For each gap the loop's C-terminal backbone is
        walked onto an ideal peptide-bond target built from the anchor, distributing the
        closure across the loop's own phi/psi dihedrals (whole residues, sidechains included).

        These gaps are floppy, solvent-exposed surface loops (that is *why* they are
        unresolved): the goal is not to recover a unique native conformation but to produce a
        *physically plausible, closed, clash-free* starting structure for MD. Each loop is
        seeded from realistic Ramachandran basins (not a straight chain) and, for each gap, a
        small **ensemble** of independent seeds is closed and the candidate that introduces the
        fewest steric clashes with the rest of the structure is kept -- preserving the one
        virtue of the old straight-chain+steer approach (it added no clashes) while giving a
        far more defensible backbone. The closed coordinates replace the loop atoms in a new
        state PDB; ``connect`` then patches the LINK bond and rebuilds the junction, and a light
        minimization (if enabled) relaxes residual strain.
        """
        import os
        import sys
        import numpy as np
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
        from scipy.spatial.distance import cdist
        from ..psfutil.loop_ccd import loop_atoms_from_pdb, close_one_loop
        from ..util.coord import pdb_replace_coords
        from ..util.progress import PestiferProgress

        self.next_basename('ccd')
        state: StateArtifacts = self.get_current_artifact('state')
        src_pdb = state.pdb.name
        gaps = self.base_molecule.protein_loop_gaps()
        if not gaps:
            logger.warning('ligate method=ccd: has_protein_loops is set but no interior gaps '
                           'were enumerated; nothing to close')
            return 0
        seed = int(ccd_specs.get('seed', 27021972))
        max_iters = int(ccd_specs.get('max_iters', 2000))
        tol = float(ccd_specs.get('tol', 0.1))
        ensemble = max(1, int(ccd_specs.get('ensemble', 10)))
        refine_iters = int(ccd_specs.get('refine', 250))
        guard = bool(ccd_specs.get('guard', False))
        on_clash = str(ccd_specs.get('on_clash', 'warn')).casefold()

        # Serials of EVERY modeled loop copy: while unclosed they hold throwaway guesscoord
        # positions, so all are excluded from every loop's closure environment.
        all_loop_serials = set()
        for gap in gaps:
            _o, _c, gser = loop_atoms_from_pdb(src_pdb, gap['loop_resids'], segname=gap['segname'])
            all_loop_serials.update(gser)
        all_loop_serials = sorted(all_loop_serials)

        def _args(i, extra_env=None):
            g = gaps[i]
            return dict(src_pdb=src_pdb, segname=g['segname'], loop_resids=g['loop_resids'],
                        n_anchor_resid=g['n_anchor_resid'], c_anchor_resid=g['c_anchor_resid'],
                        all_loop_serials=all_loop_serials, seed=seed + i, ensemble=ensemble,
                        refine=refine_iters, guard=guard, max_iters=max_iters, tol=tol,
                        extra_env=extra_env)

        # Loops that do not interfere can be closed at once (steering closes all loops
        # simultaneously too). Close every loop CONCURRENTLY against the resolved structure only;
        # a downstream minimize resolves whatever mutual (non-threaded) clashes result. The one
        # thing minimize cannot undo is one loop's backbone threaded through another's, so a
        # post-pass detects loops whose backbones interpenetrate and RE-CLOSES those sequentially
        # (each avoiding the ones before it). Deterministic: per-loop seed = base + index.
        n = len(gaps)
        results = [None] * n
        nworkers = max(1, min(n, (os.cpu_count() or 2)))

        # Determinate progress bar with a time-remaining estimate. The total work is a known
        # count -- n loops x `ensemble` candidate closures -- so `PestiferProgress(max_value=...)`
        # gives a real ETA rather than a spinner. Each worker pushes one tick per candidate onto
        # `q`; the driver drains it and, when a loop finishes, trues-up progress to that loop's full
        # ensemble share (a loop can break early after fewer than `ensemble` candidates). Only shown
        # on an interactive terminal so batch logs stay clean.
        total = n * max(1, ensemble)
        prog = {'done': 0, 'loops': 0}
        show = sys.stdout.isatty()
        bar = None
        if show:
            bar = PestiferProgress(name='ligate/ccd', color='fuchsia', max_value=total)
            bar.register_update_function(lambda: min(1.0, prog['done'] / total) if total else 1.0)

        def _drain(q):
            while True:
                try:
                    q.get_nowait()
                except Exception:
                    break
                prog['done'] += 1

        if n == 1 or nworkers == 1:
            # serial path: no mid-call draining is possible, so progress advances per loop
            for i in range(n):
                results[i] = close_one_loop(**_args(i))
                prog['loops'] += 1
                prog['done'] = prog['loops'] * ensemble
                if bar:
                    bar.go()
        else:
            with mp.Manager() as mgr:
                q = mgr.Queue()
                with ProcessPoolExecutor(max_workers=nworkers) as ex:
                    futs = {ex.submit(close_one_loop, progress_queue=q, **_args(i)): i
                            for i in range(n)}
                    pending = set(futs)
                    while pending:
                        completed, pending = wait(pending, timeout=0.25,
                                                  return_when=FIRST_COMPLETED)
                        _drain(q)
                        for fut in completed:
                            results[futs[fut]] = fut.result()
                            prog['loops'] += 1
                        # true-up: a finished loop always counts its full ensemble share
                        prog['done'] = max(prog['done'], prog['loops'] * ensemble)
                        if bar:
                            bar.go()
                    _drain(q)
        if bar:
            prog['done'] = total
            bar.finish()
        logger.debug(f'ligate/ccd: closed {n} loop(s) in parallel across {nworkers} worker(s)')

        # Optional threading repair (OFF by default). Closing all loops at once lets copies
        # converging at an assembly axis overlap; verification on the 4zmj trimer shows those
        # overlaps are ordinary clashes that the downstream minimize pushes apart cleanly (they
        # are not topologically interlinked), so no re-closure is needed and the parallelism is
        # kept. Set ccd.thread_ca > 0 to re-close, sequentially, any pair whose backbones come
        # closer than that Ca-Ca distance -- a belt-and-braces option for the (not observed here)
        # case of a genuinely interlinked pair that minimize could not undo.
        THREAD_CA = float(ccd_specs.get('thread_ca', 0.0))
        if THREAD_CA > 0:
            adj = [[] for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    if cdist(results[i]['ca'], results[j]['ca']).min() < THREAD_CA:
                        adj[i].append(j); adj[j].append(i)
            seen = [False] * n
            for i in range(n):
                if seen[i]:
                    continue
                comp, stack = [], [i]
                while stack:
                    u = stack.pop()
                    if seen[u]:
                        continue
                    seen[u] = True; comp.append(u)
                    stack.extend(v for v in adj[u] if not seen[v])
                if len(comp) > 1:
                    comp.sort()
                    logger.info(f'ligate/ccd: {len(comp)} loops interpenetrate after the parallel '
                                f'closure; re-closing them sequentially (ccd.thread_ca={THREAD_CA})')
                    placed = [results[comp[0]]['heavy']]   # keep the first; re-close the rest
                    for idx in comp[1:]:
                        results[idx] = close_one_loop(**_args(idx, extra_env=np.vstack(placed)))
                        placed.append(results[idx]['heavy'])

        # collect coordinates and emit per-loop steric diagnostics
        new_coords = {}
        broken = []
        for i, g in enumerate(gaps):
            res = results[i]
            rep = res['rep']
            tag = f"{g['segname']}:{g['loop_resids'][0]}-{g['loop_resids'][-1]} (len {len(g['loop_resids'])})"
            ndeep = rep['n_deep'] + rep['n_env_deep']
            if rep['min_ca'] < 3.0:     # a crossed backbone -- minimization cannot undo it
                broken.append((tag, rep))
                logger.warning(
                    f'ligate/ccd: loop {tag} is THREADED (min non-adjacent CA {rep["min_ca"]:.2f} A, '
                    f'worst overlap {rep["worst"]:.2f} A) -- a crossed backbone that minimization '
                    f'cannot undo. Increase ccd.ensemble/refine, or provide an externally-modeled '
                    f'loop (see docs/design/loop-modeling.md).')
            elif ndeep:
                logger.info(
                    f'ligate/ccd: loop {tag} closed with {ndeep} deep + '
                    f'{rep["n_soft"] + rep["n_env_soft"]} soft overlaps (worst {rep["worst"]:.2f} A, '
                    f'min CA {rep["min_ca"]:.2f} A) -- not threaded; a downstream minimize relaxes these.')
            elif rep['n_soft'] or rep['n_env_soft']:
                logger.info(f'ligate/ccd: loop {tag} closed with soft contacts only '
                            f'({rep["n_soft"] + rep["n_env_soft"]}; worst {rep["worst"]:.2f} A).')
            for k, serial in enumerate(res['serials']):
                new_coords[serial] = res['closed'][k]

        if broken and on_clash == 'error':
            logger.error(f'ligate/ccd: {len(broken)} loop(s) threaded; aborting (ccd.on_clash=error)')
            return 1

        out_pdb = f'{self.basename}.pdb'
        serial_list = sorted(new_coords)
        coords_arr = np.array([new_coords[s] for s in serial_list])
        row_of_serial = {s: i for i, s in enumerate(serial_list)}
        pdb_replace_coords(src_pdb, out_pdb, coords_arr, row_of_serial)
        self.pdb_to_coor(out_pdb)
        self.register(dict(psf=state.psf,
                           pdb=PDBFileArtifact(self.basename, pytestable=True),
                           coor=NAMDCoorFileArtifact(self.basename)),
                      key='state', artifact_type=StateArtifacts)
        return 0
    
    def write_gaps(self):
        """
        Write the gaps in the base molecule to a data file.
        """
        self.next_basename('gaps')
        mol: Molecule = self.get_current_artifact_data('base_molecule')
        inputfile = f'{self.basename}.inp'
        writer: GenericScripter = self.get_scripter('data')
        writer.newfile(inputfile)
        mol.write_gaps(writer)
        writer.writefile()
        self.register(self.basename, key='measure_distances_input', artifact_type=InputFileArtifact)

    def measure_distances(self, specs):
        """
        Measure the distances between loop termini.
        
        Parameters
        ----------
        specs : dict
            Specifications for the measurement, including the radius of the flexible zone around the receiver.
            This method uses the VMD scripter to create a script that measures the distances between
            the termini of the loops in the base molecule. It generates a data file containing the distances
            and saves the results in a specified output file. The method also updates the state variables
            with the results and the fixed reference structure.
        """
        comment_chars = '#!$'
        self.next_basename('measure')
        vm: VMDScripter = self.get_scripter('vmd')
        vm.newscript(self.basename)
        state: StateArtifacts = self.get_current_artifact('state')
        inputfile: Path = self.get_current_artifact_path('measure_distances_input')
        opdb: str = f'{self.basename}.pdb'
        receiver_flexible_zone_radius: float = specs.get('receiver_flexible_zone_radius', 0.0)
        resultsfile: str = f'{self.basename}.dat'
        vm.addline(f'measure_bonds {state.psf.name} {state.pdb.name} {inputfile.name} {opdb} {resultsfile} {receiver_flexible_zone_radius} ')
        vm.writescript()
        self.register(self.basename, key='measure_distances_tcl', artifact_type=VMDScriptArtifact)
        vm.runscript()
        self.register(self.basename, key='measure_distances_fixedref', artifact_type=PDBFileArtifact)
        self.register(self.basename, key='measure_distances_results', artifact_type=DataFileArtifact)
        self.register(self.basename, key='measure_distances_log', artifact_type=LogFileArtifact)
        with open(resultsfile, 'r') as f:
            datalines = f.read().split('\n')
        self.gaps = []
        for line in datalines:
            if len(line) > 0 and not line[0] in comment_chars:
                data = line.split()
                thisgap = {
                    'segname': data[0],
                    'serial_i': int(data[1]),
                    'serial_j': int(data[2]),
                    'distance': float(data[3])
                }
                self.gaps.append(thisgap)

    def do_steered_md(self, specs):
        """
        Perform steered molecular dynamics to steer the loop termini toward each other.
        """
        self.next_basename('steer')
        writer: NAMDColvarInputScripter = self.get_scripter('namd_colvar')
        writer.newfile(f'{self.basename}-cv.in')
        for i, g in enumerate(self.gaps):
            g['colvars'] = f'GAP{i:02d}'
            writer.declare_distance_cv_atoms(g)
        for i, g in enumerate(self.gaps):
            g['forceConstant'] = specs['force_constant']
            g['targ_distance'] = specs['target_distance']
            g['targ_numsteps'] = specs['nsteps']
            writer.declare_single_harmonic_distance_bias(g)
        writer.writefile()
        self.register(f'{self.basename}-cv', key='steer_colvars', artifact_type=NAMDColvarsConfigArtifact)
        savespecs = self.specs
        self.specs = specs
        result = self.namdrun(extras={
            'fixedatoms': 'on',
            'fixedatomsfile': self.get_current_artifact_path('measure_distances_fixedref'),
            'fixedatomscol': 'O',
            'colvars': 'on',
            'cv configfile': self.get_current_artifact_path('steer_colvars')
            }, single_gpu_only=True)
        # MDTask.namdrun registers the state and all output files
        self.specs = savespecs
        return result

    def connect(self, connect_specs):
        """
        Connect the loop termini using the specified ``LINK`` patch.
        """
        logger.debug(f'Connect specs: {connect_specs} (unused)')
        self.write_connect_patches()
        result = self.connect_gaps()
        return result

    def write_connect_patches(self):
        """
        Write the connect patches to a data file.
        This method generates a data file that contains the connection patches for the loop termini.
        It uses the base molecule to write the connection patches and updates the state variables with the new data file.
        The data file is named based on the current basename, which is generated by the ``next_basename`` method.
        """
        self.next_basename('gap_patches')
        mol: Molecule = self.base_molecule
        datafile = f'{self.basename}.inp'
        writer: GenericScripter = self.get_scripter('data')
        writer.newfile(datafile)
        mol.write_connect_patches(writer)
        writer.writefile()
        self.register(self.basename, key='connect_patches_input', artifact_type=InputFileArtifact)

    def connect_gaps(self):
        """
        Connect the gaps in the loop termini.
        This method uses the psfgen scripter to create a script that connects the gaps in the loop termini
        using the specified patch. It generates a new PSF file and PDB file based on the current state of the base molecule
        and the connection patches defined in the data file. The script is then executed, and if successful,
        the resulting PSF and PDB files are saved in the current state.
        
        Returns
        -------
        int
            The result of the psfgen script execution. A return value of 0 indicates success, while any other value indicates failure.
        """
        self.next_basename('heal')
        pg: PsfgenScripter = self.get_scripter('psfgen')
        pg.newscript(self.basename)
        RM: ResourceManager = self.resource_manager
        CC: CHARMMFFContent = RM.charmmff_content
        CC.copy_charmmfile_local('pestifer.top')
        charmm_topology_files: FileArtifactList = self.get_current_artifact_data('charmmff_topfiles')
        charmm_topology_files.append(CharmmffTopFileArtifact('pestifer', ext='top'))
        pg.addline(f'topology pestifer.top')
        patchfile: Path = self.get_current_artifact_path('connect_patches_input')
        state: StateArtifacts = self.get_current_artifact('state')
        pg.load_project(state.psf.name, state.pdb.name)
        pg.addline(f'source {patchfile.name}')
        pg.writescript(self.basename, guesscoord=True, regenerate=True)
        self.register(self.basename, key='connect_gaps_tcl', artifact_type=PsfgenInputScriptArtifact)
        result = pg.runscript()
        if result == 0:
            self.pdb_to_coor(f'{self.basename}.pdb')
            self.register(dict(
                psf=PSFFileArtifact(self.basename, pytestable=True), 
                pdb=PDBFileArtifact(self.basename, pytestable=True), 
                coor=NAMDCoorFileArtifact(self.basename)), key='state', artifact_type=StateArtifacts)
            self.register(self.basename, key='connect_gaps_log', artifact_type=LogFileArtifact)
        return result
    