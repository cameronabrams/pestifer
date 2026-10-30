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
import shutil

from pathlib import Path

from .loopclosure import LoopClosureMixin, _RECLOSE_ATTEMPTS  # noqa: F401 (re-exported)
from .mdtask import MDTask

from ..charmmff.charmmffcontent import CHARMMFFContent
from ..core.artifacts import *
from ..core.resourcemanager import ResourceManager
from ..molecule.molecule import Molecule
from ..scripters import GenericScripter, VMDScripter, PsfgenScripter, NAMDColvarInputScripter

logger = logging.getLogger(__name__)

class LigateTask(LoopClosureMixin, MDTask):
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
        If the task is bypassed, it logs a message and returns 0 (success) without performing
        any operations -- a structure with no protein loops is a valid input, not a failure.
        """
        self.base_molecule: Molecule = self.get_current_artifact_data('base_molecule')
        if not self.base_molecule.has_protein_loops:
            # Return success, not None.  A structure with no gaps is a legitimate input to a
            # config that carries a ligate task -- copying an Env example onto a complete
            # structure is the obvious way to get here -- and since a failed task became fatal
            # the bare `return` aborted the whole build with `result: None`.
            self.log_message('bypassed')
            self.result = 0
            return self.result
        closed_by = getattr(self.base_molecule, 'loops_closed_by', None)
        if closed_by:
            # psfgen closed and bonded them already, with this task's ccd settings
            # (LoopClosureMixin); closing again would only re-seed loops the first minimize has
            # already relaxed around.
            logger.info(f'ligate: every interior loop was closed and bonded at {closed_by}, with '
                        f'this task\'s ccd settings; nothing left to do')
            self.result = 0
            return self.result
        method = str(self.specs.get('method', 'ccd')).casefold()
        if method == 'steer':
            logger.debug('Storing sequence gaps.')
            self.write_gaps()
            logger.debug('Measuring gap distances.')
            steering_specs = self.specs.get('steer', {})
            if not steering_specs:
                logger.debug(f'No steering specifications for ligate task; this is a bug; bypassing')
                self.result = 0
                return self.result
            self.measure_distances(steering_specs)
            logger.debug('Steering loop C-termini toward their partner N-termini')
            self.result = self.do_steered_md(self.specs['steer'])
            if self.result != 0:
                return self.result
            self.result = self.finish_loop_closure(self.specs)
        else:
            self.result = self.close_and_connect_loops(self.specs)
        # Relaxation is intentionally left to an explicit downstream `minimize`/`md` task
        # (the standard pattern is to follow `ligate` with one, and `pipeline_contract`
        # REQUIRES one because `ligate` is in `_MOVES_ATOMS`). The clash-filtered ensemble
        # already hands off a closed, clash-free structure ready for that minimization.
        return self.result

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
