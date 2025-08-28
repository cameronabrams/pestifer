# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`TerminateTask` class for terminating a pestifer build.
This task is a descendant of the :class:`MDTask <pestifer.tasks.md.MDTask>` class and is used to prepare the system for termination.
It handles the copying of state files, writing chain maps, and packaging the system for NAMD runs.
The task also manages the state of the simulation, including the base molecule and various file extensions such as PSF, PDB, COOR, XSC, and VEL.
The task is designed to be used in a workflow where the simulation needs to be gracefully terminated and packaged for further analysis or continuation.
It ensures that all necessary files are collected and organized, making it easy to resume or analyze the simulation later.
"""
import logging
import shutil
import yaml
import os

from .mdtask import MDTask
from ..core.artifacts import *
from ..molecule.molecule import Molecule

logger = logging.getLogger(__name__)

class TerminateTask(MDTask):
    """
    TerminateTask class for terminating a pestifer build.
    This class inherits from the :class:`MDTask <pestifer.tasks.md.MDTask>` class and is used to prepare the system for termination.
    It handles the copying of state files, writing chain maps, and packaging the system for NAMD runs.
    """
    _yaml_header = 'terminate'
    """
    YAML header for the TerminateTask, used to identify the task in configuration files as part of a ``tasks`` list.
    This header is used to declare TerminateTask objects in YAML task lists.
    It is typically used at the end of a build workflow to finalize the state of the simulation and prepare for termination.
    """

    def do(self) -> int:
        self.next_basename()
        if 'chainmapfile' in self.specs:
            self.write_chainmaps()
        self.result = self.make_package() | self.cleanup()
        return self.result

    def write_chainmaps(self):
        """
        Write the chain maps to a YAML file.
        This method retrieves the base molecule from the state variables, gets the chain maps, and writes them to a specified YAML file.
        The chain maps are used to map chains in the molecular structure, which is useful for understanding the topology of the system.
        """
        bm: Molecule = self.get_current_artifact_data('base_molecule')
        if bm:
            maps = bm.get_chainmaps()
            with open(self.specs['chainmapfile'], 'w') as f:
                yaml.dump(maps, f)
            self.register(self.specs['chainmapfile'].replace('.yaml', ''), key='chainmapfile', artifact_type=YAMLFileArtifact)

    def make_package(self):
        """
        Create a package for a production NAMD run starting from the end of the build.
        """
        package_specs = self.specs.get('package', {})
        md_specs = package_specs.get('md', {})
        state_dir = package_specs.get('state_dir', '.')
        if not package_specs:
            logger.debug('No package specifications provided; packaging will not be performed.')
            return 0
        TarballContents = FileArtifactList()
        self.basename = self.specs.get('basename', 'my_system')
        state: StateArtifacts = self.get_current_artifact('state')
        for ext in ['psf', 'pdb', 'coor', 'xsc', 'vel']:
            fa: FileArtifact = getattr(state, ext, None)
            if fa and fa.exists():
                shutil.copy(fa.name, self.basename + '.' + ext)
                fa.data = self.basename
                assert fa.exists()
                TarballContents.append(fa)
        result = 0
        if md_specs:
            logger.debug(f'Packaging for namd using basename {self.basename}')
            save_specs = self.specs
            self.specs = md_specs
            self.specs['basename'] = package_specs.get('basename', self.basename)
            result = self.namdrun(script_only=True)
            self.specs = save_specs
            TarballContents.append(self.get_current_artifact('namd'))
            constraints = self.specs.get('constraints', {})
            if constraints:
                self.make_constraint_pdb(constraints, statekey='consref')
                TarballContents.append(self.get_current_artifact('consref'))
        else:
            logger.debug(f'No NAMD configuration is included in the package.')
        TarballContents.extend(self.get_current_artifact_data('charmmff_parfiles'))
        TarballContents.extend(self.get_current_artifact_data('charmmff_streamfiles'))
        TarballContents.make_tarball(self.basename, arcname_prefix=state_dir, unique=True)
        return result

    def cleanup(self):

        if not self.specs.get('cleanup', True):
            logger.debug('Cleanup disabled; skipping cleanup step.')
            return 0
        archive_dir = self.specs.get('archive_dir', 'archive')  

        ArtifactContents = self.pipeline.get_all_file_artifacts()

        ArtifactContents.sort(key=lambda x: x.name)
        # logger.debug(f'All artifact files: {all_artifact_files}')

        non_artifact_files = []
        cwd_files = os.listdir('.')
        # logger.debug(f'Current working directory files: {cwd_files}')
        for f in cwd_files:
            if f not in ArtifactContents:
                non_artifact_files.append(f)

        logger.debug(f'Non-artifact files in current working directory: {non_artifact_files}')

        ArtifactContents.make_tarball('artifacts', remove=True, arcname_prefix=archive_dir, unique=True)
        return 0