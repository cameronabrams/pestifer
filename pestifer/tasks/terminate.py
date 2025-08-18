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
            self.register(YAMLFileArtifact(self.specs['chainmapfile'].replace('.yaml', '')), key='chainmapfile')

    def make_package(self):
        """
        Create a package for a production NAMD run starting from the end of the build.
        """
        package_specs = self.specs.get('package', {})
        state_dir = package_specs.get('state_dir', '.')
        if not package_specs:
            logger.debug('No package specifications provided; packaging will not be performed.')
            return 0
        if 'ensemble' not in package_specs:
            raise ValueError('Ensemble must be specified in package specs for terminate task')
        TarballContents = FileArtifactList()
        self.basename = self.specs.get('basename', 'my_system')
        state: StateArtifacts = self.get_current_artifact('state')
        for ext in ['psf', 'pdb', 'coor', 'xsc', 'vel']:
            fa: FileArtifact = getattr(state, ext, None)
            if fa:
                shutil.copy(fa.name, self.basename + '.' + ext)
                self.register(type(fa)(self.basename))
                TarballContents.append(self.get_current_artifact(ext))
        logger.debug(f'Packaging for namd using basename {self.basename}')
        save_specs = self.specs
        self.specs = package_specs
        result = self.namdrun(script_only=True)
        self.specs = save_specs
        TarballContents.append(self.get_current_artifact('namd'))
        constraints = self.specs.get('constraints', {})
        if constraints:
            self.make_constraint_pdb(constraints, statekey='consref')
            TarballContents.append(self.get_current_artifact('consref'))
        TarballContents.extend(self.get_current_artifact_data('charmmff_parfiles'))
        TarballContents.extend(self.get_current_artifact_data('charmmff_streamfiles'))
        TarballContents.make_tarball(self.basename, arcname_prefix=state_dir)
        return result

    def cleanup(self):

        if not self.specs.get('cleanup', True):
            logger.debug('Cleanup disabled; skipping cleanup step.')
            return 0
        archive_dir = self.specs.get('archive_dir', 'archive')  
        ArtifactContents = FileArtifactList()

        all_my_artifacts = self.pipeline.get_artifact_collection_as_lists()
        file_artifacts = all_my_artifacts['files']
        filelist_artifacts = all_my_artifacts['filelists']

        all_artifact_files = []
        for file_artifact in file_artifacts:
            if not file_artifact.name in all_artifact_files and file_artifact.exists():
                all_artifact_files.append(file_artifact.name)
                ArtifactContents.append(file_artifact)
        for artifact in filelist_artifacts:
            for file_artifact in artifact:
                if not file_artifact.name in all_artifact_files and file_artifact.exists():
                    all_artifact_files.append(file_artifact.name)
                    ArtifactContents.append(file_artifact)

        all_artifact_files.sort()
        # logger.debug(f'All artifact files: {all_artifact_files}')

        non_artifact_files = []
        cwd_files = os.listdir('.')
        # logger.debug(f'Current working directory files: {cwd_files}')
        for f in cwd_files:
            if f not in all_artifact_files:
                non_artifact_files.append(f)

        logger.debug(f'Non-artifact files in current working directory: {non_artifact_files}')

        ArtifactContents.make_tarball('artifacts', remove=True, arcname_prefix=archive_dir)
        return 0