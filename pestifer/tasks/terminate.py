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

from .md import MDTask
from ..core.artifacts import ArtifactFileList

logger=logging.getLogger(__name__)

class TerminateTask(MDTask):
    """
    TerminateTask class for terminating a pestifer build.
    This class inherits from the :class:`MDTask <pestifer.tasks.md.MDTask>` class and is used to prepare the system for termination.
    It handles the copying of state files, writing chain maps, and packaging the system for NAMD runs.
    """
    _yaml_header='terminate'
    """
    YAML header for the TerminateTask, used to identify the task in configuration files as part of a ``tasks`` list.
    This header is used to declare TerminateTask objects in YAML task lists.
    It is typically used at the end of a build workflow to finalize the state of the simulation and prepare for termination.
    """
    def do(self):
        """
        Execute the terminate task.
        """
        self.next_basename()
        coor=self.get_current_artifact('coor')
        if not coor:
            self.pdb_to_coor()
        self.write_chainmaps()
        self.result=self.make_package()
        return self.result

    def write_chainmaps(self):
        """
        Write the chain maps to a YAML file.
        This method retrieves the base molecule from the state variables, gets the chain maps, and writes them to a specified YAML file.
        The chain maps are used to map chains in the molecular structure, which is useful for understanding the topology of the system.
        """
        bm=self.get_current_artifact_value('base_molecule')
        if bm:
            maps=bm.get_chainmaps()
            with open(self.specs['chainmapfile'],'w') as f:
                yaml.dump(maps,f)

    def make_package(self):
        """
        Create a package for the NAMD run.
        This method prepares the necessary files for a NAMD run based on the specifications provided in the task.
        It collects the required files, including the NAMD script, constraint PDB file, and any local parameters.
        The method also handles the generation of a tarball containing all the necessary files for the NAMD run.
        The tarball is named based on the basename specified in the task specifications.
        """
        TarballContents = ArtifactFileList()
        specs=self.specs.get('package',{})
        if not specs:
            logger.debug('no package specs found')
            return 0
        self.basename=specs.get('basename','my_system')
        for ext in ['psf','pdb','coor','xsc','vel']:
            fa=self.get_current_artifact(ext)
            if fa:
                shutil.copy(fa.path, self.basename + '.' + ext)
                self.register_current_artifact(type(fa)(self.basename))
                TarballContents.append(self.get_current_artifact(ext))
        # self.inherit_state()
        logger.debug(f'Packaging for namd using basename {self.basename}')
        savespecs=self.specs
        self.specs=specs
        result=self.namdrun(script_only=True)
        self.specs=savespecs
        TarballContents.append(self.get_current_artifact('namd'))
        constraints=specs.get('constraints',{})
        if constraints:
            self.make_constraint_pdb(constraints,statekey='consref')
            TarballContents.append(self.get_current_artifact('consref'))
        TarballContents.extend(self.get_current_artifact_value('charmmff_parfiles'))
        TarballContents.extend(self.get_current_artifact_value('charmmff_streamfiles'))
        TarballContents.make_tarball(self.basename)
        return result
