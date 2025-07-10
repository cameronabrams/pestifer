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
import os
import yaml

from .md import MDTask
from ..core.command import Command

logger=logging.getLogger(__name__)

class TerminateTask(MDTask):
    """
    TerminateTask class for terminating a pestifer build.
    This class inherits from the :class:`MDTask <pestifer.tasks.md.MDTask>` class and is used to prepare the system for termination.
    It handles the copying of state files, writing chain maps, and packaging the system for NAMD runs.
    """
    yaml_header='terminate'
    """
    YAML header for the TerminateTask, used to identify the task in configuration files as part of a ``tasks`` list.
    This header is used to declare TerminateTask objects in YAML task lists.
    It is typically used at the end of a build workflow to finalize the state of the simulation and prepare for termination.
    """
    def do(self):
        """
        Execute the terminate task.
        """
        self.log_message('initiated')
        self.inherit_state()
        self.next_basename()
        self.copy_state(exts=['psf','pdb','coor','xsc','vel'])
        self.write_chainmaps()
        self.write_statefile()
        self.result=self.make_package()
        self.log_message('complete')
        return self.result

    def write_chainmaps(self):
        """
        Write the chain maps to a YAML file.
        This method retrieves the base molecule from the state variables, gets the chain maps, and writes them to a specified YAML file.
        The chain maps are used to map chains in the molecular structure, which is useful for understanding the topology of the system.
        """
        bm=self.statevars.get('base_molecule',None)
        if bm:
            maps=bm.get_chainmaps()
            with open(self.specs['chainmapfile'],'w') as f:
                yaml.dump(maps,f)
            del self.statevars['base_molecule']

    def make_package(self):
        """
        Create a package for the NAMD run.
        This method prepares the necessary files for a NAMD run based on the specifications provided in the task.
        It collects the required files, including the NAMD script, constraint PDB file, and any local parameters.
        The method also handles the generation of a tarball containing all the necessary files for the NAMD run.
        The tarball is named based on the basename specified in the task specifications.
        """
        specs=self.specs.get('package',{})
        if not specs:
            logger.debug('no package specs found')
            return 0
        basename=specs.get('basename','my_system')
        # self.inherit_state()
        self.FC.clear()  # populate a file collector to make the tarball
        logger.debug(f'Packaging for namd using basename {basename}')
        savespecs=self.specs
        self.specs=specs
        result=self.namdrun(script_only=True)
        self.specs=savespecs
        self.FC.append(f'{basename}.namd')
        constraints=specs.get('constraints',{})
        if constraints:
            self.make_constraint_pdb(constraints)
            self.FC.append(self.statevars['consref'])
        local_params=self.statevars.get('charmmff_paramfiles',[])
        for n in local_params:
            self.FC.append(n)
        for ext in ['psf','pdb','coor','xsc','vel']:
            if ext in self.statevars:
                self.FC.append(self.statevars[ext])
        self.FC.tarball(specs["basename"])
        return result
    
        # if specs["topogromacs"]:
        #     logger.debug(f'running topogromacs')
        #     with open(f'{basename}_par.inp','w') as f:
        #         for pf in params['parameters']:
        #             f.write(f'parameters {pf}\n')
        #     vt=self.scripters['vmd']
        #     vt.newscript(f'{basename}_tg')
        #     vt.usescript('tg')
        #     vt.writescript()
        #     psf=self.statevars['psf']
        #     pdb=self.statevars['pdb']
        #     inputname=os.path.splitext(self.statevars['coor'])[0]
        #     vt.runscript(o=basename,pdb=pdb,psf=psf,i=inputname,parinp=f'{basename}_par.inp',ospf=f'{basename}_tg.psf',opdb=f'{basename}_tg.pdb',top=f'{basename}_topogromacs.top',cellfile=f'{basename}_cell.inp')
        #     with open(f'{basename}_cell.inp','r') as f:
        #         box=f.read().split()
        #     boxstr=' '.join(box)
        #     c=Command(f'gmx editconf -f {basename}_tg.pdb -o {basename}_topogromacs.pdb -box {boxstr}')
        #     c.run()
        #     self.FC.append(f'{basename}_topogromacs.pdb')
        #     self.FC.append(f'{basename}_topogromacs.top')
        # self.FC.tarball(specs["basename"])
        # return 0