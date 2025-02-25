# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import os
import yaml

from .md import MDTask
from ..command import Command

logger=logging.getLogger(__name__)

class TerminateTask(MDTask):  #need to inherit for namdrun() method
    yaml_header='terminate'
    
    def do(self):
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
        bm=self.statevars.get('base_molecule',None)
        if bm:
            maps=bm.get_chainmaps()
            with open(self.specs['chainmapfile'],'w') as f:
                yaml.dump(maps,f)
            del self.statevars['base_molecule']

    def make_package(self):
        specs=self.specs.get('package',{})
        # logger.debug(f'make_package specs {specs}')
        if not specs:
            return 0
        self.inherit_state()
        self.FC.clear()  # populate a file collector to make the tarball
        logger.debug(f'Packaging for namd using basename {self.basename}')
        savespecs=self.specs
        self.specs=specs
        params={}
        result=self.namdrun(script_only=True)
        self.specs=savespecs
        self.FC.append(f'{self.basename}.namd')
        constraints=specs.get('constraints',{})
        if constraints:
            self.make_constraint_pdb(constraints)
            self.FC.append(self.statevars['consref'])
        local_params=self.local_parameter_files
        for n in local_params:
            self.FC.append(n)
        for ext in ['psf','pdb','coor','xsc','vel']:
            if ext in self.statevars:
                self.FC.append(self.statevars[ext])

        if specs["topogromacs"]:
            logger.debug(f'running topogromacs')
            with open(f'{self.basename}_par.inp','w') as f:
                for pf in params['parameters']:
                    f.write(f'parameters {pf}\n')
            vt=self.writers['vmd']
            vt.newscript(f'{self.basename}_tg')
            vt.usescript('tg')
            vt.writescript()
            psf=self.statevars['psf']
            pdb=self.statevars['pdb']
            inputname=os.path.splitext(self.statevars['coor'])[0]
            vt.runscript(o=self.basename,pdb=pdb,psf=psf,i=inputname,parinp=f'{self.basename}_par.inp',ospf=f'{self.basename}_tg.psf',opdb=f'{self.basename}_tg.pdb',top=f'{self.basename}_topogromacs.top',cellfile=f'{self.basename}_cell.inp')
            with open(f'{self.basename}_cell.inp','r') as f:
                box=f.read().split()
            boxstr=' '.join(box)
            c=Command(f'gmx editconf -f {self.basename}_tg.pdb -o {self.basename}_topogromacs.pdb -box {boxstr}')
            c.run()
            self.FC.append(f'{self.basename}_topogromacs.pdb')
            self.FC.append(f'{self.basename}_topogromacs.top')
        self.FC.tarball(specs["basename"])
        return 0