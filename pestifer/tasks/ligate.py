# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import os

from .md import MDTask
from ..colvars import declare_distance_cv_atoms, declare_single_harmonic_distance_bias

logger=logging.getLogger(__name__)

class LigateTask(MDTask):
    yaml_header='ligate'

    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        self.base_molecule=self.statevars['base_molecule']
        if not self.base_molecule.has_loops(min_loop_length=self.statevars['min_loop_length']):
            self.log_message('bypassed')
            return
        logger.debug('Storing sequence gaps.')
        self.write_gaps()
        logger.debug('Measuring gap distances.')
        self.measure_distances(self.specs['steer'])
        logger.debug('Steering loop C-termini toward their partner N-termini')
        self.result=self.do_steered_md(self.specs['steer'])
        if self.result!=0:
            return self.result
        self.save_state(exts=['coor','vel'])
        logger.debug('Connecting loop C-termini to their partner N-termini')
        self.result=self.connect()
        if self.result!=0:
            return self.result
        self.save_state(exts=['psf','pdb'])
        self.log_message('complete')
        return self.result
    
    def write_gaps(self):
        self.next_basename('gaps')
        mol=self.base_molecule
        datafile=f'{self.basename}.inp'
        writer=self.writers['data']
        writer.newfile(datafile)
        mol.write_gaps(writer)
        writer.writefile()
        self.update_statevars('data',datafile,vtype='file')

    def measure_distances(self,specs):
        comment_chars='#!$'
        self.next_basename('measure')
        vm=self.writers['vmd']
        vm.newscript(self.basename)
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        datafile=self.statevars['data']
        opdb=f'{self.basename}.pdb'
        receiver_flexible_zone_radius=specs.get('receiver_flexible_zone_radius',0.0)
        resultsfile=f'{self.basename}.dat'
        vm.addline(f'measure_bonds {psf} {pdb} {datafile} {opdb} {resultsfile} {receiver_flexible_zone_radius} ')
        vm.writescript()
        vm.runscript()
        self.update_statevars('fixedref',f'{self.basename}.pdb',vtype='file')
        self.update_statevars('results',resultsfile,vtype='file')
        with open(resultsfile,'r') as f:
            datalines=f.read().split('\n')
        self.gaps=[]
        for line in datalines:
            if len(line)>0 and not line[0] in comment_chars:
                data=line.split()
                thisgap={
                    'segname':data[0],
                    'serial_i':int(data[1]),
                    'serial_j':int(data[2]),
                    'distance':float(data[3])
                }
                self.gaps.append(thisgap)

    def do_steered_md(self,specs):
        self.next_basename('steer')
        writer=self.writers['data']
        writer.newfile(f'{self.basename}-cv.inp')
        for i,g in enumerate(self.gaps):
            g['colvars']=f'GAP{i:02d}'
            declare_distance_cv_atoms(g,writer)
        for i,g in enumerate(self.gaps):
            g['forceConstant']=specs['force_constant']
            g['targ_distance']=specs['target_distance']
            g['targ_numsteps']=specs['nsteps']
            declare_single_harmonic_distance_bias(g,writer)
        writer.writefile()
        savespecs=self.specs
        self.specs=specs
        result=self.namdrun(extras={        
            'fixedatoms':'on',
            'fixedatomsfile':self.statevars['fixedref'],
            'fixedatomscol': 'O',
            'colvars': 'on',
            'colvarsconfig': f'{self.basename}-cv.inp'
            },single_gpu_only=True)
        self.specs=savespecs
        return result

    def connect(self):
        self.write_connect_patches()
        result=self.connect_gaps()
        return result

    def write_connect_patches(self):
        self.next_basename('gap_patches')
        mol=self.base_molecule
        datafile=f'{self.basename}.inp'
        writer=self.writers['data']
        writer.newfile(datafile)
        mol.write_connect_patches(writer)
        writer.writefile()
        self.update_statevars('data',datafile,vtype='file')

    def connect_gaps(self):
        self.next_basename('heal')
        pg=self.writers['psfgen']
        pg.newscript(self.basename)
        # pg.topo_aliases()
        topfile=os.path.join(self.config.charmmff_custom_path,'unter.top')
        pg.addline(f'topology {topfile}')
        # pg.usescript('loop_closure')
        patchfile=self.statevars['data']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        pg.load_project(psf,pdb)
        pg.addline(f'source {patchfile}')
        pg.writescript(self.basename,guesscoord=True,regenerate=True)
        result=pg.runscript()
        if result==0:
            self.save_state(exts=['psf','pdb'])
        return result
    