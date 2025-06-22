# Author: Cameron F. Abrams

import logging
import os

logger=logging.getLogger(__name__)

from ..basetask import BaseTask
from ..colvars import colvar_writer
from ..util.util import is_periodic

class MDTask(BaseTask):
    """ A class for handling all NAMD runs
    
    Attributes
    ----------
    `yaml_header(str)`:
        the name used to declare this task in an input yaml file

    Methods
    -------
    
    `do()`:
        Inherits state and performs namd run based on specs; updates run state

    `copy_charmmpar_local()`:abbr:
        Copies all charmm parameter files from pestifer's resource library to
        the current working directory; returns the list of file names that
        were copied

    `namdrun()`: 
        Generates the NAMD config file based on specs and then executes NAMD

    """
    yaml_header='md'

    def do(self):
        self.log_message('initiated',ensemble=self.specs.get('ensemble',None))
        self.inherit_state()            
        self.result=self.namdrun()
        if self.result==0: 
            self.save_state(exts=['coor','vel','xsc'])
            self.update_statevars('charmmff_paramfiles',self.writers['namd'].parameters)
            if os.path.exists(f'{self.basename}-energy.csv'):
                self.update_statevars('energy-csv',f'{self.basename}-energy.csv',vtype='file')
            else:
                if 'energy-csv' in self.statevars:
                    del self.statevars['energy-csv']
            if os.path.exists(f'{self.basename}-pressureprofile.csv'):
                self.update_statevars('pressureprofile-csv',f'{self.basename}-pressureprofile.csv',vtype='file')
            else:
                if 'pressureprofile-csv' in self.statevars:
                    del self.statevars['pressureprofile-csv']
        else:
            raise Exception(f'md task failed: {self.result}')
        self.log_message('complete',ensemble=self.specs.get('ensemble',None))
        return super().do()

    def namdrun(self,baselabel='',extras={},script_only=False,**kwargs):
        specs=self.specs
        addl_paramfiles=specs.get('addl_paramfiles',[])
        logger.debug(f'md task specs {specs}')
        ensemble=specs['ensemble']
        if ensemble.casefold() not in ['NPT'.casefold(),'NVT'.casefold(),'minimize'.casefold()]:
            raise Exception(f'Error: {ensemble} is not a valid ensemble type. Must be \'NPT\', \'NVT\', or \'minimize\'')
        if not baselabel:
            self.next_basename(ensemble)
        else:
            self.next_basename(baselabel)
        
        params={}
        namd_global_params=self.config['user']['namd']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        coor=self.statevars.get('coor',None)
        vel=self.statevars.get('vel',None)
        xsc=self.statevars.get('xsc',None)
        prior_paramfiles=self.statevars.get('charmmff_paramfiles',[])
        firsttimestep=self.statevars.get('firsttimestep',0)
        if specs.get('firsttimestep',None) is not None:
            firsttimestep=specs['firsttimestep']
        self.statevars['periodic']=is_periodic(xsc)

        temperature=specs['temperature']
        if ensemble.casefold()=='NPT'.casefold():
            pressure=specs['pressure']
        params['tcl']=[]
        params['tcl'].append(f'set temperature {temperature}')

        nsteps=specs['nsteps']
        dcdfreq=specs['dcdfreq']
        xstfreq=specs['xstfreq']

        constraints=specs.get('constraints',{})
        other_params=specs.get('other_parameters',{})
        colvars=specs.get('colvar_specs',{})
        params.update(namd_global_params['generic'])
        params['structure']=psf
        params['coordinates']=pdb
        params['temperature']='$temperature'

        if coor:
            params['bincoordinates']=coor
        if vel:
            params['binvelocities']=vel
            del params['temperature']
        if xsc:
            params['extendedSystem']=xsc
            if ensemble.casefold()=='NPT'.casefold() and xstfreq:
                params['xstfreq']=xstfreq
        
        if self.statevars['periodic']:
            params.update(namd_global_params['solvated'])
        else:
            params.update(namd_global_params['vacuum'])
        
        if ensemble.casefold() in ['NPT'.casefold(),'NVT'.casefold()]:
            params.update(namd_global_params['thermostat'])
            if ensemble.casefold()=='NPT'.casefold():
                if not self.statevars['periodic']:
                    raise Exception(f'Cannot use barostat on a system without PBCs')
                params['tcl'].append(f'set pressure {pressure}')
                params.update(namd_global_params['barostat'])
        params['outputName']=f'{self.basename}'

        if dcdfreq:
            params['dcdfreq']=dcdfreq
        if 'restartfreq' in specs:
            params['restartfreq']=specs['restartfreq']
        if constraints:
            self.make_constraint_pdb(constraints)
            params['constraints']='on'
            params['consref']=self.statevars['consref']
            params['conskfile']=self.statevars['consref']
            params['conskcol']='O'
        if colvars:
            writer=self.writers['data']
            writer.newfile(f'{self.basename}-cv.inp')
            colvar_writer(colvars,writer,pdb=pdb)
            writer.writefile()
            params['colvars']='on'
            params['colvarsconfig']=f'{self.basename}-cv.inp'

        params.update(other_params)
        params.update(extras)
        params['firsttimestep']=firsttimestep
        if ensemble.casefold()=='minimize'.casefold():
            assert specs['minimize']>0,f'Error: you must specify how many minimization cycles'
            params['minimize']=specs['minimize']
        elif ensemble.casefold() in ['NVT'.casefold(),'NPT'.casefold()]:
            assert nsteps>0,f'Error: you must specify how many time steps to run'
            params['run']=nsteps
        
        na=self.writers['namd']
        na.newscript(self.basename,addl_paramfiles=list(set(addl_paramfiles+prior_paramfiles)))
        cpu_override=specs.get('cpu-override',False)
        logger.debug(f'CPU-override is {cpu_override}')
        na.writescript(params,cpu_override=cpu_override)
        
        if not script_only:
            local_execution_only=not self.statevars['periodic']
            single_gpu_only=kwargs.get('single_gpu_only',False) or constraints
            result=na.runscript(single_molecule=(not self.statevars['periodic']),local_execution_only=local_execution_only,single_gpu_only=single_gpu_only,cpu_override=cpu_override)
            if result!=0:
                return -1

            if ensemble.casefold()!='minimize'.casefold():
                self.update_statevars('firsttimestep',firsttimestep+nsteps)
            else:
                self.update_statevars('firsttimestep',firsttimestep+specs['minimize'])
            if specs.get('remove_parfiles',False):
                na.remove_charmm_par(params['parameters'])
        return 0
