# Author: Cameron F. Abrams
"""
Definition of the :class:`MDTask` class for handling molecular dynamics (MD) simulations using NAMD.
This class is a descendant of the :class:`BaseTask <pestifer.core.basetask.BaseTask>` class and is used to run NAMD simulations,
manage the state of the simulation, and handle various aspects of the MD task such as ensembles, constraints, and colvars.
It manages the setup and execution of NAMD runs, including reading and writing necessary files,
updating artifacts, and handling the results of the simulation.

Usage is described in the :ref:`subs_runtasks_md` documentation.
"""
import logging
import os

logger=logging.getLogger(__name__)

from ..core.basetask import VMDTask
from ..core.pipeline import PDBFile, VELFile, XSCFile, LogFile, COORFile, DCDFile, CSVFile, NAMDConfigFile, CVFile, XSTFile
from ..util.namdcolvars import colvar_writer
from ..util.util import is_periodic

class MDTask(VMDTask):
    """ 
    A class for handling all NAMD runs
    """
    yaml_header='md'
    """
    YAML header for the MDTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """
    def do(self):
        """
        Execute the MD task.
        """
        self.log_message('initiated',ensemble=self.specs.get('ensemble',None))
        self.result=self.namdrun()
        if self.result==0:
            self.log_message('complete',ensemble=self.specs.get('ensemble',None))
        else:
            raise RuntimeError(f'md task failed: {self.result}')

    def namdrun(self,baselabel='',extras={},script_only=False,**kwargs):
        """
        Run a NAMD simulation based on the specifications provided in the task.
        This method prepares the necessary parameters, writes the NAMD script, and executes it.
        It handles different ensembles (NPT, NVT, minimize), sets up constraints and colvars, and manages the simulation state.
        If `script_only` is True, it only writes the script without executing it.
        
        Parameters
        ----------
        baselabel : str, optional
            The base label for the NAMD run. If not provided, it will be generated based on the ensemble type.
        extras : dict, optional
            Additional parameters to be included in the NAMD script.
        script_only : bool, optional
            If True, only the script will be written without executing it.
        **kwargs : dict, optional
            Additional keyword arguments that may include execution options such as `single_gpu_only`.
        """
        specs=self.specs
        addl_paramfiles=specs.get('addl_paramfiles',[])
        logger.debug(f'md task specs {specs}')
        ensemble=specs['ensemble']
        if ensemble.casefold() not in ['NPAT'.casefold(),'NPT'.casefold(),'NVT'.casefold(),'minimize'.casefold()]:
            raise Exception(f'Error: {ensemble} is not a valid ensemble type. Must be \'NPAT\', \'NPT\', \'NVT\', or \'minimize\'')
        if not baselabel:
            self.next_basename(ensemble)
        else:
            self.next_basename(baselabel)
        
        params={}
        namd_global_params=self.config['user']['namd']
        psf=self.get_current_artifact_value('psf')
        pdb=self.get_current_artifact_value('pdb')
        coor=self.get_current_artifact_value('coor')
        vel=self.get_current_artifact_value('vel')
        xsc=self.get_current_artifact_value('xsc')
        prior_paramfiles=[x.name for x in self.get_current_artifact_value('charmmff_paramfiles')]
        firsttimestep=self.get_current_artifact_value('firsttimestep')
        if not firsttimestep:
            firsttimestep=0
        if specs.get('firsttimestep',None) is not None:
            firsttimestep=specs['firsttimestep']
        if xsc is not None:
            self.register_current_artifact('periodic',is_periodic(xsc.path))

        temperature=specs['temperature']
        if ensemble.casefold()=='NPT'.casefold() or ensemble.casefold()=='NPAT'.casefold():
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
        params['structure']=psf.path
        params['coordinates']=pdb.path
        params['temperature']='$temperature'

        if coor:
            params['bincoordinates']=coor.path
        if vel:
            params['binvelocities']=vel.path
            del params['temperature']
        if xsc:
            params['extendedSystem']=xsc.path
            if (ensemble.casefold()=='NPT'.casefold() or ensemble.casefold()=='NPAT'.casefold()) and xstfreq:
                params['xstfreq']=xstfreq

        if self.get_current_artifact_value('periodic'):
            params.update(namd_global_params['solvated'])
        else:
            params.update(namd_global_params['vacuum'])
        
        if ensemble.casefold() in ['NPT'.casefold(),'NVT'.casefold(), 'NPAT'.casefold()]:
            params.update(namd_global_params['thermostat'])
            if ensemble.casefold()=='NPT'.casefold():
                if not self.get_current_artifact_value('periodic'):
                    raise Exception(f'Cannot use barostat on a system without PBCs')
                params['tcl'].append(f'set pressure {pressure}')
                params.update(namd_global_params['barostat'])
            elif ensemble.casefold()=='NPAT'.casefold():
                params['tcl'].append(f'set pressure {pressure}')
                params.update(namd_global_params['membrane'])
        params['outputName']=f'{self.basename}'

        if dcdfreq:
            params['dcdfreq']=dcdfreq
        if 'restartfreq' in specs:
            params['restartfreq']=specs['restartfreq']
        if constraints:
            self.make_constraint_pdb(constraints)
            consref=self.get_current_artifact_value('consref')
            if not consref:
                raise RuntimeError(f'No constraint reference PDB file found for task {self.taskname}')
            params['constraints']='on'
            params['consref']=consref.path
            params['conskfile']=consref.path
            params['conskcol']='O'
        if colvars:
            writer=self.scripters['data']
            writer.newfile(f'{self.basename}-cv.inp')
            colvar_writer(colvars,writer,pdb=pdb)
            writer.writefile()
            self.register_current_artifact('colvarsconfig',CVFile(f'{self.basename}-cv.inp'))
            params['colvars']='on'
            params['colvarsconfig']=f'{self.basename}-cv.inp'

        params.update(other_params)
        params.update(extras)
        params['firsttimestep']=firsttimestep
        if ensemble.casefold()=='minimize'.casefold():
            assert specs['minimize']>0,f'Error: you must specify how many minimization cycles'
            params['minimize']=specs['minimize']
        elif ensemble.casefold() in ['NVT'.casefold(),'NPT'.casefold(),'NPAT'.casefold()]:
            assert nsteps>0,f'Error: you must specify how many time steps to run'
            params['run']=nsteps
        
        na=self.scripters['namd']
        na.newscript(self.basename,addl_paramfiles=list(set(addl_paramfiles+prior_paramfiles)))
        self.register_current_artifact('charmmff_params',na.parameters)
        cpu_override=specs.get('cpu-override',False)
        logger.debug(f'CPU-override is {cpu_override}')
        na.writescript(params,cpu_override=cpu_override)
        self.register_current_artifact('namdscript',NAMDConfigFile(f'{self.basename}.namd'))
        if not script_only:
            local_execution_only=not self.get_current_artifact_value('periodic')
            single_gpu_only=kwargs.get('single_gpu_only',False) or constraints
            result=na.runscript(single_molecule=local_execution_only,local_execution_only=local_execution_only,single_gpu_only=single_gpu_only,cpu_override=cpu_override)
            namd_output_exts={e:t for e,t in zip(['coor','vel','xsc','log','dcd','xst'],[COORFile,VELFile,XSCFile,LogFile,DCDFile,XSTFile])}
            for ext,artifact_type in namd_output_exts.items():
                if os.path.exists(f'{self.basename}.{ext}'):
                    self.register_current_artifact(ext,artifact_type(f'{self.basename}.{ext}'))
            self.coor_to_pdb()
            self.register_current_artifact('pdb',PDBFile(f'{self.basename}.pdb'))
            if hasattr(na.logparser,'dataframes'):
                for key in na.logparser.dataframes:
                    self.register_current_artifact(f'{key}-csv',CSVFile(f'{self.basename}-{key}.csv'))
            if result!=0:
                return -1

            if ensemble.casefold()!='minimize'.casefold():
                self.register_current_artifact('firsttimestep',firsttimestep+nsteps)
            else:
                self.register_current_artifact('firsttimestep',firsttimestep+specs['minimize'])

        return 0
