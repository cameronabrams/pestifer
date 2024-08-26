# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import os
import logging
logger=logging.getLogger(__name__)
from collections import UserDict
from . import Resources
from .command import CondaCheck
from .stringthings import my_logger
from ycleptic.yclept import Yclept
from importlib.metadata import version

segtype_of_resname={}
charmm_resname_of_pdb_resname={}
res_321={}
res_123={}

class ResourceManager(UserDict):
    excludes=['__pycache__','_archive','bash'] # exclude from top-level resources
    def __init__(self):
        data={}
        data['root']=os.path.dirname(Resources.__file__)
        roottk=data['root'].split(os.sep)
        rootdepth=len(roottk)
        for root,dirs,files in os.walk(data['root']):
            tk=root.split(os.sep)
            absdepth=len(tk)
            if rootdepth<absdepth<rootdepth+3:
                rn=tk[rootdepth:absdepth]
                bn=tk[-1]
                if not bn in ResourceManager.excludes:
                    myset(data,rn,root)
        super().__init__(data)

def myset(a_dict,keylist,val):
    if len(keylist)==1:
        a_dict[keylist[0]]=val
    else:
        key=keylist.pop(0)
        if not key in a_dict or type(a_dict[key])==type(val):
            a_dict[key]={}
        myset(a_dict[key],keylist,val)

def injest_packmol_memgen_databases():
    logging.getLogger("pmmg_log").setLevel(logging.CRITICAL)
    import packmol_memgen as pmmg
    import packmol_memgen.lib.charmmlipid2amber as cla
    import pandas as pd
    pmmfile=os.path.join(os.path.split(pmmg.__file__)[0],'data','memgen.parm')
    if os.path.exists(pmmfile):
        lipdf=pd.read_csv(pmmfile,sep=r'\s+',header=22)
        logger.debug(f'memgen_parm database has {lipdf.shape[0]} entries')

    cladir=os.path.split(cla.__file__)[0]
    clafile=os.path.join(cladir,'charmmlipid2amber.csv')
    if os.path.exists(clafile):
        cladf=pd.read_csv(clafile,header=1,skiprows=0)
        logger.debug(f'charmmlipid2amber database has {cladf.shape[0]} entries')
    return lipdf,cladf

class Config(Yclept):
    def __init__(self,userfile='',**kwargs):
        vrep=f"""ycleptic v. {version("ycleptic")}
pidibble v. {version("pidibble")}"""
        quiet=kwargs.get('quiet',False)
        if not quiet:
            my_logger(vrep,logger.info,just='<',frame='*',fill='')
        r=ResourceManager()
        logger.debug(f'Resources:')
        for k,v in r.items():
            if type(v)!=dict:
                logger.debug(f'  {k}: {v}')
            else:
                logger.debug(f'  {k}:')
                for kk,vv in v.items():
                    logger.debug(f'      {kk}: {vv}')
        # resolve full pathname of YCleptic base config for this application
        basefile=os.path.join(r['config'],'base.yaml')
        assert os.path.exists(basefile)
        super().__init__(basefile,userfile=userfile)
        self['Resources']=r
        self['Conda']=CondaCheck()
        conda_info=self['Conda'].info()
        cpu_info=self.cpu_info()
        if not quiet:
            my_logger(conda_info,logger.info,just='<',frame='*',fill='')
            my_logger(cpu_info,logger.info,just='<',frame='*',fill='')
        self._set_internal_shortcuts(**kwargs)
        self._set_external_apps(verify_access=(userfile!=''))

    def cpu_info(self):
        self.slurmvars={k:os.environ[k] for k in os.environ if 'SLURM' in k}
        self.local_ncpus=os.cpu_count()
        retstr=''
        if self.slurmvars:
            nnodes=int(self.slurmvars['SLURM_NNODES'])
            ntaskspernode=int(self.slurmvars['SLURM_NTASKS_PER_NODE'])
            ncpus=nnodes*ntaskspernode
            retstr+=f'SLURM: #nodes {nnodes}; total number of cpus {ncpus}'
        else:
            retstr+=f'Local number of CPUs: {self.local_ncpus}'
            ncpus=self.local_ncpus
        self.ncpus=ncpus
        return retstr

    def _set_external_apps(self,verify_access=True):
        self.namd2=self['user']['paths']['namd2']
        self.charmrun=self['user']['paths']['charmrun']
        self.vmd=self['user']['paths']['vmd']
        if verify_access:
            assert os.access(self.charmrun,os.X_OK)
            assert os.access(self.namd2,os.X_OK)
            assert os.access(self.vmd,os.X_OK)

    def _set_internal_shortcuts(self,**kwargs):
        self.progress=len(self.slurmvars)==0
        self.tcl_root=os.path.join(self['Resources']['root'],'tcl')
        assert os.path.exists(self.tcl_root)
        self.tcl_pkg_path=self['Resources']['tcl']['pkg']
        assert os.path.exists(self.tcl_pkg_path)
        self.tcl_script_path=self['Resources']['tcl']['scripts']
        assert os.path.exists(self.tcl_script_path)
        self.vmd_startup_script=os.path.join(self.tcl_root,'vmdrc.tcl')
        assert os.path.exists(self.vmd_startup_script)
        self.charmmff_toppar_path=self['Resources']['charmmff']['toppar']
        assert os.path.exists(self.charmmff_toppar_path)
        self.charmmff_custom_path=self['Resources']['charmmff']['custom']
        assert os.path.exists(self.charmmff_custom_path)
        self.charmmff_pdb_path=self['Resources']['charmmff']['pdb']
        assert os.path.exists(self.charmmff_pdb_path)
        self.user_charmmff_toppar_path=''
        if hasattr(self,'user'):
            self.user_charmmff_toppar_path=os.path.join(self['user']['charmff'],'toppar')
            assert os.path.exists(self.user_charmmff_toppar_path)
        self.namd2_config_defaults=self['user']['namd2']
        self.segtypes=self['user']['psfgen']['segtypes']
        for stn,stspec in self.segtypes.items():
            if stspec and 'rescodes' in stspec:
                self['user']['psfgen']['segtypes'][stn]['invrescodes']={v:k for k,v in stspec['rescodes'].items()}
            if stspec and 'resnames' in stspec:
                initresnames=stspec['resnames']
                for r in initresnames:
                    if len(r)>4 and r.isalnum():
                        rabbrv=r[:4]
                        # logger.debug(r,rabbrv)
                        if not rabbrv in stspec['resnames']:
                            stspec['resnames'].append(rabbrv)
                            # logger.debug(f'Adding abbreviated resname {rabbrv} to resnames for segtype {stn}')
        self.pdb_to_charmm_resnames={}
        for alias in self['user']['psfgen']['aliases']:
            tok=alias.split()
            if tok[0]=='residue':
                self.pdb_to_charmm_resnames[tok[1]]=tok[2]
        self.segtype_of_resname={}
        for st in self['user']['psfgen']['segtypes']:
            res=self['user']['psfgen']['segtypes'][st].get('resnames',[])
            for r in res:
                self.segtype_of_resname[r]=st

        # globals
        charmm_resname_of_pdb_resname.update(self.pdb_to_charmm_resnames)
        for p,c in charmm_resname_of_pdb_resname.items():
            self.segtype_of_resname[c]=self.segtype_of_resname[p]
        segtype_of_resname.update(self.segtype_of_resname)
        res_123.update(self['user']['psfgen']['segtypes']['protein']['invrescodes'])
        res_321.update(self['user']['psfgen']['segtypes']['protein']['rescodes'])

        packmol_memgen_task_specs='packmol_memgen' in [list(x.keys())[0] for x in self['user']['tasks']]
        if packmol_memgen_task_specs or kwargs.get('memgen_test',False):
            # let's supress the irritating message that pops up when packmol_memgen is imported
            atools=self['user']['ambertools']
            atools['memgen_parm_df'],atools['charmmlipid2amber_df']=injest_packmol_memgen_databases()
            if not self.check_ambertools():
                raise Exception(f'You have a "packmol_memgen" task but ambertools is not available.')
        else:
            self['user']['ambertools']['is_unnecessary']=True

    def check_ambertools(self):
        def error_message():
            my_logger(f'Based on your inputs, this instance of {__package__}\nrequires the ambertools package.\nTo use ambertools, {__package__} requires it to be installed via conda.\nWe recommend you create a new conda environment\nand install ambertools and {__package__} in it, and run {__package__} with it activated.',logger.debug)

        if not self['Conda'].conda_root:
            error_message()
            return False
        at_specs=self['user']['ambertools']
        at_specs['available']=False
        active_env=self['Conda'].active_env
        explicit_env_ambertools_version=self['Conda'].get_package_version('ambertools',env=active_env,from_list=True)
        if explicit_env_ambertools_version!=None:
            at_specs['available']=True
            at_specs['local']=True
            logger.debug(f'ambertools v. {explicit_env_ambertools_version}')
            return True
        else:
            error_message()
            return False
        # logger.debug(f'No ambertools package found in active environment {active_env}.')
        # # if active environment can't run ambertools and the explicit env is specified
        # else:
        #     logger.debug(f'Checking explicitly named conda env {explicit_ambertools_venv} for ambertools...')
        #     if self['Conda'].env_exists(explicit_ambertools_venv):
        #         explicit_env_ambertools_version=self['Conda'].get_package_version(explicit_ambertools_venv,'ambertools',from_list=True)
        #         if explicit_env_ambertools_version!=None:
        #             at_specs['available']=True
        #             at_specs['local']=False
        #             logger.debug(f'...ambertools v. {explicit_env_ambertools_version} detected')
        #             return True
        #         else:
        #             logger.debug(f'No ambertools package found in specified environment {explicit_ambertools_venv}.')
        #     else:
        #         logger.debug(f'Warning: No specified conda environment "{explicit_ambertools_venv}" found.')
        # logger.debug(f'No ambertools found after checking active/specified environment, so')
        # logger.debug(f'searching all available non-active conda envs for ambertools...')
        # for env in self['Conda'].conda_envs:
        #     if env!=self['Conda'].active_env: # since we already examined it
        #         logger.debug(f'Checking non-active conda env {env} for ambertools...')
        #         ambertools_version=self['Conda'].get_package_version(env,'ambertools',from_list=True)
        #         if ambertools_version!=None:
        #             logger.debug(f'Non-active env {env} ambertools version "{ambertools_version}"')
        #             at_specs['available']=True
        #             at_specs['local']=False
        #             at_specs['venv']=env
        #             return True
        # logger.debug(f'No ambertools found in any conda environment.')
        # error_message()
        # return False