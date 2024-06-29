
import os
import logging
logger=logging.getLogger(__name__)
from collections import UserDict
from pestifer import PestiferResources
from .command import CondaCheck
from .stringthings import my_logger
from ycleptic.yclept import Yclept
from ycleptic.yclept import __version__ as __ycleptic_version__
from pidibble.pdbparse import __version__ as __pidibble_version__

segtype_of_resname={}
charmm_resname_of_pdb_resname={}
res_321={}
res_123={}

class ResourceManager(UserDict):
    excludes=['__pycache__','_archive','bash'] # exclude from top-level resources
    def __init__(self):
        data={}
        data['root']=os.path.dirname(PestiferResources.__file__)
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

class Config(Yclept):
    def __init__(self,userfile=''):
        logger.info(f'ycleptic {__ycleptic_version__}')
        logger.info(f'pidibble {__pidibble_version__}')
        r=ResourceManager()
        logger.debug(f'Resources {r}')
        # resolve full pathname of YCleptic base config for this application
        basefile=os.path.join(r['config'],'base.yaml')
        assert os.path.exists(basefile)
        super().__init__(basefile,userfile=userfile)

        self['Resources']=r
        self['Conda']=CondaCheck()
        my_logger(self['Conda'].info(),logger.info,just='<',frame='*',fill='')
        self._set_shortcuts()

    def _set_shortcuts(self):

        self.namd2=self['user']['paths']['namd2']
        assert os.access(self.namd2,os.X_OK)
        self.charmrun=self['user']['paths']['charmrun']
        assert os.access(self.charmrun,os.X_OK)
        self.vmd=self['user']['paths']['vmd']
        assert os.access(self.vmd,os.X_OK)

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
                            logger.debug(f'Adding abbreviated resname {rabbrv} to resnames for segtype {stn}')
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
        if packmol_memgen_task_specs:
            if not self.check_ambertools():
                raise Exception(f'You have a "packmol_memgen" task but ambertools is not available.')
        else:
            self['user']['ambertools']['is_unnecessary']=True

    def check_ambertools(self):
        def error_message():
            my_logger(f'Based on your inputs, this instance of {__package__}\nrequires the ambertools package,\nbut no conda environments were detected.\nTo use ambertools, {__package__} requires it to be installed via conda.\nWe recommend you create a new conda environment\nand install ambertools and {__package__} in it, and run {__package__} with it activated.\nIn your future {__package__} config files, indicate this environment under "ambertools"->"venv".',logger.debug)

        if not self['Conda'].conda_root:
            error_message()
            return False
        at_specs=self['user']['ambertools']
        at_specs['available']=False
        explicit_ambertools_venv=at_specs.get('venv',None)
        # patch_ambertools=at_specs.get('patch_packmol_memgen_pdbremix',False)
        # if no explicit env for running ambertools is specified, 
        # check the active environment first
        if explicit_ambertools_venv==None:
            active_env=self['Conda'].active_env
            explicit_env_ambertools_version=self['Conda'].get_package_version(active_env,'ambertools',from_list=True)
            if explicit_env_ambertools_version!=None:
                at_specs['available']=True
                at_specs['local']=True
                logger.debug(f'ambertools v. {explicit_env_ambertools_version} detected')
                # if patch_ambertools:
                #     if self['Conda'].patch_packmol_memgen_pdbremix(active_env):
                #         logger.warning(f'I just now patched packmol_memgen pdbremix in the active environment {active_env}!!')
                return True
            logger.debug(f'No ambertools package found in active environment {active_env}.')
        # if active environment can't run ambertools and the explicit env is specified
        else:
            logger.debug(f'Checking explicitly named conda env {explicit_ambertools_venv} for ambertools...')
            if self['Conda'].env_exists(explicit_ambertools_venv):
                explicit_env_ambertools_version=self['Conda'].get_package_version(explicit_ambertools_venv,'ambertools',from_list=True)
                if explicit_env_ambertools_version!=None:
                    at_specs['available']=True
                    at_specs['local']=False
                    logger.debug(f'...ambertools v. {explicit_env_ambertools_version} detected')
                    # if patch_ambertools:
                    #     self['Conda'].patch_packmol_memgen_pdbremix(explicit_ambertools_venv)
                    #     logger.warning(f'I just now patched packmol_memgen pdbremix in the specified environment {explicit_ambertools_venv}!!')
                    return True
                else:
                    logger.debug(f'No ambertools package found in specified environment {explicit_ambertools_venv}.')
            else:
                logger.debug(f'Warning: No specified conda environment "{explicit_ambertools_venv}" found.')
        logger.debug(f'No ambertools found after checking active/specified environment, so')
        logger.debug(f'searching all available non-active conda envs for ambertools...')
        for env in self['Conda'].conda_envs:
            if env!=self['Conda'].active_env: # since we already examined it
                logger.debug(f'Checking non-active conda env {env} for ambertools...')
                ambertools_version=self['Conda'].get_package_version(env,'ambertools',from_list=True)
                if ambertools_version!=None:
                    logger.debug(f'Non-active env {env} ambertools version "{ambertools_version}"')
                    at_specs['available']=True
                    at_specs['local']=False
                    at_specs['venv']=env
                    # if patch_ambertools:
                    #     self['Conda'].patch_packmol_memgen_pdbremix(env)
                    #     logger.warning(f'I just now patched packmol_memgen pdbremix in the environment {env}!!')
                    return True
        logger.debug(f'No ambertools found in any conda environment.')
        error_message()
        return False