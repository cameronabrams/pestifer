
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
        basefile=os.path.join(r['config'],'base.yaml')
        super().__init__(basefile,userfile=userfile)
        self['Resources']=r
        self['Conda']=CondaCheck()
        my_logger(self['Conda'].info(),logger.info,just='<',frame='*',fill='')
        self._set_shortcuts()

    def _set_shortcuts(self):
        self.namd2=self['user']['paths']['namd2']
        self.charmrun=self['user']['paths']['charmrun']
        self.vmd=self['user']['paths']['vmd']
        # self.tcl_path=self['Resources']['tcl']
        self.tcl_root=os.path.join(self['Resources']['root'],'tcl')
        self.tcl_pkg_path=self['Resources']['tcl']['pkg']
        self.tcl_proc_path=self['Resources']['tcl']['proc']
        self.tcl_script_path=self['Resources']['tcl']['scripts']
        self.vmd_startup_script=os.path.join(self.tcl_script_path,'pestifer-vmd.tcl')
        self.charmmff_toppar_path=self['Resources']['charmmff']['toppar']
        self.charmmff_custom_path=self['Resources']['charmmff']['custom']
        self.user_charmmff_toppar_path=''
        if hasattr(self,'user'):
            self.user_charmmff_toppar_path=os.path.join(self['user']['charmff'],'toppar')
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
        self['user']['ambertools']['available']=False
        amberhome=os.environ.get('AMBERHOME',None)
        if amberhome!=None:
            logger.debug(f'AMBERHOME is set to "{amberhome}". It is assumed you have a working AmberTools installation.')
            self['user']['ambertools']['available']=True
            self['user']['ambertools']['local']=True
            return True
        ok_numpy_version=self['user']['ambertools'].get('ok_numpy_version','1.19.5') # for unpatched packmol-memgen
        if not self['Conda'].conda_root:
            logger.debug(f'You do not have any conda environments, and AMBERHOME is not set.')
            logger.debug(f'We recommend you create a new conda environment in which numpy<={ok_numpy_version} and ambertools==21.12')
            return False
        explicit_env=self['user']['ambertools'].get('venv',None)
        if explicit_env!=None:
            logger.debug(f'Checking explicitly named conda env {explicit_env} for numpy<={ok_numpy_version} and ambertools')
            if self['Conda'].env_exists(explicit_env):
                explicit_env_numpy_version=self['Conda'].get_package_version(explicit_env,'numpy')
                explicit_env_ambertools_version=self['Conda'].get_package_version(explicit_env,'ambertools',from_list=True)
                if explicit_env_numpy_version!=None and explicit_env_numpy_version<='1.19.5' and explicit_env_ambertools_version!=None:
                    self['user']['ambertools']['available']=True
                    self['user']['ambertools']['local']=False
                    return True
                else:
                    logger.debug(f'Conda env {explicit_env} is not compatible with ambertools')
            else:
                logger.debug(f'Warning: Explicitly named conda env {explicit_env}  not found.')
        logger.debug(f'Searching all conda envs for an ambertools-compatible env...')
        local_numpy_version=self['Conda'].get_package_version(self['Conda'].active_env,'numpy')
        if local_numpy_version<=ok_numpy_version:
            logger.debug(f'Active env {self["Conda"].active_env} has numpy version {local_numpy_version}')
            ambertools_version=self['Conda'].get_package_version(self['Conda'].active_env,'ambertools',from_list=True)
            if ambertools_version!=None:
                self['user']['ambertools']['available']=True
                self['user']['ambertools']['local']=True
                return True
            else:
                logger.debug(f'Ambertools is not installed in active conda env {self["Conda"].active_env}')
        else:
            for env in self['Conda'].conda_envs:
                if env!=self['Conda'].active_env: # since we already examined it
                    logger.debug(f'Checking non-active conda env {env} for numpy<={ok_numpy_version} and ambertools')
                    remote_numpy_version=self['Conda'].get_package_version(env,'numpy')
                    if remote_numpy_version!=None and remote_numpy_version<=ok_numpy_version:
                        logger.debug(f'Non-active env {env} has numpy version {local_numpy_version}')
                        ambertools_version=self['Conda'].get_package_version(env,'ambertools',from_list=True)
                        logger.debug(f'Non-active env {env} ambertools version "{ambertools_version}"')
                        if ambertools_version!=None:
                            self['user']['ambertools']['available']=True
                            self['user']['ambertools']['local']=False
                            self['user']['ambertools']['venv']=env
                            return True
        logger.debug(f'Ambertools is not available.')
        logger.debug(f'We recommend you create a new conda environment in which numpy<={ok_numpy_version} and ambertools==21.12')
        logger.debug(f'In your future config files, indicate this environment under "ambertools"->"venv".')
        return False