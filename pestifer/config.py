
import os
import logging
logger=logging.getLogger(__name__)
from collections import UserDict
from pestifer import PestiferResources
from ycleptic.yclept import Yclept
from ycleptic.yclept import __version__ as __ycleptic_version__

segtype_of_resname={}
charmm_resname_of_pdb_resname={}
res_321={}
res_123={}

class ResourceManager(UserDict):
    excludes=['__pycache__','_archive','bash']
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
        logger.info(f'Ycleptic v {__ycleptic_version__}')
        r=ResourceManager()
        logger.debug(f'Resources {r}')
        basefile=os.path.join(r['config'],'base.yaml')
        super().__init__(basefile,userfile=userfile)
        self['Resources']=r
        self._set_shortcuts()

    def _set_shortcuts(self):
        self.namd2=self['user']['paths']['namd2']
        self.charmrun=self['user']['paths']['charmrun']
        self.vmd=self['user']['paths']['vmd']
        self.tcl_path=self['Resources']['tcl']
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
