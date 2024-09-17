# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import os
import logging
from collections import UserDict
from importlib.metadata import version
from ycleptic.yclept import Yclept

from . import Resources
from .stringthings import my_logger

logger=logging.getLogger(__name__)
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
        cpu_info=self.cpu_info()
        if not quiet:
            my_logger(cpu_info,logger.info,just='<',frame='*',fill='')
        self._set_internal_shortcuts(**kwargs)
        self._set_external_apps(verify_access=(userfile!=''))

    def cpu_info(self):
        self.slurmvars={k:os.environ[k] for k in os.environ if 'SLURM' in k}
        self.local_ncpus=os.cpu_count()
        retstr=''
        if self.slurmvars and 'SLURM_NNODES' in self.slurmvars and 'SLURM_NTASKS_PER_NODE' in self.slurmvars:
            # we are in a batch execution
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
        self.packmol=self['user']['paths']['packmol']
        if verify_access:
            assert os.access(self.charmrun,os.X_OK)
            assert os.access(self.namd2,os.X_OK)
            assert os.access(self.vmd,os.X_OK)
            assert os.access(self.packmol,os.X_OK)

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

