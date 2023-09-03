"""

.. module:: config
   :synopsis: Manages all configuration input and parsing, and access to installed resources
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import logging
logger=logging.getLogger(__name__)
import yaml
import glob
from collections import UserDict
# from .util import special_update, replace, is_tool
from pestifer import PestiferResources

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
                    # logger.debug(f'{rn} {root}')
                    myset(data,rn,root)
                    
        # self._walk(data,[x for x in glob.glob(data['root']+'/*') if os.path.isdir(x)])
        super().__init__(data)

def myset(a_dict,keylist,val):
    if len(keylist)==1:
        a_dict[keylist[0]]=val
    else:
        key=keylist.pop(0)
        if not key in a_dict or type(a_dict[key])==type(val):
            a_dict[key]={}
        myset(a_dict[key],keylist,val)

class Config(UserDict):
    def __init__(self,userconfigfile=''):
        data={}
        data['Resources']=ResourceManager()
        logger.debug(f'Resources {data["Resources"]}')
        l=glob.glob(data['Resources']['config']+'/*.yaml')
        logger.debug(f'{l}')
        for fn in l:
            with open(fn,'r') as f:
                bn=os.path.basename(fn)
                bn=os.path.splitext(bn)[0]
                logger.debug(f'Pestifer uses {bn}: {fn}')
                data[bn]=yaml.safe_load(f)
        assert 'base' in data
        assert 'help' in data
        if userconfigfile:
            if os.path.exists(userconfigfile):
                logger.debug(f'Pestifer uses {userconfigfile}: {fn}')
                with open(userconfigfile,'r') as f:
                    data['user']=yaml.safe_load(f)
                    # self._user_defaults()
        super().__init__(data)
        self._set_shortcuts()

    # def _user_defaults(self):
    #     D=self['help']
    #     U=self['user']
    #     dwalk(D,U)

    def _set_shortcuts(self):
        self.tcl_proc_path=self['Resources']['tcl']['proc']
        self.tcl_script_path=self['Resources']['tcl']['scripts']
        self.vmd_startup_script=os.path.join(self.tcl_script_path,'pestifer-vmd.tcl')
        self.charmmff_toppar_path=self['Resources']['charmmff']['toppar']
        self.charmmff_custom_path=self['Resources']['charmmff']['custom']
        self.user_charmmff_toppar_path=''
        if hasattr(self,'user'):
            self.user_charmmff_toppar_path=os.path.join(self['user']['charmff'],'toppar')
        self.namd2_config_defaults=self['base']['namd2']
        for stn,stspec in self['base']['psfgen']['segtypes'].items():
            if stspec and 'rescodes' in stspec:
                self['base']['psfgen']['segtypes'][stn]['invrescodes']={v:k for k,v in stspec['rescodes'].items()}
        self.pdb_to_charmm_resnames={}
        for alias in self['base']['psfgen']['aliases']:
            tok=alias.split()
            if tok[0]=='residue':
                self.pdb_to_charmm_resnames[tok[1]]=tok[2]
    
    def segtype_by_resname(self,resname):
        sts=self['base']['psfgen']['segtypes']
        for s,rl in sts.items():
            if resname in rl['resnames']:
                return s
        return 'other'
    
    def charmmify_resname(self,pdbname):
        return self.pdb_to_charmm_resnames.get(pdbname,pdbname)
    
    def res_321(self,resname):
        return self['base']['psfgen']['segtypes']['protein']['rescodes'][resname]

    def res_123(self,rescode):
        return self['base']['psfgen']['segtypes']['protein']['invrescodes'][rescode]
    
def dwalk(D,I):
    tld=[x['name'] for x in D['directives']]
    # logger.debug(f'dwalk along {tld} for {I}')
    for d in tld:
        tidx=tld.index(d)
        dx=D['directives'][tidx]
        # logger.debug(f' d {d}')
        typ=dx['type']
        # logger.debug(f'- {d} typ {typ} I {I}')
        if not d in I:
            # logger.debug(f' -> not found {d}')
            if typ in ['str','int','float', 'bool']:
                if 'default' in dx:
                    I[d]=dx['default']
                    # logger.debug(f' ->-> default {d} {I[d]}')
                elif 'required' in dx:
                    if dx['required']:
                        raise Exception(f'Directive {d} requires a value in {I}')
            elif typ=='dict':
                if 'required' in dx:
                    if not dx['required']:
                        continue
                I[d]={}
                dwalk(dx,I[d])
            elif typ=='list':
                if 'required' in dx:
                    if not dx['required']:
                        continue
                I[d]=[] # and do nothing else
        else:
            if typ=='str' and 'choices' in dx:
                assert I[d] in dx['choices'],f'Directive {d} must be one of {", ".join(dx["choices"])}'
            if typ=='dict':
                dwalk(dx,I[d])
            elif typ=='list':
                lwalk(dx,I[d])

def lwalk(D,L):
    if 'directives' in D:
        tld=[x['name'] for x in D['directives']]
        # logger.debug(f'lwalk on {tld}')
        for item in L:
            # check this item against its directive
            itemname=list(item.keys())[0]
            # logger.debug(f' - item {item}')
            assert itemname in tld
            tidx=tld.index(itemname)
            dx=D['directives'][tidx]
            typ=dx['type']
            if typ in ['str','int','float']:
                pass # don't really expect a list of scalars?
            elif typ=='dict':
                if not item[itemname]:
                    item[itemname]={}
                dwalk(dx,item[itemname])
    

