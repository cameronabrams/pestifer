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
import textwrap
from collections import UserDict
from pestifer import PestiferResources

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
                logger.debug(f'Pestifer reads user config {userconfigfile}')
                with open(userconfigfile,'r') as f:
                    data['user']=yaml.safe_load(f)
                    # self._user_defaults()
        else:
            data['user']={}
        super().__init__(data)
        self._user_defaults()
        self._set_shortcuts()

    def _user_defaults(self):
        D=self['help']
        U=self['user']
        dwalk(D,U)

    def dump_user(self,filename='complete-user.yaml'):
        with open(filename,'w') as f:
            f.write('# Pestifer -- Cameron F. Abrams -- cfa22@drexel.edu\n')
            f.write('# Dump of complete user config file\n')
            yaml.dump(self['user'],f)

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
        self.namd2_config_defaults=self['base']['namd2']
        self.segtypes=self['base']['psfgen']['segtypes']
        for stn,stspec in self.segtypes.items():
            if stspec and 'rescodes' in stspec:
                self['base']['psfgen']['segtypes'][stn]['invrescodes']={v:k for k,v in stspec['rescodes'].items()}
        self.pdb_to_charmm_resnames={}
        for alias in self['base']['psfgen']['aliases']:
            tok=alias.split()
            if tok[0]=='residue':
                self.pdb_to_charmm_resnames[tok[1]]=tok[2]
        self.segtype_of_resname={}
        for st in self['base']['psfgen']['segtypes']:
            res=self['base']['psfgen']['segtypes'][st].get('resnames',[])
            for r in res:
                self.segtype_of_resname[r]=st

        # globals
        charmm_resname_of_pdb_resname.update(self.pdb_to_charmm_resnames)
        for p,c in charmm_resname_of_pdb_resname.items():
            self.segtype_of_resname[c]=self.segtype_of_resname[p]
        segtype_of_resname.update(self.segtype_of_resname)
        res_123.update(self['base']['psfgen']['segtypes']['protein']['invrescodes'])
        res_321.update(self['base']['psfgen']['segtypes']['protein']['rescodes'])


    def make_default_specs(self,*args):
        holder={}
        make_def(self['help']['directives'],holder,*args)
        return holder

def make_def(L,H,*args):
    if len(args)==1:
        name=args[0]
        try:
            item_idx=[x["name"] for x in L].index(name)
        except:
            raise ValueError(f'{name} is not a recognized directive')
        item=L[item_idx]
        for d in item.get("directives",[]):
            if "default" in d:
                H[d["name"]]=d["default"]
            else:
                H[d["name"]]=None
        if not "directives" in item:
            if "default" in item:
                H[item["name"]]=item["default"]
            else:
                H[item["name"]]=None
    elif len(args)>1:
        arglist=list(args)
        nextarg=arglist.pop(0)
        args=tuple(arglist)
        try:
            item_idx=[x["name"] for x in L].index(nextarg)
        except:
            raise ValueError(f'{nextarg} is not a recognized directive')
        item=L[item_idx]
        make_def(item["directives"],H,*args)

def userhelp(L,logf,*args,end=''):
    if len(args)==0:
        logf(f'    Help available for {", ".join([dspec["name"] for dspec in L])}{end}')
    elif len(args)==1:
        name=args[0]
        try:
            item_idx=[x["name"] for x in L].index(name)
        except:
            raise ValueError(f'{name} is not a recognized directive')
        item=L[item_idx]
        logf(f'{item["name"]}:{end}')
        logf(f'    {textwrap.fill(item["text"],subsequent_indent="      ")}{end}')
        logf(f'    type: {item["type"]}{end}')
        if "default" in item:
            logf(f'    default: {item["default"]}{end}')
        if "choices" in item:
            logf(f'    allowed values: {", ".join(item["choices"])}{end}')
        if item.get("required",False):
            logf(f'    A value is required.{end}')
        if "directives" in item:
            userhelp(item["directives"],logf,end=end)
    else:
        arglist=list(args)
        nextarg=arglist.pop(0)
        args=tuple(arglist)
        try:
            item_idx=[x["name"] for x in L].index(nextarg)
        except:
            raise ValueError(f'{nextarg} is not a recognized directive')
        item=L[item_idx]
        logf(f'{nextarg}->{end}')
        userhelp(item['directives'],logf,*args,end=end)

def dwalk(D,I):
    """Process the user's config-dict I by walking recursively through it 
       along with the default config-specification dict D
       
       I is the dict yaml-read from the user input
       D is thd config-specification dict yaml-read from the package resources
    """
    # get the name of each config directive at this level in this block
    tld=[x['name'] for x in D['directives']]
    if I==None:
        raise ValueError(f'Null dictionary found; expected a dict with key(s) {tld} under \'{D["name"]}\'.')
    ud=list(I.keys())
    for u in ud:
        if not u in tld:
            raise ValueError(f'Directive \'{u}\' invalid; expecting one of {tld} under \'{D["name"]}\'.')
    # logger.debug(f'dwalk along {tld} for {I}')
    # for each directive name
    for d in tld:
        # get its index in the list of directive names
        tidx=tld.index(d)
        # get its dictionary
        dx=D['directives'][tidx]
        # logger.debug(f' d {d}')
        # get its type
        typ=dx['type']
        # logger.debug(f'- {d} typ {typ} I {I}')
        # if this directive name does not already have a key in the result
        if not d in I:
            # logger.debug(f' -> not found {d}')
            # if it is a scalar
            if typ in ['str','int','float', 'bool']:
                # if it has a default, set it
                if 'default' in dx:
                    I[d]=dx['default']
                    # logger.debug(f' ->-> default {d} {I[d]}')
                # if it is flagged as required, die since it is not in the read-in
                elif 'required' in dx:
                    if dx['required']:
                        raise Exception(f'Directive \'{d}\' of \'{D["name"]}\' requires a value.')
            # if it is a dict
            elif typ=='dict':
                # if it is explicitly tagged as not required, do nothing
                if 'required' in dx:
                    if not dx['required']:
                        continue
                # whether required or not, set it as empty and continue the walk,
                # which will set defaults for all descendants
                I[d]={}
                dwalk(dx,I[d])
            elif typ=='list':
                if 'required' in dx:
                    if not dx['required']:
                        continue
                I[d]=[] # and do nothing else
        # this directive does appear in I
        else:
            if typ=='str' and 'choices' in dx:
                assert I[d] in dx['choices'],f'Directive \'{d}\' of \'{dx["name"]}\' must be one of {", ".join(dx["choices"])}'
            if typ=='dict':
                # process descendants
                dwalk(dx,I[d])
            elif typ=='list':
                # process list-item children
                lwalk(dx,I[d])

def lwalk(D,L):
    if 'directives' in D:
        tld=[x['name'] for x in D['directives']]
        # logger.debug(f'lwalk on {tld}')
        for item in L:
            # check this item against its directive
            itemname=list(item.keys())[0]
            # logger.debug(f' - item {item}')
            if not itemname in tld:
                raise ValueError(f'Element \'{itemname}\' of list \'{D["name"]}\' is not valid; expected one of {tld}')
            tidx=tld.index(itemname)
            dx=D['directives'][tidx]
            typ=dx['type']
            if typ in ['str','int','float']:
                # because a list directive indicates an ordered sequence of tasks and we expect each
                # task to be a dictionary specifying the task and not a single scalar value,
                # we will ignore this one
                logger.debug(f'Warning: Scalar list-element-directive \'{dx}\' in \'{dx["name"]}\' ignored.')
            elif typ=='dict':
                if not item[itemname]:
                    item[itemname]={}
                dwalk(dx,item[itemname])
            else:
                logger.debug(f'Warning: List-element-directive \'{itemname}\' in \'{dx["name"]}\' ignored.')


