# Author: Cameron F. Abrams <cfa22@drexel.edu>.
""" A dictionary object for containing all mods from a task specification
"""
from .util import inspect_classes
from .mods import ModTypes
from collections import UserDict
import logging
logger=logging.getLogger(__name__)

class ModManager(UserDict):
    """A class for initializing and collecting all mods into 
       a single, organized object 
       
    """

    mod_classes,modlist_classes=inspect_classes('pestifer.mods','List')
    
    def __init__(self,input_specs={}):
        """ 
        Arguments
        ---------
        input_specs: dict
          dictionary of mod shortcode specifications
        """
        self.used={}
        super().__init__({})
        self.injest(input_specs)

    def filter_copy(self,modnames=[],**fields):
        result=ModManager()
        # logger.debug(f'ModManager filter_copy modnames {modnames} fields {fields}')
        self.counts()
        for name,Cls in self.mod_classes.items():
            modtype=Cls.modtype
            header=Cls.yaml_header
            # logger.debug(f'passing {modtype} {header}...')
            # logger.debug(f'1-{header in modnames} 2-{modtype in self}')
            # if modtype in self and header in self[modtype]:
            #     logger.debug(f'   3-{header in self[modtype]}')
            if header in modnames and modtype in self and header in self[modtype]:
                modlist=self[modtype][header].filter(**fields)
                # logger.debug(f'{len(modlist)} filtered from {len(self[modtype][header])} mods')
                if len(modlist)>0:
                    if not modtype in result:
                        result[modtype]={}
                    result[modtype][header]=modlist
        return result

    def injest(self,input_obj,overwrite=False):
        if type(input_obj) in self.mod_classes.values():
            self._injest_mod(input_obj)
        elif type(input_obj) in self.modlist_classes.values():
            return self._injest_modlist(input_obj,overwrite=overwrite)
        elif type(input_obj)==dict:
            self._injest_moddict(input_obj,overwrite=overwrite)
        elif type(input_obj)==list and len(input_obj)==0: # a blank call
            pass
        else:
            raise TypeError(f'Cannot injest object of type {type(input_obj)} into modmanager')
    
    def _injest_mod(self,a_mod):
        Cls=type(a_mod)
        modclassident=[x for x,y in self.mod_classes.items() if y==Cls][0]
        LCls=self.modlist_classes.get(f'{modclassident}List',list)
        modtype=Cls.modtype
        assert modtype in ModTypes,f'Modtype {modtype} is not recognized'
        header=Cls.yaml_header
        if not modtype in self:
            self[modtype]={}
        if not header in self[modtype]:
            self[modtype][header]=LCls([])
        self[modtype][header].append(a_mod)

    def _injest_modlist(self,a_modlist,overwrite=False):
        if len(a_modlist)==0: # can handle an empty list...
            return a_modlist  # ...by returning it
        LCls=type(a_modlist)
        Cls=type(a_modlist[0])
        modtype=Cls.modtype
        header=Cls.yaml_header
        if not modtype in self:
            self[modtype]={}
        if overwrite or not header in self[modtype]:
            self[modtype][header]=LCls([])
        self[modtype][header].extend(a_modlist)
        return self[modtype][header]
        
    def _injest_moddict(self,moddict,overwrite=False):
        # does nothing if moddict is empty, but let's just be sure
        if len(moddict)==0:
            return
        for name,Cls in self.mod_classes.items():
            modtype=Cls.modtype
            header=Cls.yaml_header
            if header in moddict:
                LCls=self.modlist_classes.get(f'{name}List',list)
                if not modtype in self:
                    self[modtype]={}
                if overwrite or not header in self[modtype]:
                    self[modtype][header]=LCls([])
                for entry in moddict[header]:
                    assert type(entry) in [str,dict],f'Error: expected mod specification of type str or dict, got {type(entry)}'
                    self[modtype][header].append(Cls(entry))

    def retire(self,modtype):
        if modtype in self:
            self.used[modtype]=self[modtype].copy()
            del self[modtype]
    
    def expel(self,expelled_residues):
        ex=ModManager()
        for name,Cls in self.mod_classes.items():
            LCls=self.modlist_classes.get(f'{name}List',list)
            modtype=Cls.modtype
            if modtype in self:
                exl=LCls([])
                for mod in self[modtype]:
                    for r in expelled_residues:
                        matches=mod.wildmatch(resseqnum=r.resseqnum,insertion=r.insertion)
                        if matches:
                            exl.append(mod)
                for mod in exl:
                    self[modtype].remove(mod)
                ex.injest(exl)
        return ex
    
    def counts(self):
        for name,Cls in self.mod_classes.items():
            modtype=Cls.modtype
            header=Cls.yaml_header
            this_count=0
            if modtype in self:
                if header in self[modtype]:
                    this_count=len(self[modtype][header])
            logger.debug(f'{header} {this_count}')