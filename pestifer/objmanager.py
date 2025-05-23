# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" A dictionary object for containing all objs from a task specification
"""
from .util.util import inspect_package_dir
from collections import UserDict
import logging
logger=logging.getLogger(__name__)
import os

from . import objs

ObjCats=['seq','topol','coord','generic']

class ObjManager(UserDict):
    """A class for initializing and collecting all objs into 
       a single, organized object 
       
    """
    obj_classes,objlist_classes=inspect_package_dir(os.path.dirname(objs.__file__),key='List')
    
    def __init__(self,input_specs={}):
        """ 
        Arguments
        ---------
        input_specs: dict
          dictionary of obj shortcode specifications
        """
        self.used={}
        super().__init__({})
        self.injest(input_specs)

    def filter_copy(self,objnames=[],**fields):
        result=ObjManager()
        # logger.debug(f'ObjManager filter_copy objnames {objnames} fields {fields}')
        self.counts()
        for name,Cls in self.obj_classes.items():
            objcat=Cls.objcat
            header=Cls.yaml_header
            # logger.debug(f'passing {objcat} {header}...')
            # logger.debug(f'1-{header in objnames} 2-{objcat in self}')
            # if objcat in self and header in self[objcat]:
            #     logger.debug(f'   3-{header in self[objcat]}')
            if header in objnames and objcat in self and header in self[objcat]:
                objlist=self[objcat][header].filter(**fields)
                # logger.debug(f'{len(objlist)} filtered from {len(self[objcat][header])} objs')
                if len(objlist)>0:
                    if not objcat in result:
                        result[objcat]={}
                    result[objcat][header]=objlist
        return result

    def injest(self,input_obj,overwrite=False):
        if type(input_obj) in self.obj_classes.values():
            self._injest_obj(input_obj)
        elif type(input_obj) in self.objlist_classes.values():
            return self._injest_objlist(input_obj,overwrite=overwrite)
        elif type(input_obj)==dict:
            self._injest_objdict(input_obj,overwrite=overwrite)
        elif type(input_obj)==list and len(input_obj)==0: # a blank call
            pass
        else:
            raise TypeError(f'Cannot injest object of type {type(input_obj)} into objmanager')
    
    def _injest_obj(self,a_obj):
        Cls=type(a_obj)
        objclassident=[x for x,y in self.obj_classes.items() if y==Cls][0]
        LCls=self.objlist_classes.get(f'{objclassident}List',list)
        objcat=Cls.objcat
        assert objcat in ObjCats,f'Object category {objcat} is not recognized'
        header=Cls.yaml_header
        if not objcat in self:
            self[objcat]={}
        if not header in self[objcat]:
            self[objcat][header]=LCls([])
        self[objcat][header].append(a_obj)

    def _injest_objlist(self,a_objlist,overwrite=False):
        if len(a_objlist)==0: # can handle an empty list...
            return a_objlist  # ...by returning it
        LCls=type(a_objlist)
        Cls=type(a_objlist[0])
        objcat=Cls.objcat
        header=Cls.yaml_header
        if not objcat in self:
            self[objcat]={}
        if overwrite or not header in self[objcat]:
            self[objcat][header]=LCls([])
        self[objcat][header].extend(a_objlist)
        return self[objcat][header]
        
    def _injest_objdict(self,objdict,overwrite=False):
        # does nothing if objdict is empty, but let's just be sure
        logger.debug(f'objdict {objdict}')
        if len(objdict)==0:
            return
        for name,Cls in self.obj_classes.items():
            objcat=Cls.objcat
            header=Cls.yaml_header
            if header in objdict:
                LCls=self.objlist_classes.get(f'{name}List',list)
                if not objcat in self:
                    self[objcat]={}
                if overwrite or not header in self[objcat]:
                    self[objcat][header]=LCls([])
                for entry in objdict[header]:
                    assert type(entry) in [str,dict],f'Error: expected obj specification of type str or dict, got {type(entry)}'
                    logger.debug(f'entry {entry} objcat {objcat} header {header}')
                    self[objcat][header].append(Cls(entry))

    def retire(self,objcat):
        if objcat in self:
            self.used[objcat]=self[objcat].copy()
            del self[objcat]
    
    def expel(self,expelled_residues):
        ex=ObjManager()
        for name,Cls in self.obj_classes.items():
            LCls=self.objlist_classes.get(f'{name}List',list)
            objcat=Cls.objcat
            if objcat in self:
                exl=LCls([])
                for obj in self[objcat]:
                    for r in expelled_residues:
                        matches=obj.wildmatch(resseqnum=r.resseqnum,insertion=r.insertion)
                        if matches:
                            exl.append(obj)
                for obj in exl:
                    self[objcat].remove(obj)
                ex.injest(exl)
        return ex
    
    def counts(self):
        for name,Cls in self.obj_classes.items():
            objcat=Cls.objcat
            header=Cls.yaml_header
            this_count=0
            if objcat in self:
                if header in self[objcat]:
                    this_count=len(self[objcat][header])
            logger.debug(f'{header} {this_count}')