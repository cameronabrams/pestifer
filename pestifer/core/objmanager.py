# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A dictionary object for containing all objs from a task specification.  Objs are objects that handle various modifications
to a molecular structure, such as cleavages, translations, rotations, and so on.  These are all defined in the :mod:`pestifer.objs` subpackage.
This module defines the :class:`ObjManager` class, which is responsible for managing these objects and their organization.
It provides methods for ingesting objects, filtering them, and managing their categories.
The :class:`ObjManager` inherits from :class:`collections.UserDict`, allowing it to behave like a dictionary
while providing additional functionality specific to managing molecular objects.
It also includes methods for counting objects, retiring categories, and expelling objects based on specific criteria.
"""
from __future__ import annotations

import logging

from collections import UserDict
from typing import Callable
from pestifer.objs.resid import ResID

from .baseobj import BaseObj, BaseObjList

from ..objs.cfusion import Cfusion, CfusionList
from ..objs.cleavagesite import CleavageSite, CleavageSiteList
from ..objs.crot import Crot, CrotList
from ..objs.deletion import Deletion, DeletionList
from ..objs.insertion import Insertion, InsertionList
from ..objs.graft import Graft, GraftList
from ..objs.link import Link, LinkList
from ..objs.mutation import Mutation, MutationList
from ..objs.orient import Orient, OrientList
from ..objs.patch import Patch, PatchList
from ..objs.rottrans import RotTrans, RotTransList
from ..objs.seqadv import Seqadv, SeqadvList
from ..objs.ssbond import SSBond, SSBondList
from ..objs.ssbonddelete import SSBondDelete, SSBondDeleteList
from ..objs.substitution import Substitution, SubstitutionList
from ..objs.ter import Ter, TerList

from ..molecule.residue import Residue, ResidueList

logger = logging.getLogger(__name__)

_ObjCats = {'seq', 'topol', 'coord', 'generic'}
"""
Object categories used in the ObjManager.
This set defines the categories of objects that can be managed by the `ObjManager`.
Each category corresponds to a specific type of modification or information related to molecular structures.

- ``seq``: Sequence-related objects
- ``topol``: Topology-related objects
- ``coord``: Coordinate-related objects
- ``generic``: Generic objects that do not fit into the other categories

"""
class ObjManager(UserDict[str, UserDict[str, BaseObjList]]):
    """
    A class for initializing and collecting all objs into 
    a single, organized object.  Outermost key is the object category
    and the innermost is the object type name.

    Parameters
    ----------
    input_specs : dict, optional
        A dictionary of object specifications to be ingested into the ObjManager.
        If not provided, an empty dictionary is used.

    Attributes
    ----------
    used : dict
        A dictionary that stores retired objects, allowing for retrieval of previously used objects.
    """

    _obj_classes = (
        (Cfusion, CfusionList),
        (CleavageSite, CleavageSiteList),
        (Crot, CrotList),
        (Deletion, DeletionList),
        (Insertion, InsertionList),
        (Graft, GraftList),
        (Link, LinkList),
        (Mutation, MutationList),
        (Orient, OrientList),
        (Patch, PatchList),
        (RotTrans, RotTransList),
        (Seqadv, SeqadvList),
        (SSBond, SSBondList),
        (SSBondDelete, SSBondDeleteList),
        (Substitution, SubstitutionList),
        (Ter, TerList)
    )
    _obj_classes_byYAML = {cls[0]._yaml_header: cls[0] for cls in _obj_classes}
    _objlist_classes_byYAML = {cls[0]._yaml_header: cls[1] for cls in _obj_classes}    

    def __init__(self, input_specs={}):
        """ 
        Initializes the ObjManager with a dictionary of object specifications.
        This method sets up the object classes and their corresponding list classes,
        and prepares the ObjManager to ingest objects based on the provided specifications.

        Parameters
        ----------  
        input_specs : dict
            dictionary of obj shortcode specifications
        """
        self.data = {}
        super().__init__({})
        self.used = {}
        self.ingest(input_specs)

    def filter_copy(self, condition: Callable[[BaseObj], bool], objnames: list[str] = []) -> ObjManager:
        """
        Returns a copy of the ObjManager with only the objects that match the given fields.
        
        Parameters
        ----------
        objnames : list
          list of object names to filter by; if empty, all objects are included
        fields : dict
          dictionary of fields to filter by; only objects that match all fields are included

        Returns
        -------
        ObjManager
            a new ObjManager containing only the objects that match the given fields. If no object names are provided, an empty ObjManager is returned.
        """
        result = ObjManager()
        # self.counts()
        for objcat, catdict in self.data.items():
            for header, objlist in catdict.items():
                if len(objnames) > 0 and header not in objnames:
                    continue
                objlist = objlist.filter(condition)
                if len(objlist) > 0:
                    if not objcat in result:
                        result[objcat] = {}
                    result[objcat][header] = objlist
        return result

    def ingest(self, input_obj, overwrite=False):
        """
        Ingest an object or list of objects into the ObjManager.

        Parameters
        ----------
        input_obj : object, list, or dict
            an object of type Obj, a list of objects, or a dictionary of objects to be ingested
        overwrite : bool
            if True, overwrite existing objects with the same header; if False, append to existing objects
        """
        if input_obj is None:
            return
        if type(input_obj) in [tup[0] for tup in self._obj_classes]:  # an obj class
            self._ingest_obj(input_obj)
        elif type(input_obj) in [tup[1] for tup in self._obj_classes]:  # a list class
            if len(input_obj) > 0:
                logger.debug(f'Ingesting {type(input_obj)} of length {len(input_obj)}')
            self._ingest_objlist(input_obj, overwrite=overwrite)
        elif type(input_obj) == dict:
            self._ingest_objdict(input_obj, overwrite=overwrite)
        elif type(input_obj) == list and len(input_obj) == 0:  # a blank call
            pass
        else:
            raise TypeError(f'Cannot ingest object of type {type(input_obj)} into objmanager')

    def _ingest_obj(self, a_obj):
        Cls = type(a_obj)
        LCls = self._objlist_classes_byYAML[Cls._yaml_header]
        objcat = Cls._objcat
        assert objcat in _ObjCats, f'Object category {objcat} is not recognized'
        header = Cls._yaml_header
        if not objcat in self:
            self[objcat] = {}
        if not header in self[objcat]:
            self[objcat][header] = LCls()
        self[objcat][header].append(a_obj)
        logger.debug(f'Ingested {str(a_obj)} into {objcat} {header}')

    def _ingest_objlist(self, a_objlist, overwrite=False):
        LCls = type(a_objlist)
        Cls = LCls._item_type
        objcat = Cls._objcat
        header = Cls._yaml_header
        logger.debug(f'Ingesting {LCls.__name__} with header {header} into {objcat}')
        if not objcat in self:
            self[objcat] = {}
        if overwrite or not header in self[objcat]:
            self[objcat][header] = LCls()
        self[objcat][header].extend(a_objlist)

    def _ingest_objdict(self, objdict: dict, overwrite=False):
        # does nothing if objdict is empty, but let's just be sure
        for objYAMLname, objlist in objdict.items():
            Cls = self._obj_classes_byYAML.get(objYAMLname, None)
            LCls = self._objlist_classes_byYAML.get(objYAMLname, None)
            objcat = Cls._objcat
            header = Cls._yaml_header
            if not objcat in self:
                self[objcat] = {}
            if overwrite or not header in self[objcat]:
                self[objcat][header] = LCls([])
            for entry in objlist:
                # logger.debug(f'entry {entry.__class__.__name__} objcat {objcat} header {header}')
                if not isinstance(entry, Cls):
                    # logger.debug(f'Converting {entry} to {Cls.__name__}')
                    self[objcat][header].append(Cls(entry))
                else:
                    self[objcat][header].append(entry)

    def retire(self, objcat):
        """
        Retire an object category from the ObjManager.
        
        Parameters
        ----------
        objcat : str
          the object category to retire; if it exists, it will be removed from the ObjManager and stored in the used dictionary
        """
        if objcat in self:
            self.used[objcat] = self[objcat].copy()
            del self[objcat]

    def expel(self, expelled_residues: ResidueList) -> 'ObjManager':
        """
        Expel all objects that match the given residues from the ObjManager.

        Parameters
        ----------
        expelled_residues : list
          list of residues to expel; each residue should have a resseqnum and insertion attribute
        
        Returns
        -------
        ObjManager
            a new ObjManager containing the expelled objects; the original ObjManager is modified to remove these objects
        """
        ex = ObjManager()
        for objcat, catdict in self.items():
            for header, objlist in catdict.items():
                if len(objlist) == 0:
                    continue
                LCls = self._objlist_classes_byYAML[header]
                expelled_objs = LCls([])
                for obj in objlist:
                    for r in expelled_residues:
                        assert isinstance(r, Residue), f'Expected Residue, got {type(r)}'
                        assert isinstance(obj, LCls._item_type), f'Expected {LCls._item_type}, got {type(obj)}'
                        assert isinstance(r.resid, ResID)
                        matches = obj.wildmatch(resid=r.resid)
                        if matches:
                            logger.debug(f'{obj.__class__.__name__} {obj} matches residue {repr(r)}')
                            if not obj in expelled_objs:
                                expelled_objs.append(obj)
                        # else:
                        #     logger.debug(f'{obj.__class__.__name__} {obj} does not match residue {r.describe()}')
                if len(expelled_objs) > 0:
                    ex.ingest(expelled_objs)
                    for obj in expelled_objs:
                        # logger.debug(f'Expelling {obj.__class__.__name__} {obj} from {objcat} {header}')
                        objlist.remove(obj)
        return ex

    def counts_by_header(self) -> dict:
        """
        Counts the number of objects in each category and header.
        This method logs the counts of objects in each category and header.
        It iterates through the object classes and their corresponding categories and headers,
        counting the number of objects in each header within each category.
        """
        counts = {} 
        for objcat, catdict in self.items():
            for header, objlist in catdict.items():
                counts[header] = len(objlist)
        for YAMLname, Cls in self._obj_classes_byYAML.items():
            if YAMLname not in counts:
                counts[YAMLname] = 0
        return counts
