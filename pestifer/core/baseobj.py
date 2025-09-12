# Author: Cameron F. Abrams, <cfa22@drexel.edu>
# with extensive contributions from ChatGPT 4o and 5
"""
:class:`pydantic.BaseModel` objects with attribute controls

The :class:`BaseObj` class provides a framework for creating objects with controlled attributes.  It is an abstract class that requires subclasses to implement the `_adapt()` static method that returns a dict that sets attribute values.  It also provides a number of utility methods for creating objects from various input types, validating attributes, and comparing objects.

Any subclass that defines Field objects will automatically have their attributes validated according to the rules specified in the BaseObj class.  Subclasses can also define their own ClassVar attributes that are outside the control logic of BaseObj.

Two types of Fields are supported:

1. Required fields, which must be present in the input data and are validated by Pydantic.
2. Optional fields, which are not required but can be included in the input data, and have default values defined.

Pairs of optional fields can be defined as mutually exclusive, meaning that if one field is present, the other cannot be. This is useful for cases where only one of a set of related fields should be specified.

Any field can also be constrained to a set of allowed values, and dependencies can be defined such that if one field is present, another must also be present.

Any field can also be assigned a list of dependent fields depending on its value, which must also be present if the field value is set. This allows for complex relationships between fields to be defined.

Any field can also be declared as ignored when comparing objects, meaning that it will not be included in the comparison logic.

"""
from __future__ import annotations
from argparse import Namespace
from pydantic import BaseModel, model_validator, ConfigDict, model_serializer
from typing import ClassVar, Any, Iterable, Iterator, Self, TypeVar, Generic, get_args, get_origin, Callable
import hashlib
import operator
import yaml
import logging
from collections import UserList

logger = logging.getLogger(__name__)
from abc import abstractmethod, ABCMeta

class BaseObj(BaseModel):
    _required_fields:    ClassVar[set[str]]                       = set()
    _optional_fields:    ClassVar[set[str]]                       = set()
    _mutually_exclusive: ClassVar[set[frozenset[str]]]            = set()
    _attr_choices:       ClassVar[dict[str, set[Any]]]            = {}
    _attr_dependencies:  ClassVar[dict[str, dict[Any, set[Any]]]] = {}
    _ignore_fields:      ClassVar[set[str]]                       = set()

    _objcat: ClassVar[str] = ''

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        extra='forbid',
        frozen=False,
    )
    """
    Configuration for :class:`pydantic.BaseModel`.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(**self.__class__._adapt(*args, **kwargs))

    @property
    def objcat(self):
        return self._objcat

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """Default: supports dict, Namespace, or plain kwargs."""
        if args:
            arg = args[0]
            if isinstance(arg, dict):
                return arg
            elif isinstance(arg, Namespace):
                return vars(arg)
        return kwargs

    def __repr__(self) -> str:
        cls_name = self.__class__.__name__
        parts = []

        # Track fields required due to dependencies
        dependency_required_fields = set()

        for field, conditions in self._attr_dependencies.items():
            if field in self.__class__.__pydantic_fields__:
                field_value = getattr(self, field, None)
                if '*' in conditions:
                    dependency_required_fields |= conditions['*']
                if field_value in conditions:
                    dependency_required_fields |= conditions[field_value]

        # Build output
        for name, field in self.__class__.__pydantic_fields__.items():
            value = getattr(self, name, None)
            if name in self._required_fields:
                parts.append(f"{name}={value!r}")
            elif name in self._optional_fields:
                if value is not None or name in dependency_required_fields:
                    parts.append(f"{name}={value!r}")

        return f"{cls_name}({', '.join(parts)})"

    def copy(self, deep=False) -> Self:
        """
        Uses :meth:`pydantic.BaseModel.model_copy` to create a copy of the instance.
        If deep is True, performs a deep copy of nested objects.
        """
        return self.model_copy(deep=deep)

    @model_validator(mode="before")
    @classmethod
    def _validate_schema(cls, values: dict[str, Any]) -> dict[str, Any]:
        declared_fields = cls._required_fields | cls._optional_fields
        received_fields = set(values)

        # Ensure all fields are explicitly declared
        unexpected = received_fields - declared_fields
        if unexpected:
            raise ValueError(
                f"Fields not declared in _required_fields or _optional_fields: {sorted(unexpected)}"
            )

        # Required field check
        missing = [f for f in cls._required_fields if f not in values]
        if missing:
            raise ValueError(f"Missing required fields: {sorted(missing)}")

        # Mutual exclusivity
        for pair in cls._mutually_exclusive:
            present = [f for f in pair if f in values]
            if len(present) > 1:
                raise ValueError(f"Fields {pair} are mutually exclusive: got {present}")

        # Constrained choices
        for field, allowed in cls._attr_choices.items():
            if field in values and values[field] not in allowed:
                raise ValueError(f"Field '{field}' must be one of {allowed}, got '{values[field]}'")

        # Conditional dependencies
        for field, value_dependencies in cls._attr_dependencies.items():
            if field in values:
                value = values[field]
                # Handle wildcard '*' as a catch-all dependency
                if '*' in value_dependencies:
                    for dep in value_dependencies['*']:
                        if dep not in values:
                            raise ValueError(f"If '{field}' is set, '{dep}' must also be set.")
                if value in value_dependencies:
                    for dep in value_dependencies[value]:
                        if dep not in values:
                            raise ValueError(f"If '{field}' is set to '{value}', '{dep}' must also be set.")

        return values

    @model_validator(mode="after")
    def _validate_nested(self):
        for name in self.__pydantic_fields__:
            val = getattr(self, name)

            if isinstance(val, BaseObj):
                val._validate_nested()

            elif isinstance(val, list):
                for item in val:
                    if isinstance(item, BaseObj):
                        item._validate_nested()

            elif isinstance(val, dict):
                for item in val.values():
                    if isinstance(item, BaseObj):
                        item._validate_nested()

        return self

    @model_serializer(mode='plain')
    def serialize(self) -> dict[str, Any]:
        data = {}
        for field_name, value in self.__dict__.items():
            if isinstance(value, BaseObj):
                # logger.debug(f"Serializing nested BaseObj field '{field_name}' in {self.__class__.__name__}")
                data[field_name] = value  # preserve as BaseObj instance
            else:
                data[field_name] = value  # normal field
        # logger.debug('Serialized data:')
        # for k, v in data.items():
        #     logger.debug(f"  {k}: {v} (type {type(v)})")
        return data

    def set(self, shallow=False, **fields):
        """
        Sets the values of the calling instance's attributes
        to the values in fields. If shallow is False, and if there are
        nested objects, their attributes will also be set.
        If shallow is True, only the attributes of the calling instance
        will be set.
        """
        # logger.debug(f"Setting fields {fields} on {self.__class__.__name__} instance; shallow={shallow}")
        if not shallow:
            self.set_nested(**fields)
        else:
            for key, value in fields.items():
                setattr(self, key, value)  # will fail if key is not in annotations

    def set_nested(self, **fields):
        """
        Recursively sets attributes of nested BaseObj instances.
        This method is called by set() when shallow is False.
        """
        for key in self.__pydantic_fields__:
            value = getattr(self, key)
            # logger.debug(f"Processing field '{key}' with value-type '{type(value)}' in {self.__class__.__name__}")
            if key in fields:
                # logger.debug(f"Setting field '{key}' to '{fields[key]}' in {self.__class__.__name__}")
                setattr(self, key, fields[key])  # will fail if key is not in annotations
            elif isinstance(value, BaseObj):
                # if the value is a BaseObj, set its attributes
                value.set_nested(**fields)
            elif isinstance(value, list):
                # if the value is a list of BaseObj, set its items
                for item in value:
                    if isinstance(item, BaseObj):
                        item.set_nested(**fields)
            elif isinstance(value, dict):
                # if the value is a dict of BaseObj, set its items
                for item in value.values():
                    if isinstance(item, BaseObj):
                        item.set_nested(**fields)

    # Comparison helpers
    def _dict_for_comparison(self) -> dict:
        return {
            k: v for k, v in self.model_dump().items()
            if k not in self._ignore_fields
        }

    def __eq__(self, other: BaseObj | dict) -> bool:
        if self is other: return True
        if isinstance(other, dict):
            return self._dict_for_comparison() == other
        if not isinstance(other, self.__class__) and not isinstance(other, dict):
            return False
        return self._dict_for_comparison() == other._dict_for_comparison()

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented
        keys = self._dict_for_comparison().keys() & other._dict_for_comparison().keys()
        return all(getattr(self, k) <= getattr(other, k) for k in keys) and any(
            getattr(self, k) < getattr(other, k) for k in keys
        )

    def __le__(self, other: object) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented
        keys = self._dict_for_comparison().keys() & other._dict_for_comparison().keys()
        return all(getattr(self, k) <= getattr(other, k) for k in keys)
    
    def __setattr__(self, name, value):
        if name.startswith('_') and name not in self.__private_attributes__:
            # This lets internal tools assign attributes like __class__ etc.
            object.__setattr__(self, name, value)
        else:
            super().__setattr__(name, value)

    # Utility matchers
    def strhash(
        self,
        fields: list[str] | None = None,
        sep: str = "-",
        digest: str | None = None,  # e.g., "md5", "sha256"
        exclude_none: bool = True,
        exclude_ignored: bool = True,
    ) -> str | bytes:
        """
        Generate a hash-like string from selected fields in the model.

        Parameters
        ----------
        fields : list of field names to include; if None, include all.
        sep : string separator for components (default: "-")
        digest : optional name of hashlib algorithm to use ("md5", "sha1", "sha256", etc.)
        exclude_none : if True, skip fields with value None
        exclude_ignored : if True, exclude fields in _ignore_fields

        Returns
        -------
        str or bytes:
            Joined string or digest string (hex or binary, depending on digest)
        """
        data = self.model_dump()
        keys = fields or sorted(data.keys())

        if exclude_ignored:
            keys = [k for k in keys if k not in self._ignore_fields]
        if exclude_none:
            keys = [k for k in keys if data[k] is not None]

        parts = [f"{k}={data[k]}" for k in keys]
        raw_string = sep.join(parts)

        if digest:
            h = hashlib.new(digest)
            h.update(raw_string.encode("utf-8"))
            return h.hexdigest()  # or h.digest() for raw bytes

        return raw_string

    def wildmatch(self, **fields) -> bool:
        """
        Check if the object matches any of the provided fields with wildcards.
        This allows for substring matching in the field names.
        """
        for substr, target_val in fields.items():
            logger.debug(f"Checking requested {substr}={target_val} of type {type(target_val)} in {self.__class__.__name__}")
            match_found = False
            for k, v in self.model_dump().items():
                if substr in k and v == target_val:
                    match_found = True
                    break
            if not match_found:
                return False
        return True
    
    def weak_lt(self,other: object, attr=[]) -> bool:
        """
        Weak less than comparison that only checks fields listed in the attr parameter.
        """
        if not isinstance(other, self.__class__):
            return NotImplemented
        self_dict = self._dict_for_comparison()
        other_dict = other._dict_for_comparison()
        if not self_dict or not other_dict:
            return False
        for key in attr:
            if key in self_dict and key in other_dict:
                if operator.lt(self_dict[key], other_dict[key]):
                    return True
        return False
    
    def dump(self):
        """
        Simple dump of this item's __dict__ in YAML format.

        Returns
        -------
        str :
            yaml-format dump of calling instance
        """
        return yaml.dump(self.model_dump(), default_flow_style=False, sort_keys=False)
    
    def inlist(self,a_list):
        """
        Checks if the calling instance is a member of a_list
        This is a simple membership test that checks if the calling instance
        is equal to any element in a_list.

        Parameters
        ----------
        a_list : list
            the list of objects checked to see if object is a member

        Returns
        -------
        bool :
            True if calling instance is a member of a_list
        """
        for s in a_list:
            if s==self:
                return True
        return False

    def map_attr(self,mapped_attr,key_attr,map):
        """
        Simple cross-attribute mapper. 

        Parameters
        ----------
        mapped_attr : str
            name of attribute to which the mapping result 
            will be applied
        key_attr : str
            name of attribute whose value is mapped
        map : dict
            the map

        Returns
        -------
        None :
            This method modifies the calling instance in place.
            It sets the mapped_attr to the value from the map
            corresponding to the key_attr's value.
            If the key_attr's value is not in the map, it does nothing.
        """
        if not hasattr(self, mapped_attr) or not hasattr(self, key_attr):
            return

        key_value = getattr(self, key_attr)
        if key_value in map:
            setattr(self, mapped_attr, map[key_value])

    def swap_attr(self,attr1,attr2):
        """
        Simple attribute value swapper
        
        Parameters
        ----------
        attr1 : str
            name of first attribute in pair to swap
        attr2 : str
            name of second attribute in pair to swap
        """
        if not hasattr(self, attr1) or not hasattr(self, attr2):
            return

        value1 = getattr(self, attr1)
        value2 = getattr(self, attr2)

        setattr(self, attr1, value2)
        setattr(self, attr2, value1)

    def copy_attr(self, recv_attr, src_attr):
        """
        Simple attribute copier 
        
        Parameters
        ----------
        recv_attr : str
            name of receiver attribute; must exist
        src_attr : str
            name of source attribute; must exist
        """
        if not hasattr(self, recv_attr) or not hasattr(self, src_attr):
            return

        value = getattr(self, src_attr)
        setattr(self, recv_attr, value)

    def assign_obj_to_attr(self, attr, objList, **matchattr):
        """
        Given the dictionary matchattr that maps attribute names in self to corresponding attribute names in the objects in objList, pluck out the element in objList that matches self and assign it to the attr attribute of self.

        Parameters
        ----------
        attr : str
            attribute name that will receive the matched object from objList
        objList : list
            list of objects that is searched
        matchattr : dict
            attribute_in_searched_objects:attribute_in_self
        """
        adict = {k:getattr(self, v) for k, v in matchattr.items()}
        pluckedObj = objList.get(objList.dict_to_condition(adict)) # get returns None, a list of matches, or a single match
        # the only case where we assign is when we get a single match
        if pluckedObj is not None and type(pluckedObj) != objList:
            setattr(self, attr, pluckedObj)

    def update_attr_from_obj_attr(self, attr, obj_attr, attr_of_obj_attr):
        """
        Update an attribute of self from an attribute of one of self's object attributes.

        Parameters
        ----------
        attr : str
            attribute name of self that will be set to value of self.obj_attr.attr_of_obj_attr
        obj_attr : str
            name of object attribute from which the attribute value is taken
        attr_of_obj_attr : str
            name of attribute in obj_attr that is set to
            attr of caller
        """
        setattr(self, attr, getattr(getattr(self, obj_attr), attr_of_obj_attr))

    def update_attr_from_objlist_elem_attr(self, attr, objlist_attr, index_of_obj_in_objlist_attr, attr_of_obj_attr):
        """
        Set value of caller's attribute from an attribute
        of an object in a list of objects
        
        Parameters
        ----------
        attr : str
            attribute name of self that will be set to value of self.objlist_attr[index_of_obj_in_objlist_attr].attr_of_obj_attr
        objlist_attr : str
            name of the caller's attribute that is a list of objects
        index_of_obj_in_objlist_attr : int
            index of the object in the list
        attr_of_obj_attr : str
            name of the attribute in the object that is set to
            attr of caller
        """
        setattr(self, attr, getattr(getattr(self, objlist_attr)[index_of_obj_in_objlist_attr], attr_of_obj_attr))

class GenericListMeta(ABCMeta):
    """ 
    Metaclass for abstract :class:`BaseObjList` to handle generic type arguments.
    
    A subclass of :class:`BaseObjList` must be declared as class MyList(BaseObjList[MyType])
    where MyType is a subclass of :class:`BaseObj`.

    Ultimately, this enables subclasses of :class:`BaseObj` to have Fields that are subclasses of :class:`BaseObjList`.
    We need this because certain objects when put into lists acquire new attributes that are not present in the original object.
    """
    def __new__(mcls, name, bases, namespace, **kwargs):
        cls = super().__new__(mcls, name, bases, namespace, **kwargs)
        # Use __orig_bases__ to introspect generic base with type args
        orig_bases = namespace.get('__orig_bases__', ())
        for base in orig_bases:
            origin = get_origin(base)
            args = get_args(base)
            if origin and args: # and not isinstance(args[0], TypeVar):
                cls._item_type = args[0]
                return cls
        else:
            # Fallback: walk MRO and look for _item_type
            for base in cls.__mro__[1:]:
                inherited_type = getattr(base, "_item_type", None)
                if inherited_type is not None and not isinstance(inherited_type, TypeVar):
                    cls._item_type = inherited_type
                    return cls
        if getattr(cls, "_item_type", None) is None or isinstance(cls._item_type, TypeVar):
            raise TypeError(
                f"{name} must subclass BaseObjList with a concrete Pydantic model type parameter, "
                f"e.g., `class {name}(BaseObjList[MyModel])` or subclass something that does."
            )

T = TypeVar("T", bound=BaseObj)

class BaseObjList(UserList[T], Generic[T], metaclass=GenericListMeta):
    def __init__(self, initlist: Iterable[T] = ()):
        self.data = []
        for item in initlist:
            self.data.append(self._validate_item(item))

    def _validate_item(self, item: T):
        if isinstance(item, dict):
            item = self._item_type(**item)
        elif not isinstance(item, self._item_type):
            raise TypeError(f"{item} is not an instance of {self._item_type}")
        return item

    @classmethod
    def validate(cls, v):
        if isinstance(v, cls):
            return v
        elif isinstance(v, Iterable):
            return cls(v)
        else:
            raise TypeError(f"Expected {cls.__name__} or iterable of {cls.__args__[0].__name__}")

    def append(self, item: BaseObj) -> None:
        self._validate_item(item)
        self.data.append(item)

    def extend(self, other: Iterable[BaseObj]) -> None:
        for item in other:
            self._validate_item(item)
            self.append(item)

    def insert(self, index: int, item: BaseObj) -> None:
        self._validate_item(item)
        self.data.insert(index, item)

    def remove_instance(self, obj_to_remove: BaseObj):
        """
        Remove the exact instance `obj_to_remove` from the list,
        by checking identity (not equality), avoiding __eq__ overhead.
        Raises ValueError if not found.
        """
        for i, item in enumerate(self.data):
            if item is obj_to_remove:
                del self.data[i]
                return
        raise ValueError("Object not found in list")

    def __add__(self, other: Iterable[BaseObj]) -> BaseObjList:
        new = self.__class__(self)
        new.extend(other)
        return new

    def __getitem__(self, index: int) -> BaseObj | BaseObjList:
        result = super().__getitem__(index)
        if isinstance(index, slice):
            return type(self)(result)
        return result
    
    def __iter__(self) -> Iterator[BaseObj]:
        return iter(self.data)

    def __setitem__(self, index: int, item: BaseObj) -> None:
        self._validate_item(item)
        self.data[index] = item

    def __iadd__(self, other: Iterable[BaseObj]) -> BaseObjList:
        self.extend(other)
        return self
        
    @abstractmethod
    def describe(self) -> str:
        """
        Abstract method to describe the contents of the BaseObjList.
        Subclasses should implement this method to provide a meaningful description.
        """
        raise NotImplementedError("Subclasses must implement the describe method.")

    def filter(self, func: Callable[[BaseObj], bool]) -> BaseObjList:
        """
        Filters the list of BaseObj instances based on a provided function.

        Parameters
        ----------
        func : Callable[[BaseObj], bool]
            A function that takes a BaseObj instance and returns True if it should be included in the filtered list.

        Returns
        -------
        BaseObjList
            A new BaseObjList containing only the objects for which func returned True.
        """
        return self.__class__(list(filter(func, self.data)))

    def ifilter(self, func: Callable[[BaseObj], bool]) -> list[int]:
        """
        Filters the list of BaseObj instances based on a provided function and returns the list of indices for items passing the filter

        Parameters
        ----------
        func : Callable[[BaseObj], bool]
            A function that takes a BaseObj instance and returns True if it should be included in the filtered list.

        Returns
        -------
        BaseObjList
            A new BaseObjList containing only the objects for which func returned True.
        """
        return [i for i, item in enumerate(self.data) if func(item)]

    def get(self, func: Callable[[BaseObj], bool]) -> BaseObj | BaseObjList | None:
        """
        Special implementation of filter that returns a single object if there is only one match;
        if there are multiple matches, all are returned in a list;
        if there are no matches, None is returned.
        """
        matches = self.filter(func=func)
        if not matches:
            return None
        if len(matches) == 1:
            return matches[0]
        return matches

    def iget(self, func: Callable[[BaseObj], bool]) -> int | list[int] | None:
        """
        Special implementation of ifilter that returns the index of a single index if there is only one match;
        if there are multiple matching indices, all are returned in a list;
        if there are no matches, None is returned.
        """
        matches = self.ifilter(func=func)
        if not matches:
            return None
        if len(matches) == 1:
            return matches[0]
        return matches

    @staticmethod
    def dict_to_condition(condition_dict: dict[str, Callable[[BaseObj], bool]], conjunction: str = 'and') -> Callable[[BaseObj], bool]:
        """
        Converts a dictionary of conditions into a lambda function.
        The dictionary keys are attribute names, and the values are the conditions.
        
        Returns a callable function that can be used in the filter method.
        """
        # Create a lambda function that checks each key-value condition
        if conjunction == 'or':
            return lambda obj: any(
                (getattr(obj, key) in value if isinstance(value, list) else (getattr(obj, key) == value if callable(value) is False else value(getattr(obj, key))))
                for key, value in condition_dict.items()
            )
        return lambda obj: all(
            # Check if the value is a list and perform an "in" check if it is
            (getattr(obj, key) in value if isinstance(value, list) else (getattr(obj, key) == value if callable(value) is False else value(getattr(obj, key))))
            for key, value in condition_dict.items()
        )

    def remove_and_return(self, item: BaseObj) -> BaseObj | None:
        """
        Remove an item from the list and return it.

        Parameters
        ----------
        item : BaseObj
            The item to remove from the list.

        Returns
        -------
        BaseObj | None :
            The removed item if it was found, otherwise None.
        """
        if item in self.data:
            self.data.remove(item)
            return item
        
    def give_item(self, item: BaseObj, target_list: BaseObjList) -> bool:
        """
        Transfer an item from this list to another BaseObjList.

        Parameters
        ----------
        item : BaseObj
            The item to transfer.
        target_list : BaseObjList
            The list to transfer the item to.

        Returns
        -------
        bool
            True if the transfer was successful, False otherwise.
        """
        if item in self.data:
            self.remove(item)
            target_list.append(item)
            return True
        return False
    
    def take_item(self, item: BaseObj, source_list: BaseObjList) -> bool:
        """
        Take an item from another BaseObjList and add it to this list.

        Parameters
        ----------
        item : BaseObj
            The item to take.
        source_list : BaseObjList
            The list to take the item from.

        Returns
        -------
        bool
            True if the transfer was successful, False otherwise.
        """
        if item in source_list:
            source_list.remove(item)
            self.append(item)
            return True
        return False

    def set(self, shallow=False, **fields):
        """
        Element attribute-setter
        
        Parameters
        ----------
        fields : dict
            attribute-name:value pairs used to set the attributes
            of the calling instance
        """
        for item in self.data:
            item.set(shallow=shallow, **fields)

    def prune_exclusions(self, exclude_func: Callable[[BaseObj], bool]) -> BaseObjList:
        """
        Prune elements from the calling instance based on a custom exclusion function.

        Parameters
        ----------
        exclude_func : Callable[[BaseObj], bool]
            A function that takes a BaseObj and returns True if it should be excluded.

        Returns
        -------
        BaseObjList :
            items removed from calling instance
        """
        remove_us = self.filter(exclude_func)
        for item in remove_us:
            self.remove(item)
        logger.debug(f"Pruning {len(remove_us)} items from {self.__class__.__name__} with custom criteria")
        return remove_us

    def sort(self, by=None, reverse=False):
        """
        ObjList sort function, a simple override of UserList.sort()

        Parameters
        ----------
        by : list, optional
            names of attributes used to execute __lt__ comparisons between
            elements of calling instance
        reverse : bool, optional
            if true, reverse the sense of the sort
        """
        if not by:
            self.data.sort(reverse=reverse)
        else:
            key = operator.attrgetter(*by)
            self.data.sort(key=key, reverse=reverse)

    def uniqattrs(self, attrs: list[str] = [], with_counts: bool = False) -> dict[str, list]:
        """
        Generates a dictionary of list of unique values for each 
        attribute 
        
        Parameters
        ----------
        attrs : list
            names of attributes for which uniqueness is requested
        with_counts : bool
            if true, include the counts of occurences of each unique
            value such that each dictionary value is now a list of
            tuples, each of which is a concatentation of the unique 
            value the count of its occurence among all elements
            in the calling instance
        
        Returns
        -------
        dict :
            lists of unique values for each attribute key, or, 
            if with_counts is true, list of (value,count)
        """
        uattrs = {k: [] for k in attrs}
        for item in self:
            for a in attrs:
                v = getattr(item, a)
                if with_counts:
                    try:
                        idx = [x[0] for x in uattrs[a]].index(v)
                    except:
                        uattrs[a].append([v, 0])
                        idx = -1
                    uattrs[a][idx][1] += 1
                else:
                    if not v in uattrs[a]:
                        uattrs[a].append(v)
        return uattrs

    def binnify(self, fields: list[str] = []) -> dict[str, BaseObjList]:
        """
        Simple binning of all elements by unique hashes of values of fields

        Parameters
        ----------
        fields : list
            attribute names used to build the hash to test for uniqueness
        
        Returns
        -------
        dict :
            lists of items for each unique key

        """
        bins: dict[str, BaseObjList] = {}
        for item in self:
            key = item.strhash(fields=fields)
            if not key in bins:
                bins[key] = self.__class__([])
            bins[key].append(item)
        # for k,v in bins.items():
        #     if len(v)>1:
        #         logger.debug(f'binnify: bin {k} has {len(v)} items: {v}')
        return bins

    def puniq(self, fields: list[str] = []) -> bool:
        """
        Simple test that all elements of self are unique among fields 
        
        Parameters
        ----------
        fields : list, optional
            attribute names used to build the hash to test for uniqueness
        
        Returns
        -------
        bool : True if all elements are unique
        
        """
        bins = self.binnify(fields=fields)
        for b, v in bins.items():
            if len(v) > 1:
                logger.debug(f'puniquify: bin {b} has {len(v)} items: {v}')
        return len(bins) == len(self)

    def puniquify(self, attrs: list[str] = [], stash_attr_name: str = 'ORIGINAL_ATTRIBUTES') -> None:
        """
        Systematic attribute altering to make all elements unique
        
        There may be a set of attributes for which no two elements may
        have the exact same set of respective values.  This method 
        scans the calling instance for such collisions and, if any
        is found, it adds one to the value of the attribute named
        in the first element of the 'fields' list (assumes this
        attribute is numeric!).  This could lead to other collisions
        so multiple passes through the calling instance are made
        until there are no more collisions.  Each such value 
        change results in storing the original values in a new
        attribute.

        Parameters
        ----------
        fields : list, optional
            attribute names used to build the hash to test for uniqueness;
            if unset, all attributes are used
        stash_attr_name : str, optional
            name given to a new dict attribute used to store all original 
            attribute name:value pairs

        """
        local_attr_copy = attrs.copy() if attrs else []
        logger.debug(f'puniquify: Checking uniqueness for {local_attr_copy}; first binning:')
        bins = self.binnify(fields=local_attr_copy)
        logger.debug(f' # bins {len(bins)} vs # items {len(self.data)}')
        if len(bins) == len(self.data):
            # all items are unique, nothing to do
            return
        stillworking = True
        num_iter = 0
        while stillworking and num_iter < 10:
            stillworking = False
            num_iter += 1
            for v in bins.values():
                if len(v) > 1:  # this key has more than one item
                    stillworking = True
                    for d in v[1:]:  # skip the first one, it is the original
                        # stash the original attribute values in a dict under stash_attr_name
                        # but only do this ONCE per object so we are only saving its original values
                        if not hasattr(d, stash_attr_name):
                            logger.debug(f'puniquify: Stashing original attributes {local_attr_copy} for {str(d)} under {stash_attr_name}')
                            setattr(d, stash_attr_name, {k: getattr(d, k) for k in local_attr_copy})
                        while d.strhash(local_attr_copy) in bins:
                            logger.debug(f'puniquify: Collision detected for {d.strhash(local_attr_copy)}; incrementing {repr(local_attr_copy[0])}')
                            # increment the first value until the hash is unique
                            # this assumes the first field in local_attr_copy is numeric
                            value_to_increment = getattr(d, local_attr_copy[0])
                            try:
                                value_to_increment = operator.add(value_to_increment, 1)
                                setattr(d, local_attr_copy[0], value_to_increment)
                            except TypeError:
                                raise TypeError(f"Field '{local_attr_copy[0]}' must be addable to int for uniquification.")
            if stillworking:
                # re-bin the items
                logger.debug(f'Rebinning at iteration {num_iter}')
                bins = self.binnify(fields=local_attr_copy)
        if num_iter == 10:
            raise RuntimeError(f"puniquify: Unable to uniquify after {num_iter} iterations; giving up.")
        assert(self.puniq(fields=local_attr_copy))

    def map_attr(self, mapped_attr: str, key_attr: str, map: dict) -> None:
        """
        Simple cross-attribute mapper. 

        Parameters
        ----------
        mapped_attr : str
            name of attribute to which the mapping result 
            will be applied
        key_attr : str
            name of attribute whose value is mapped
        map : dict
            the map

        Returns
        -------
        None :
            This method modifies the calling instance in place.
            It sets the mapped_attr to the value from the map
            corresponding to the key_attr's value.
            If the key_attr's value is not in the map, it does nothing.
        """
        for item in self:
            item.map_attr(mapped_attr, key_attr, map)

    def remove_duplicates(self, fields: list[str] = []) -> None:
        """
        Removes duplicates from the calling instance
        The duplicates are determined by the hash of the values
        of the attributes named in the fields list. 

        Parameters
        ----------
        fields : list, optional
            attribute names used to build the hash to test for uniqueness;
            if unset, all attributes are used
        """
        bins=self.binnify(fields=fields)
        # logger.debug(f'remove_duplicates: bin counts {[len(b) for b in bins.values()]}')
        self.clear()
        assert len(self)==0
        for b in bins.values():
            # logger.debug(f'preserving {str(b[0])}')
            self.append(b[0])
            # for c in b[1:]:
            #     logger.debug(f'discarding {str(c)}')

    def assign_objs_to_attr(self, attr: str, objList: BaseObjList, **matchattr) -> BaseObjList:
        """
        Assigns the single object from objList whose
        attributes match the matchattr dict to the 
        calling instances attr attribute, but only if
        the result of the match-search is valid,
        otherwise assign None to attr

        Parameters
        ----------
        attr : str
            attribute name
        objList : list
            list of objects that is searched
        matchattr : dict
            attribute:values used in searching the 
            list of objects
        """
        for item in self:
            item.assign_obj_to_attr(attr, objList, **matchattr)
        delete_us = [item for item in self if getattr(item, attr) is None]
        for item in delete_us:
            self.remove(item)
        return self.__class__(delete_us)

    def update_attr_from_obj_attr(self, attr: str, obj_attr: str, attr_of_obj_attr: str) -> None:
        """
        Set value of attributes of all elements of caller
        from another attribute of a separate object

        Parameters
        ----------
        attr : str
            attribute name
        obj_attr : str
            name of object attribute from which the attribute value is taken
        attr_of_obj_attr : str
            name of attribute in obj_attr that is set to
            attr of caller
        """
        for item in self:
            item.update_attr_from_obj_attr(attr, obj_attr, attr_of_obj_attr)

    def update_attr_from_objlist_elem_attr(self, attr: str, objlist_attr: str, index_of_obj_in_objlist_attr: int, attr_of_obj_attr: str) -> None:
        """
        Set value of caller's attribute from an attribute
        of an object in a list of objects
        
        Parameters
        ----------
        attr : str
            attribute name to be updated
        objlist_attr : str
            name of the caller's attribute that is a list of objects
        index_of_obj_in_objlist_attr : int
            index of the object in the list
        attr_of_obj_attr : str
            name of the attribute in the object that is set to
            attr of caller
        """
        for item in self:
            item.update_attr_from_objlist_elem_attr(attr, objlist_attr, index_of_obj_in_objlist_attr, attr_of_obj_attr)