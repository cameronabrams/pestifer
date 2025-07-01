# Author: Cameron F. Abrams <cfa22@drexel.edu>.
""" Namespace objects with attribute controls

The BaseObj class is a derivative of argparse's
very nice Namespace class, and the ObjList class as a derivative
of collection's UserList.

Attribute controls are enforced using class attributes that 
allow specification of 

    * attributes that are "required" or "optional";
    * pairs of attributes that are mutually exclusive;
    * attributes whose values must be selected from
      some predefined set of choices;
    * attributes that depend on other attributes;
    * attributes that are ignored for purposes of comparing
      the overall "value" of any instances

A class with methods useful for manipulating lists of instances 
of BaseObj's or their derivatives is also defined here.

"""
import operator
import yaml
import logging
import inspect
from collections import UserList
logger=logging.getLogger(__name__)
from argparse import Namespace
from functools import singledispatchmethod

class BaseObj(Namespace):
    """
    A class defining a namespace with custom attribute controls.

    Required attributes are those that must have a value assigned
    upon instance creation.  Optional attributes *may* have a 
    value assigned but need not.  Any pair of optional attributes
    can be designated as "mutually exlusive", meaning that if one
    is assigned a value at instance creation, the other is not
    allowed to exist.  Values for any attribute can be restricted 
    to one of given set of values.  Any optional attribute may
    by assigned other optional attributes that must also have 
    values in order to have a value.  Finally, any attribute may
    be designated as silent when relative comparisons between
    object instances are made.
   
    """

    req_attr=[]
    """
    List of required attributes as strings. These attributes must have values assigned at instance creation time, otherwise an AssertionError is raised
    """

    opt_attr=[]
    """
        List of optional attributes as strings. These attributes may have values assigned at instance creation time, but need not; if they do not have values, they are simply not present in the instance's __dict__.
    """

    alt_attr=[]
    """
        Lists of pairs of mutually exclusive attributes.
        Each pair is a list of two strings, each of which is an attribute label. If one of the attributes in the pair has a value assigned at instance creation time, the other is not allowed to have a value. If neither attribute in the pair has a value, then both are simply
        not present in the instance's __dict__.
        e.g., [['is_edible','is_flavorable'], ...]
        This means that if 'is_edible' has a value, 'is_flavorable'
        cannot have a value, and vice versa.
        If neither attribute in the pair has a value, then both are simply not present in the
    """

    attr_choices={}
    """
        Dictionary of attributes with limited set of possible values.
        The keys are attribute labels, and the values are lists of possible
        values for the respective attributes. If an attribute is present
        in the instance's __dict__, its value must be one of the values in the
        list. If the attribute is not present, it is simply not present
        in the instance's __dict__.
        e.g., {'flavor':['sweet','sour','bitter'], ...}
        This means that if 'flavor' has a value, it must be one of the
        values in the list, e.g., 'sweet', 'sour', or 'bitter'.
        If 'flavor' is not present, it is simply not present in the instance's __dict__.
    """

    opt_attr_deps={}
    """ 
        Dictionary of optional attributes that have required attributes
        that must also have values assigned at instance creation time.
        The keys are optional attribute labels, and the values are lists of
        required attribute labels that must also have values assigned
        at instance creation time if the optional attribute has a value.
        e.g., {'has_sauce':['sauce_type'], ...}
        This means that if 'has_sauce' has a value, 'sauce_type'
        must also have a value assigned at instance creation time.
        If 'has_sauce' is not present, it is simply not present in the instance's __dict__.
    """
    
    ignore_attr=[]
    """
        List of attributes that are ignored when comparing instances.
        This is a list of attribute labels that are not considered when
        comparing instances of this class. If an attribute is in this list,
        its value is not considered when comparing instances of this class.
        e.g., ['timestamp','id']
        This means that if 'timestamp' or 'id' are present, their values
        are not considered when comparing instances of this class.
        If an attribute is not present, it is simply not present in the instance's __dict__.
    """
    
    @singledispatchmethod
    def __init__(self,input_obj):
        """
        Default constructor fails; all constructions are argument-type dependent 
        """
        msg=f'Cannot initialize {self.__class__} from object type {type(input_obj)}'
        logger.error(msg)
        raise TypeError(msg)

    @__init__.register(Namespace)
    def _from_namespace(self,ns):
        super().__init__(**(ns.__dict__))

    @__init__.register(dict)
    def _from_dict(self,input_dict):
        """
        BaseObj constructor when single positional argument is a dict
        Parameters
        ----------
        input_dict : dict
            dictionary of input attribute:value pairs
        """
        # mutually exclusive attribute labels should not 
        # appear in the list of required attributes
        assert all([(not x in self.req_attr and not y in self.req_attr) for x,y in self.alt_attr]),"Mutually exclusive attributes should not appear in the required list"
        prep_dict={}
        for k,v in input_dict.items():
            if not self.req_attr and not self.opt_attr:
                # let's treat an undifferentiated BaseObj like a Namespace
                prep_dict[k]=v
            elif k in self.req_attr or k in self.opt_attr:
                prep_dict[k]=v
        # all required attributes have values
        assert all([att in prep_dict.keys() for att in self.req_attr]),f"Not all required attributes have values\nRequired: {self.req_attr}\nProvided: {list(prep_dict.keys())}"
        # all mutual-exclusivity requirements are met
        assert all([
            ((x[0] in prep_dict.keys() and not x[1] in prep_dict.keys()) or 
             (x[1] in prep_dict.keys() and not x[0] in prep_dict.keys()) or not (x[0] in prep_dict.keys() or x[1] in prep_dict.keys()))
            for x in self.alt_attr]),"Mutual exclusivity requirements unmet"
        # all attributes with limited set of possible values have valid values
        try:
            assert all([prep_dict[x] in choices for x,choices in self.attr_choices.items()]),"Invalid choices in one or more attributes"
        except:
            raise AssertionError(f'Trouble initializing instance of {self.__class__}')
        # for each optional attribute present, all required dependent attributes are also present
        for oa,dp in self.opt_attr_deps.items():
            if oa in prep_dict.keys():
                assert all([x in prep_dict.keys() for x in dp]),"Dependent required attributes of optional attributes not present"
        super().__init__(**prep_dict)

    def __eq__(self,other):
        """
        Defines the equality operator for BaseObj instances.
        The two instances are considered equal if their attribute values are equal.
        This includes all required attributes and any common optional attributes.

        Parameters
        ----------
        other : BaseObj
            other BaseObj instance
        """
        if not other:
            return False
        if id(self)==id(other):
            return True
        if not self.__class__==other.__class__:
            return False
        attr_list=self.req_attr+[x for x in self.opt_attr if (x in other.__dict__ and x in self.__dict__)]
        for x in self.ignore_attr:
            if x in attr_list:
                attr_list.remove(x)
        test_list=[]
        for k in attr_list:
            test_list.append(self.__dict__[k]==other.__dict__[k])
        return all(test_list)
    
    def __lt__(self,other):
        """
        Defines the less-than operator for BaseObj instances.

        Parameters
        ----------
        other : BaseObj
            other BaseObj instance
        """
        if not other:
            return False
        if not self.__class__==other.__class__:
            return False
        attr_list=self.req_attr+[x for x in self.opt_attr if (x in other.__dict__ and x in self.__dict__)]
        for x in self.ignore_attr:
            if x in attr_list:
                attr_list.remove(x)
        lt_list=[self.__dict__[k]<other.__dict__[k] for k in attr_list if hasattr(self.__dict__[k],'__lt__')]
        le_list=[self.__dict__[k]<=other.__dict__[k] for k in attr_list if hasattr(self.__dict__[k],'__lt__')]
        return all(le_list) and any(lt_list)
    
    def __le__(self,other):
        """
        Defines the less-than-or-equal-to operator for BaseObj instances.
        The two instances are considered less than or equal to each other if their attribute values are equal or if one is less than the other.

        Parameters
        ----------
        other : BaseObj
            other BaseObj instance
        """
        if not other:
            return False
        if not self.__class__==other.__class__:
            return False
        attr_list=self.req_attr+[x for x in self.opt_attr if (x in other.__dict__ and x in self.__dict__)]
        for x in self.ignore_attr:
            if x in attr_list:
                attr_list.remove(x)
        le_list=[self.__dict__[k]<=other.__dict__[k] for k in attr_list if hasattr(self.__dict__[k],'__lt__')]
        return all(le_list)
    
    def weak_lt(self,other,attr=[]):
        """
        Defines a weak less-than operator for BaseObj instances.  Only those attributes listed in the attr
        parameter are considered in the comparison

        Parameters
        ----------
        other : BaseObj
            other BaseObj instance
        attr : list, optional
            attributes used for the less-than comparison
        """
        lt_list=[self.__dict__[k]<other.__dict__[k] for k in attr if hasattr(self.__dict__[k],'__lt__')]
        le_list=[self.__dict__[k]<=other.__dict__[k] for k in attr if hasattr(self.__dict__[k],'__lt__')]
        return all(le_list) and any(lt_list)

    def strhash(self,fields=[]):
        """
        Generates a string hash of the calling instance's attributes.
        The hash is a string of the values of the attributes sorted by their names.
        
        Parameters
        ----------
        fields : list, optional
            attributes to be included in the identifier
        """
        if fields:
            items=sorted([(k,x) for k,x in self.__dict__.items() if k in fields], key=lambda it: it[0])
        else:
            items=sorted(self.__dict__.items(), key=lambda it: it[0])
        result='-'.join([str(x[1]) for x in items if x!=[]])
        return result
    
    def matches(self,**fields):
        """
        Defines a non-strict match operator for BaseObj instances.
        The two instances are considered a match if all key:val pairs in fields
        match the corresponding key:val pairs in the calling instance's attributes.

        Parameters
        ----------
        field : dict, optional
            attribute:value pairs that are checked against those of the object
        
        Returns
        -------
        bool :
            True if calling instance attribute values match the fields input
        """
        for k,v in fields.items():
            if not k in self.__dict__:
                # logger.debug(f'missing key {k}')
                return False
            if v!=self.__dict__[k]:
                # logger.debug(f'wrong value ({self.__dict__[k]}) for key {k}; wanted {v}')
                return False
        return True
    
    def allneg(self,**fields):
        """ 
        Determines if the calling instance has NO attributes
        that match the key:val pairs in fields.
        This is a strict no-match operator, meaning that if any
        key:val pair in fields matches the corresponding key:val pair
        in the calling instance's attributes, the method returns False.
        
        Parameters
        ----------
        field : dict, optional
            attribute:value pairs that are checked against those of the object
        
        Returns
        -------
        bool :
            True if EXACTLY NONE of the calling instance attribute values match 
            the fields input
        """
        goodfields={}
        for k,v in fields.items():
            if k in self.__dict__:
                goodfields[k]=v
        if len(goodfields)==0: return True
        return all([x!=y for x,y in [(goodfields[k],self.__dict__[k]) for k in goodfields.keys()]])

    def wildmatch(self,**fields):
        """ 
        Defines a wildcard match operator for BaseObj instances.
        The two instances are considered a match if all key:val pairs in fields
        match the corresponding key:val pairs in the calling instance's attributes,
        such that the keys in fields are substrings of the attribute names

        Parameters
        ----------
        field : dict, optional
            attribute:value pairs that are checked against those of the object, except
            each attribute name here can match as a substring to any object attribute
            name
        
        Returns
        -------
        bool :
            True if calling instance attribute values match the fields input
        """
        for k,v in fields.items():
            for tk in self.__dict__.keys():
                if k in tk:  # passed in key is a substring of objects attribute name
                    logger.debug(f'wildmatch matches key {k} to attr name {tk}')
                    break
            else:
                return False
            if v!=self.__dict__[tk]:
                # logger.debug(f'wrong value ({self.__dict__[k]}) for key {k}; wanted {v}')
                return False
        return True

    def dump(self):
        """
        Simple dump of this item's __dict__ in YAML format.

        Returns
        -------
        str :
            yaml-format dump of calling instance
        """
        retdict={}
        retdict['instanceOf']=type(self).__name__
        retdict.update(self.__dict__)
        return yaml.dump(retdict)
    
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
        """
        if map:
            key=self.__dict__[key_attr]
            if key in map:
                val=map[key]
                self.__dict__[mapped_attr]=val

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
        if hasattr(self,attr1) and hasattr(self,attr2):
            val1=getattr(self,attr1)
            val2=getattr(self,attr2)
            setattr(self,attr2,val1)
            setattr(self,attr1,val2)
            if not hasattr(self,'swapped'):
                self.swapped={}
            self.swapped[attr1]=attr2
        else:
            if not hasattr(self,attr1):
                raise AttributeError(f'attribute 1 {attr1} not found.')
            if not hasattr(self,attr2):
                raise AttributeError(f'attribute 2 {attr2} not found.')

    def copy_attr(self,recv_attr,src_attr):
        """
        Simple attribute copier 
        
        Parameters
        ----------
        recv_attr : str
            name of receiver attribute; must exist
        src_attr : str
            name of source attribute; must exist
        """
        if hasattr(self,recv_attr) and hasattr(self,src_attr):
            if not hasattr(self,'copied_over'):
                self.copied_over={}
            self.copied_over[recv_attr]=getattr(self,recv_attr)
            setattr(self,recv_attr,getattr(self,src_attr))
        else:
            if not hasattr(self,recv_attr):
                raise AttributeError(f'receiving attribute {recv_attr} not found.')
            if not hasattr(self,src_attr):
                raise AttributeError(f'source attribute {src_attr} not found.')

    def set(self,shallow=False,**fields):
        """
        Sets the values of the calling instance's attributes
        to the values in fields.  If shallow is False, and if there are
        nested objects, their attributes will also be set.
        If shallow is True, only the attributes of the calling instance
        will be set.

        Parameters
        ----------
        fields : dict, optional
            dictionary of attribute-name:value pairs used to set 
            object attribute values respectively
        """
        for k,v in fields.items():
            # first if k is an attribute, set its value
            if k in self.__dict__:
                self.__dict__[k]=v
            # if k is an attr of any object of self
            if not shallow:
                for sk,sv in self.__dict__.items():
                    if hasattr(sv,'req_attr') and hasattr(sv,k):
                        sv.set(**fields)
                    elif hasattr(sv,'objliststamp'):
                        for subitem in sv:
                            if hasattr(subitem,'req_attr') and hasattr(subitem,k):
                                subitem.set(**fields)

    def assign_obj_to_attr(self,attr,objList,**matchattr):
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
        # assert getattr(self,attr)==None
        adict={k:getattr(self,v) for k,v in matchattr.items()}
        myObj=objList.get(**adict)
        if myObj!=None and myObj!=[]:
            setattr(self,attr,myObj)
        else:
            logger.debug(f'Unable assign attribute \'{attr}\' of \'{type(self).__name__}\'')
            logger.debug(f'   \'{type(self).__name__}\' instance attributes:')
            for k,v in self.__dict__.items():
                logger.debug(f'        {k}: {v} ({type(v).__name__})')
            logger.debug(f'   matchattr:')
            for k,v in matchattr.items():
                logger.debug(f'        {k}: {v} ({type(v).__name__})')
            logger.debug(f'   adict:')
            for k,v in adict.items():
                logger.debug(f'        {k}: {v} ({type(v).__name__})')
            logger.debug(f'from list of type {type(objList).__name__}')
            tmp=objList[1]
            logger.debug(f'    Attributes of type {type(tmp).__name__}:')
            for k,v in tmp.__dict__.items():
                logger.debug(f'        {k} ({type(v).__name__})')
            # raise ValueError(f'stop')
    
    def update_attr_from_obj_attr(self,attr,obj,obj_attr):
        """
        Set value of caller's attribute from another
        attribute of a separate object

        Parameters
        ----------
        attr : str
            attribute name
        obj : str
            name of object from which the attribute value is taken
        obj_attr : str
            name of attribute in obj that is set to 
            attr of caller
        """
        # attr_obj=getattr(self,obj)
        # logger.debug(f'accepting {obj_attr} from {obj} {str(attr_obj)} into {type(self)}({id(self)})')
        setattr(self,attr,getattr(getattr(self,obj),obj_attr))
    
    def update_attr_from_objlist_elem_attr(self,attr,objlist,index,obj_attr):
        """
        Set value of caller's attribute from an attribute
        of an object in a list of objects
        
        Parameters
        ----------
        attr : str
            attribute name
        objlist : list
            list of objects
        index : int
            index of object in list
        obj_attr : str
            attribute name in object
        """
        setattr(self,attr,getattr(getattr(self,objlist)[index],obj_attr))

class CloneableObj(BaseObj):
    """
    A class defining a custom namespace that can be cloned 
    
    This BaseObj specialization gives instances the ability to
    be cloned and to remember their source object.  Attribute
    values can optionally be changed during the cloning
    operation.  The attribute ``clone_of`` is added as an optional
    attribute; any BaseObj instance with an attribute named
    ``clone_of`` is interpreted as a clone.

    """
    opt_attr=BaseObj.opt_attr+['clone_of']
    def clone(self,**options):
        """
        Generates and returns a clone of the calling instance
        
        Parameters
        ----------
        options : dict
            dictionary of attribute-name:value pairs that are 
            explicitly assigned in the clone
        
        Returns
        -------
        obj :
            The new cloned object
        """
        input_dict={k:v for k,v in self.__dict__.items() if k in self.req_attr}
        input_dict.update({k:v for k,v in self.__dict__.items() if k in self.opt_attr})
        input_dict.update(options)
        input_dict['clone_of']=self
        x=self.__class__(input_dict)
        return x
    
    def is_clone(self):
        """
        Determines if the calling instance is a clone
        A clone is defined as an instance that has a 'clone_of' attribute

        Returns
        -------
        bool :
            True if calling instance is a clone
        """
        return 'clone_of' in self.__dict__
    def get_original(self):
        """
        Returns identifier of object the calling instance is a clone of 
        or None if calling instance is not a clone 
        
        Returns
        -------
        obj :
            identifier of object the calling instance is a clone of
        None :
            calling instance is not a clone
        """
        if self.is_clone():
            return self.clone_of
        else:
            return None
        
class AncestorAwareObj(CloneableObj):
    """
    A class defining custom, cloneable namespaces that can exists in a hierarchy 
    
    In cases where an object instance attribute is *another* object, like a BaseObj or
    any derivative thereof, one often would like to allow the "child" object to have
    direct knowledge of its "ancestry".  This class introduces the optional
    attribute 'ancestor_obj' to the CloneableObj class, and defines methods that
    allow calling instances to set this attribute.
    
    """
    opt_attr=CloneableObj.req_attr+['ancestor_obj']
    ignore_attr=CloneableObj.ignore_attr+['ancestor_obj']
    
    def claim_self(self,stamp):
        """
        Applies the stamp to the calling instance's ``ancestor_obj``
        attribute 
        
        Parameters
        ----------
        stamp : obj
            object assigned to the ancestor_obj attribute
        """
        self.ancestor_obj=stamp
        # self.__dict__['ancestor_obj']=stamp

    def claim_descendants(self,stamp):
        """
        Recursively instructs calling instance to stamp all of 
        its descendant objects 
        
        Parameters
        ----------
        stamp : obj
            object assigned to all ancestor_obj attributes of any
            ancestor-aware attributes
        """
        self.claim_self(stamp)
        for attr_name,obj in self.__dict__.items():
            if not attr_name=='ancestor_obj': # very important! or, could stamp self last?
                classes=inspect.getmro(type(obj))
                if AncestorAwareObj in classes or AncestorAwareObjList in classes:
                    obj.claim_descendants(stamp)

class StateInterval(AncestorAwareObj):
    """
    A class defining a custom, ancestor-aware namespace that encodes the 
    state between two bounding values
    
    In reference to an arbitrary list, this class defines a construction
    that assigns a value for ``state`` to items in the list inclusively between
    two indices.  This is an ancestor-aware namespace with the additional
    attributes ``state`` and ``bounds``

    """
    req_attr=AncestorAwareObj.req_attr+['state','bounds']
    opt_attr=AncestorAwareObj.opt_attr+['build']
    def __init__(self,input_dict):
        if not 'build' in input_dict:
            input_dict['build']=False
        super().__init__(input_dict)
    
    def declare_buildable(self):
        """
        Sets the build attribute to True
        """
        self.build=True

    def increment_rightbound(self):
        """
        Increments the position-1 element of the bounds attribute
        of the calling instance
        """
        self.bounds[1]+=1
    
    def num_items(self):
        """
        Returns a calculated value of the number of items inclusively
        between the two bounds, assuming we are referring to a container
        object with sequential indicies 
        
        Returns
        -------
        int :
           rightbound minus leftbound plus 1
        """
        return self.bounds[1]-self.bounds[0]+1
    
    def pstr(self):
        """
        Returns a pretty string version
        """
        return f'{self.state}({self.bounds[1]-self.bounds[0]+1})'

class ObjList(UserList):
    """
    List of BaseObjs or derivatives thereof
    
    When collected into lists, BaseObj instances can acquire collective
    importance that needs to be handled.  This class allows for filtering,
    sorting, and other functionalities on such lists.

    """
    objliststamp=True
    def __init__(self,data):
        """
        Standard initialization of the UserList
        """
        super().__init__(data)

    def filter(self,**fields):
        """
        Creates a returns a new list containing objects
        whose attributes match the fields dictionary

        This method uses the ``matches()`` instance method
        of the BaseObj class

        Parameters
        ----------
        fields : dict
            attribute-name:value pairs used to search for
            hits in the calling instance
        
        Returns
        -------
        list :
             list of elements matching fields
        """
        retlist=self.__class__([])
        for r in self:
            if r.matches(**fields): retlist.append(r)
        return retlist
    
    def negfilter(self,**fields):
        """
        Applies a negated filter to the calling instance
        Creates a returns a new list containing objects
        whose attributes do NOT match the fields dictionary
        This method uses the ``allneg()`` instance method
        of the BaseObj class
        
        Parameters
        ----------
        fields : dict
            attribute-name:value pairs used to search for
            hits in the calling instance

        Returns
        -------
        list :
            list of elements whose attributes do NOT match fields
        """
        retlist=self.__class__([])
        for r in self:
            if r.allneg(**fields): retlist.append(r)
        return retlist
        

    def ifilter(self,**fields):
        """
        Creates a returns a new list containing indices
        of objects whose attributes match the fields dictionary
        This method uses the ``matches()`` instance method
        of the BaseObj class
        
        Parameters
        ----------
        fields : dict
            attribute-name:value pairs used to search for
            'hits' in the calling instance
        
        Returns
        -------
        list
            list of indices matching fields
        """
        retlist=[]
        for i,r in enumerate(self):
            if r.matches(**fields): retlist.append(i)
        return retlist
    
    def get(self,**fields):
        """
        Special implementation of filter
        
        get returns a single object if there is only one match;
        if there are multiple matches, all are returned in a list;
        if there are no matches, an empty list is returned.
        
        Parameters
        ----------
        fields : dict
            attribute-name:value pairs used to search for
            'hits' in the calling instance

        Returns
        -------
        None :
            if no elements matching fields are found in calling instance
        obj :
            if one element matching fields is found, this is it
        list :
            list of all elements matching fields if more than one matches
        """
        R=self.filter(**fields)
        if len(R)==0:
            return type(self)([])
        elif len(R)==1:
            return R[0]
        else:
            return R
        
    def iget(self,**fields):
        """
        Special implementation of ifilter
        iget returns a single index if there is only one match;
        if there are multiple matches, all indices are returned in a list;
        if there are no matches, None is returned.

        Parameters
        ----------
        fields : dict
            attribute-name:value pairs used to search for
            hits in the calling instance

        Returns
        -------
        None :
            if no elements matching fields are found in calling instance
        int :
            if one element matching fields is found, this is its index
        list :
            list of all indices matching fields if more than one matches
        """
        I=self.ifilter(**fields)
        if len(I)==0:
            return None
        elif len(I)==1:
            return I[0]
        else:
            return I

    def set(self,shallow=False,**fields):
        """
        Element attribute-setter
        
        Parameters
        ----------
        fields : dict
            attribute-name:value pairs used to set the attributes
            of the calling instance
        """
        for item in self:
            item.set(shallow=shallow,**fields)

    def prune(self,objlist=[],attr_maps=[]):
        """
        Attribute-based pruning by referencing and mapping another list

        Parameters
        ----------
        objlist : list
            "reference" list of objects whose attributes are consulted to 
            decide what to prune out of the calling instance
        attr_maps : list
            dictionaries used independently to map attribute values of 
            elements of the reference list to those of the calling instance
        
        Returns
        -------

        list :
            items removed from calling instance
        """
        acc_list=self.__class__([])
        for ref_obj in objlist:
            for m in attr_maps:
                # each map is *independently* used to gather hits
                thru_dict={}
                for myattrlabel,objattrlabel in m.items():
                    thru_dict[myattrlabel]=ref_obj.__dict__[objattrlabel]
                acc_list.extend(self.filter(**thru_dict))
        for item in acc_list:
            self.remove(item)
        return acc_list
    
    def prune_exclusions(self,**kwargs):
        """
        Attribute-based list pruning

        Parameters
        ----------
        kwargs : dict
            dictionary of attribute-name:value pairs that define elements
            to be pruned from the calling instance
        
        Returns
        -------

        list :
            items removed from calling instance
        """
        acc_list=self.__class__([])
        for k,v in kwargs.items():
            for item in v:
                thru_dict={k:item}
                acc_list.extend(self.filter(**thru_dict))
        logger.debug(f'pruning out {len(acc_list)} items')
        for item in acc_list:
            if item in self:
                self.remove(item)
        return acc_list
    
    def get_attr(self,S:tuple,**fields):
        """
        Return list of S-matching attributes from among elements matching fields
         
        Parameters
        ----------
        S : tuple
            the name of an object-list attribute and dictionary of attribute-name:value
            pairs used to filter that object-list attribute
        fields : dict
            attribute-name:value pairs used to find attributes
            of the calling instance to which the S-filter is applied

        Returns
        -------
        obj, list :
            result of applying get() to object-list attributes
        """
        for r in self:
            if r.matches(**fields): break
        else:
            return self.__class__([])
        attr,subfields=S
        assert type(attr)==str
        assert type(subfields)==dict
        if not attr in r.__dict__:
            return self.__class__([])
        L=r.__dict__[attr]
        return L.get(**subfields)
    
    def sort(self,by=None,reverse=False):
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
            key=operator.attrgetter(*by)
            self.data.sort(key=key,reverse=reverse)

    def uniqattrs(self,attrs=[],with_counts=False):
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
            if with_counts ins true, list of (value,count)
        """
        uattrs={k:[] for k in attrs}
        for item in self:
            for a in attrs:
                v=getattr(item,a)
                if with_counts:
                    try:
                        idx=[x[0] for x in uattrs[a]].index(v)
                    except:
                        uattrs[a].append([v,0])
                        idx=-1
                    uattrs[a][idx][1]+=1
                else:
                    if not v in uattrs[a]:
                        uattrs[a].append(v)
        return uattrs
    
    def binnify(self,fields=[]):
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
        bins={}
        for item in self:
            key=item.strhash(fields=fields)
            if not key in bins:
                bins[key]=[]
            bins[key].append(item)
        for k,v in bins.items():
            if len(v)>1:
                logger.debug(f'binnify: bin {k} has {len(v)} items: {v}')
        return bins
    
    def puniq(self,fields=[]):
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
        bins=self.binnify(fields=fields)
        return len(bins)==len(self)
    
    def puniquify(self,fields=[],new_attr_name='_ORIGINAL_',make_common=[]):
        """
        Systematic attribute altering to make all elements unique
        
        There may be a set of attributes for which no two elements may
        have the exact same set of respective values.  This method 
        scans the calling instance for such collisions and, if any
        is found, it adds one to the value of the attributed named
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
        new_attr_name : str, optional
            name given to a new attribute used to store all original 
            attribute name:value pairs
        make_common : list
            list of attribute names that are set to a single set
            of common values *after* uniquifying

        """
        assert not any([x in fields for x in make_common])
        bins=self.binnify(fields=fields)
        tmp_fields=fields
        if not fields:
            tmp_fields=[x for x in __dict__.keys()]
        stillworking=True
        while stillworking:
            stillworking=False
            for k,v in bins.items():
                if len(v)>1:
                    stillworking=True
                    for d in v[1:]:
                        if not new_attr_name in d.__dict__:
                            d.__dict__[new_attr_name]={k:d.__dict__[k] for k in tmp_fields}
                        while d.strhash(tmp_fields) in bins:
                            d.__dict__[tmp_fields[0]]+=1
            if stillworking: 
                bins=self.binnify(fields=tmp_fields)
        assert(self.puniq(fields=tmp_fields))
        use_common={k:self[0].__dict__[k] for k in make_common}
        if use_common:
            for s in self:
                # check to see if any updates are needed
                valcheck=all([s.__dict__[k]==v for k,v in use_common.items()])
                if not valcheck:
                    if not new_attr_name in s.__dict__:
                        s.__dict__[new_attr_name]={}
                    s.__dict__[new_attr_name].update({k:s.__dict__[k] for k in make_common})
                    s.__dict__.update(use_common)
    
    def state_bounds(self,state_func):
        """
        Reduces calling instance to a list of state intervals 
        
        Parameters
        ----------
        state_func : method
            A function that returns the state value of a single 
            element of the calling instance

        Returns
        -------
        list : state intervals
        
        """
        slices=StateIntervalList([])
        if len(self)==0:
            return slices
        for i,item in enumerate(self):
            if not slices:
                slices.append(StateInterval({'state':state_func(item),'bounds':[i,i]}))
                continue
            state=state_func(item)
            last_state=slices[-1].state
            if state==last_state:
                slices[-1].increment_rightbound()
            else:
                slices.append(StateInterval({'state':state_func(item),'bounds':[i,i]}))
        return slices
    
    def map_attr(self,mapped_attr,key_attr,map):
        """
        Simple cross-attribute mapper, applied to each element
        of the calling instance
        
        Parameters
        ----------
        mapped_attr : str
            name of attribute to which the mapping result 
            will be applied
        key_attr : str
            name of attribute whose value is mapped
        map : dict
            the map
        """
        if map:
            for item in self:
                item.map_attr(mapped_attr,key_attr,map)
    
    def swap_attr(self,attr1,attr2):
        """
        Simple attribute value swapper, applied to each element
        of caller
        
        Parameters
        ----------
        attr1 : str
            name of first attribute in pair to swap
        attr2 : str
            name of second attribute in pair to swap
        """        
        for item in self.data:
            item.swap_attr(attr1,attr2)

    def copy_attr(self,recv_attr,src_attr):
        """
        Simple attribute copier, applied to each element of caller
        
        Parameters
        ----------
        recv_attr : str
            name of receiver attribute; must exist
        src_attr : str
            name of source attribute; must exist
        """
        for item in self.data:
            item.copy_attr(recv_attr,src_attr)

    def assign_objs_to_attr(self,attr,objList,**matchattr):
        """
        Assigns the single object from objList whose
        attributes match the matchattr dict to the 
        attr attribute of every element in the calling
        instance; any elements that 
         
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
        for s in self:
            s.assign_obj_to_attr(attr,objList,**matchattr)
        delete_us=[s for s in self if getattr(s,attr)==None]
        for s in delete_us:
            self.remove(s)
        return self.__class__(delete_us)
    
    def update_attr_from_obj_attr(self,attr,obj,obj_attr):
        """
        Set value of attribues of all elements of caller
        from another attribute of a separate object

        Parameters
        ----------
        attr : str
            attribute name
        obj : object
            object from which the attribute value is taken
        obj_attr : str
            name of attribute in obj that is set to 
            attr of caller
        """
        # logger.debug(f'update_attr_from_obj_attr considers type {type(self)} ({len(self)}) attr {attr} obj {obj} obj_attr {obj_attr}')
        for item in self:
            # logger.debug(f'-> item {item}')
            item.update_attr_from_obj_attr(attr,obj,obj_attr)

    def update_attr_from_objlist_elem_attr(self,attr,objlist,index,obj_attr):
        """
        Set value of attribues of all elements of caller
        from another attribute of index'th element objectlist

        Parameters
        ----------
        attr : str
            attribute name
        objlist : str
            name of list of objects from which the attribute value of the 
            index'th element is taken
        index : int
            index of element in objlist from which attribute value is taken
        obj_attr : str
            name of attribute in objlist[index] that is set to 
            attr of caller
        """
        for item in self: 
            item.update_attr_from_objlist_elem_attr(attr,objlist,index,obj_attr)

    def remove_duplicates(self,fields=[]):
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

class CloneableObjList(ObjList):
    """
    A class for lists of cloneable objs 
    
    """
    def clone(self,**options):
        """
        List cloner

        Parameters
        ----------
        options : dict
            dictionary of attribute-name:value pairs that are 
            explicitly assigned in the clones
        
        Returns
        -------
        list :
            The new list of new cloned objects   
        """
        R=self.__class__([])
        for item in self:
            R.append(item.clone(**options))
        return R

class AncestorAwareObjList(CloneableObjList):
    """
    A class for lists of ancestor-aware objs 

    """
    def claim_descendants(self,stamp):
        """
        Instructs caller to stamp all elements' descendant objects 
        
        Parameters
        ----------
        stamp : obj
            object assigned to all ancestor_obj attributes of any
            ancestor-aware attributes
        """
        for obj in self:
            obj.claim_descendants(stamp)

class StateIntervalList(AncestorAwareObjList):
    """
    A class for lists of StateIntervals
    """

    def insert(self,alien:StateInterval):
        """
        Inserts a StateInterval into the calling instance
        The insertion is done in such a way that the calling instance
        remains sorted by the bounds of the intervals, and that
        the new interval is inserted in such a way that it does not
        overlap with any existing intervals.  If the new interval
        is completely contained within an existing interval, it is
        not inserted.  If the new interval overlaps with an existing
        interval, the existing interval is split into two intervals
        to accommodate the new interval.

        Parameters
        ----------
        alien : StateInterval
            the StateInterval to be inserted into the calling instance

        """
        for cont in self:
            if alien.bounds[0]>cont.bounds[0] and alien.bounds[1]<cont.bounds[1]:
                break
        else:
            cont=None
        if cont:
            i=self.index(cont)+1
            self.insert(i,alien)
            self.insert(i+1,StateInterval({'bounds':[alien.bounds[1]+1,cont.bounds[1]],'state':cont.state,'build':cont.build}))
            cont.bounds[1]=alien.bounds[0]-1
        else:
            logger.debug(f'No compatible insertion point for alien interval {str(alien)}')