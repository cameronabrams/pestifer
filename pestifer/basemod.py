# Author: Cameron F. Abrams <cfa22@drexel.edu>.
"""Simple namespace objects with attribute controls

This module defines the BaseMod class as a derivative of argparse's
very nice Namespace class, and the ModList class as a derivative
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

"""
import operator
import yaml
import logging
import inspect
from collections import UserList
logger=logging.getLogger(__name__)
from argparse import Namespace

class BaseMod(Namespace):
    """A class defining a namespace with custom attribute controls.

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
    
    ...
    
    Attributes
    ----------
    req_attr : list
        labels of required attributes
    opt_attr : list
        labels of optional attributes
    alt_attr : list
        pairs of attribute labels in which each pair declares the 
        two attributes to be mutually exclusive
    attr_choices : dict
        lists of value options for any attribute;
        e.g., {'response':['yes','no','maybe'], ...}
    opt_attr_deps : dict
        lists of attributes that are required by optional attributes;
        e.g., {'flavor':['is_edible','is_flavorable'], ...}
    ignore_att : list
        attributes whose values are ignored when comparing instances


    Methods
    -------

    __eq__(other)
        equality operator for BaseMod instances.  Two instances are considered
        equal if
        * they are actually the same instance; or 
        * they have the same values for all required attributes and for any 
          common optional attributes.
    
    __lt__(other)
       less than operator for BaseMod instances.  Instance A is less than B if
       * all required and common optional attributes of A are less than or equal 
         to those of B; AND
       * at least one such attribute of A is strictly less than that of B  

    """
    req_attr=[]
    opt_attr=[]
    alt_attr=[]
    attr_choices={}
    opt_attr_deps={}
    ignore_attr=[]
    def __init__(self,input_dict={}):
        """
        Parameters
        ----------
        input_dict : dict, optional
            dictionary of input attribute:value pairs
        """
        # mutually exclusive attribute labels should not appear in the list of required attributes
        assert all([(not x in self.req_attr and not y in self.req_attr) for x,y in self.alt_attr]),"Mutually exclusive attributes should not appear in the required list"
        prep_dict={}
        for k,v in input_dict.items():
            if k in self.req_attr or k in self.opt_attr:
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
            raise Exception(f'Trouble initializing instance of {self.__class__}')
        # for each optional attribute present, all required dependent attributes are also present
        for oa,dp in self.opt_attr_deps.items():
            if oa in prep_dict.keys():
                assert all([x in prep_dict.keys() for x in dp]),"Dependent required attributes of optional attributes not present"
        super().__init__(**prep_dict)

    def __eq__(self,other):
        """
        Parameters
        ----------
        other : BaseMod
            other BaseMod instance
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
        Parameters
        ----------
        other : BaseMod
            other BaseMod instance
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
        Parameters
        ----------
        other : BaseMod
            other BaseMod instance
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
        Parameters
        ----------
        other : BaseMod
            other BaseMod instance
        attr : list, optional
            attributes used for the less-than comparison
        """
        lt_list=[self.__dict__[k]<other.__dict__[k] for k in attr if hasattr(self.__dict__[k],'__lt__')]
        le_list=[self.__dict__[k]<=other.__dict__[k] for k in attr if hasattr(self.__dict__[k],'__lt__')]
        return all(le_list) and any(lt_list)
    def __pkey(self,fields=[]):
        if not fields:
            return tuple([self.__dict__.values()])
        return tuple([x for k,x in self.__dict__.items() if k in fields])
    def phash(self,fields=[]):
        return hash(self.__pkey(fields))
    def matches(self,**fields):
        # non strict match -- returns true if all key:val pairs in fields
        # match the corresponding key:val pairs in self.__dict__
        for k,v in fields.items():
            if not k in self.__dict__:
                # logger.debug(f'missing key {k}')
                return False
            if v!=self.__dict__[k]:
                # logger.debug(f'wrong value ({self.__dict__[k]}) for key {k}; wanted {v}')
                return False
        return True
    def dump(self):
        retdict={}
        retdict['instanceOf']=type(self).__name__
        retdict.update(self.__dict__)
        return yaml.dump(retdict)
    def inlist(self,a_list):
        for s in a_list:
            if s==self:
                return True
        return False
    def map_attr(self,mapped_attr,key_attr,map):
        if map:
            key=self.__dict__[key_attr]
            val=map[key]
            self.__dict__[mapped_attr]=val
    def set(self,**fields):
        for k,v in fields.items():
            if k in self.__dict__:
                self.__dict__[k]=v
    def assign_obj_to_attr(self,attr,objList,**matchattr):
        assert getattr(self,attr)==None
        adict={k:getattr(self,v) for k,v in matchattr.items()}
        myObj=objList.get(**adict)
        if myObj!=None:
            setattr(self,attr,myObj)
    def update_attr_from_obj_attr(self,attr,obj,obj_attr):
        setattr(self,attr,getattr(getattr(self,obj),obj_attr))

class CloneableMod(BaseMod):
    opt_attr=BaseMod.opt_attr+['clone_of']
    def clone(self,**options):
        input_dict={k:v for k,v in self.__dict__.items() if k in self.req_attr}
        input_dict.update({k:v for k,v in self.__dict__.items() if k in self.opt_attr})
        input_dict.update(options)
        input_dict['clone_of']=self
        x=self.__class__(input_dict)
        return x
    def is_clone(self):
        return 'clone_of' in self.__dict__
    def get_original(self):
        return self.clone_of
    
class AncestorAwareMod(CloneableMod):
    opt_attr=CloneableMod.req_attr+['ancestor_obj']
    ignore_attr=CloneableMod.ignore_attr+['ancestor_obj']
    
    def claim_self(self,stamp):
        self.__dict__['ancestor_obj']=stamp

    def claim_descendants(self,stamp):
        self.claim_self(stamp)
        for attr_name,obj in self.__dict__.items():
            if not attr_name=='ancestor_obj': # very important! or, could stamp self last?
                classes=inspect.getmro(type(obj))
                if AncestorAwareMod in classes or AncestorAwareModList in classes:
                    obj.claim_descendants(stamp)

class StateInterval(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['state','bounds']
    def increment_rightbound(self):
        self.bounds[1]+=1
    def num_residues(self):
        return self.bounds[1]-self.bounds[0]
    def pstr(self):
        return f'{self.state}({self.bounds[1]-self.bounds[0]+1})'

class ModList(UserList):
    def __init__(self,data):
        super().__init__(data)

    def filter(self,**fields):
        retlist=self.__class__([])
        for r in self:
            if r.matches(**fields): retlist.append(r)
        return retlist
    def get(self,**fields):
        R=self.filter(**fields)
        if len(R)==0:
            return None
        elif len(R)==1:
            return R[0]
        else:
            return R
    def set(self,**fields):
        for item in self:
            item.set(**fields)

    def prune(self,objlist=[],attr_maps=[]):
        """ given a list of objects and a list of mappings of my element attributes to 
            the foreign object attributes, delete any of my elements whose attribute
            values match the mapped foreign object attributes """
        acc_list=self.__class__([])
        for obj in objlist:
            for m in attr_maps:
                thru_dict={}
                for myattrlabel,objattrlabel in m.items():
                    thru_dict[myattrlabel]=obj.__dict__[objattrlabel]
                acc_list.extend(self.filter(**thru_dict))
        for item in acc_list:
            if item in self:
                self.remove(item)
    def prune_exclusions(self,**kwargs):
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
        if not by:
            self.data.sort(reverse=reverse)
        else:
            key=operator.attrgetter(*by)
            self.data.sort(key=key,reverse=reverse)
    def uniqattrs(self,attrs=[],with_counts=False):
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
        bins={}
        for item in self:
            key=item.phash(fields)
            if not key in bins:
                bins[key]=[]
            bins[key].append(item)
        return bins
    def puniq(self,fields):
        bins=self.binnify(fields)
        return len(bins)==len(self)
    def puniquify(self,fields,new_attr_name='_ORIGINAL_',make_common=[]):
        # each element must be unique in the sense that no two 
        # elements have the same values of attributes listed in 
        # fields[]; returns a map keyed by unique index whose
        # values are the original values of the attributes listed
        # in fields; uniqueness is achieved by changing the 
        # value of the leading field only
        assert not any([x in fields for x in make_common])
        bins=self.binnify(fields)
        stillworking=True
        while stillworking:
            stillworking=False
            for k,v in bins.items():
                if len(v)>1:
                    stillworking=True
                    for d in v[1:]:
                        if not new_attr_name in d.__dict__:
                            d.__dict__[new_attr_name]={k:d.__dict__[k] for k in fields}
                        while d.phash(fields) in bins:
                            d.__dict__[fields[0]]+=1
            if stillworking: 
                bins=self.binnify(fields)
        assert(self.puniq(fields))
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
        if map:
            for item in self:
                item.map_attr(mapped_attr,key_attr,map)
    def assign_objs_to_attr(self,attr,objList,**matchattr):
        for s in self:
            s.assign_obj_to_attr(attr,objList,**matchattr)
        delete_us=[s for s in self if getattr(s,attr)==None]
        for s in delete_us:
            self.remove(s)
        return self.__class__(delete_us)
    def update_attr_from_obj_attr(self,attr,obj,obj_attr):
        for item in self:
            item.update_attr_from_obj_attr(attr,obj,obj_attr)

class CloneableModList(ModList):
    def clone(self,**options):
        R=self.__class__([])
        for item in self:
            R.append(item.clone(**options))
        return R

class AncestorAwareModList(CloneableModList):
    def claim_descendants(self,stamp):
        for obj in self:
            obj.claim_descendants(stamp)

class StateIntervalList(AncestorAwareModList):
    pass
