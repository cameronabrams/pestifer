"""

.. module:: basemod
   :synopsis: defines the BaseMod class and ModList class
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import operator
import yaml
import logging
from collections import UserList
logger=logging.getLogger(__name__)
class BaseMod:
    # list of required attribute labels
    req_attr=[]
    # list of optional attribute labels
    opt_attr=[]
    # list of 2-tuples of attribute labels for mutually exclusive attributes
    alt_attr=[]
    # dictionary of lists of value choices for any attribute, keyed on attribute label
    attr_choices={}
    # dictionary of lists of labels of attributes that are required if the optional attributes whose labels key this dictionary are present 
    opt_attr_deps={}
    # list of attributes that are IGNORED for all comparisons
    ignore_attr=[]
    def __init__(self,input_dict={}):
        # print(f'b {input_dict}')
        # mutually exclusive attribute labels should not appear in the list of required attributes
        assert all([(not x in self.req_attr and not y in self.req_attr) for x,y in self.alt_attr]),"Mutually exclusive attributes should not appear in the required list"
        for k,v in input_dict.items():
            if k in self.req_attr or k in self.opt_attr:
                self.__dict__[k]=v
        # all required attributes have values
        assert all([att in self.__dict__.keys() for att in self.req_attr]),f"Not all required attributes have values\nRequired: {self.req_attr}\nProvided: {list(self.__dict__.keys())}"
        # all mutual-exclusivity requirements are met
        assert all([
            ((x[0] in self.__dict__.keys() and not x[1] in self.__dict__.keys()) or 
             (x[1] in self.__dict__.keys() and not x[0] in self.__dict__.keys()) or not (x[0] in self.__dict__.keys() or x[1] in self.__dict__.keys()))
            for x in self.alt_attr]),"Mutual exclusivity requirements unmet"
        # all attributes with limited set of possible values have valid values
        assert all([self.__dict__[x] in choices for x,choices in self.attr_choices.items()]),"Invalid choices in one or more attributes"
        # for each optional attribute present, all required dependent attributes are also present
        for oa,dp in self.opt_attr_deps.items():
            if oa in self.__dict__.keys():
                assert all([x in self.__dict__.keys() for x in dp]),"Dependent required attributes of optional attributes not present"

    def __eq__(self,other):
        """__eq__ == operator for BaseMods.  A and B are considered equal if they have the same values for all required attributes and for any *common* optional attributes.
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
        # logger.debug(f'Equality test for class {str(self.__class__)}: attributes considered are {attr_list}')
        test_list=[]
        for k in attr_list:
            # if type(self.__dict__[k])==str or not hasattr(self.__dict__[k],'__setitem__'):
            test_list.append(self.__dict__[k]==other.__dict__[k])
            # else:
            #     logger.debug(f'Ignoring container-like attribute {k} in comparison of two {str(self.__class__)} instances')
        return all(test_list)
    def __lt__(self,other):
        """__lt__ < operator for BaseMods.  A<B if all required and *common* optional attributes of A are less than or equal to those of B AND *at least one* such attribute of A is *strictly less than* its counterpart in B.
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
        """ __le__ <= operator for BaseMods """
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
    # def claim_object(self,obj,stamp):
    #     if hasattr(obj,'__dict__'):
    #         obj.__dict__['ancestor_obj']=stamp
    #         obj.claim_descendants(stamp)
    #     elif hasattr(obj,'len'): # a dict with AncestorAwareMods in the values or a list of AncestorAwareMods
    #         if type(obj)==dict:
    #             for subattr,subobj in obj.__dict__.items():
    #                 self.claim_object(subobj,stamp)
    #         else: #assume list
    #             for subobj in obj:
    #                 self.claim_object(subobj,stamp)
    
    def claim_self(self,stamp,level):
        self.__dict__['ancestor_obj']=stamp

    def claim_descendants(self,stamp,level):
        self.claim_self(stamp,level+1)
        for attr,obj in self.__dict__.items():
            if attr!='ancestor_obj':
                if hasattr(obj,'__dict__'):
                    obj.claim_descendants(stamp,level+1)
                    

class StateInterval(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['state','bounds']
    def increment_rightbound(self):
        self.bounds[1]+=1
    def pstr(self):
        return f'{self.state}({self.bounds[1]-self.bounds[0]+1})'

class ModList(UserList):
    def __init__(self,data):
        super().__init__(data)
    # def __iter__(self):
    #     for r in self.L:
    #         yield r
    # def __getitem__(self,index):
    #     return self.L[index]
    # def __setitem__(self,index,value):
    #     self.L[index]=value
    # def __add__(self,other):
    #     self.L.extend(other.L)
    #     return self
    # def __delitem__(self,index):
    #     del self.L[index]
    # def __contains__(self,item):
    #     return item in self.L
    # def __eq__(self,other):
    #     return all([x==y for x,y in zip(self.L,other.L)])
    # def __len__(self):
    #     return len(self.L)
    # def index(self,item):
    #     for i,x in enumerate(self.L):
    #         if item==x:
    #             return i
    #     raise Exception(ValueError,f'Mod not found in list')
    # def append(self,item):
    #     self.L.append(item)
    # def extend(self,items):
    #     self.L.extend(items)
    # def remove(self,item):
    #     self.L.remove(item)
    # def clear(self):
    #     self.L.clear()
    # def pop(self,idx):
    #     item=self.L.pop(idx)
    #     return item
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
        # logger.debug(f'pruning out {len(acc_list)} items')
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


class CloneableModList(ModList):
    def clone(self,**options):
        R=self.__class__([])
        for item in self:
            R.append(item.clone(**options))
        return R

class AncestorAwareModList(CloneableModList):
    def claim_descendants(self,stamp,level):
        for obj in self:
            if hasattr(obj,'__dict__'):
                obj.__dict__['ancestor_obj']=stamp
                obj.claim_descendants(stamp,level)

class StateIntervalList(AncestorAwareModList):
    pass
