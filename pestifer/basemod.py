"""

.. module:: basemod
   :synopsis: defines the BaseMod class and ModList class
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import operator
import yaml
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
    def __init__(self,input_dict={}):
        # mutually exclusive attribute labels should not appear in the list of required attributes
        assert all([(not x in self.req_attr and not y in self.req_attr) for x,y in self.alt_attr]),"Mutually exclusive attributes should not appear in the required list"
        for k,v in input_dict.items():
            if k in self.req_attr or k in self.opt_attr:
                self.__dict__[k]=v
        # all required attributes have values
        assert all([att in self.__dict__.keys() for att in self.req_attr]),f"Not all required attributes have values {self.__dict__}"
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
        if not self.__class__==other.__class__:
            return False
        attr_list=self.req_attr+[x for x in self.opt_attr if (x in other.__dict__ and x in self.__dict__)]
        test_list=[self.__dict__[k]==other.__dict__[k] for k in attr_list]
        return all(test_list)
    def __lt__(self,other):
        """__lt__ < operator for BaseMods.  A<B if all required and *common* optional attributes of A are less than or equal to those of B AND *at least one* such attribute of A is *strictly less than* its counterpart in B.
        """
        if not other:
            return False
        if not self.__class__==other.__class__:
            return False
        attr_list=self.req_attr+[x for x in self.opt_attr if (x in other.__dict__ and x in self.__dict__)]
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
        le_list=[self.__dict__[k]<=other.__dict__[k] for k in attr_list if hasattr(self.__dict__[k],'__lt__')]
        return all(le_list)
    def weak_lt(self,other,attr=[]):
        lt_list=[self.__dict__[k]<other.__dict__[k] for k in attr if hasattr(self.__dict__[k],'__lt__')]
        le_list=[self.__dict__[k]<=other.__dict__[k] for k in attr if hasattr(self.__dict__[k],'__lt__')]
        return all(le_list) and any(lt_list)

    def matches(self,**fields):
        # non strict match -- returns true if all key:val pairs in fields
        # match the corresponding key:val pairs in self.__dict__
        acc=True
        for k,v in fields.items():
            if not k in self.__dict__:
                return False
            acc&=(v==self.__dict__[k])
        return acc
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

class CloneableMod(BaseMod):
    def clone(self,**options):
        input_dict={k:v for k,v in self.__dict__.items() if k in self.req_attr}
        input_dict.update({k:v for k,v in self.__dict__.items() if k in self.opt_attr})
        input_dict.update(options)
        x=self.__class__(input_dict)
        return x

class ModList(list):
    def __init__(self,data):
        self.L=data
    def __iter__(self):
        for r in self.L:
            yield r
    def __getitem__(self,index):
        return self.L[index]
    def __setitem__(self,index,value):
        self.L[index]=value
    def __add__(self,other):
        self.L.extend(other.L)
        return self
    def __delitem__(self,index):
        del self.L[index]
    def __eq__(self,other):
        return all([x==y for x,y in zip(self.L,other.L)])
    def __len__(self):
        return len(self.L)
    def index(self,item):
        for i,x in enumerate(self.L):
            if item==x:
                return i
        raise Exception(ValueError,f'Mod not found in list')
    def append(self,item):
        self.L.append(item)
    def extend(self,items):
        self.L.extend(items)
    def remove(self,item):
        self.L.remove(item)
    def clear(self):
        self.L.clear()
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
    def get_attr(self,S:tuple,**fields):
        for r in self:
            if r.matches(**fields): break
        else:
            return None
        attr,subfields=S
        assert type(attr)==str
        assert type(subfields)==dict
        if not attr in r.__dict__:
            return None
        L=r.__dict__[attr]
        return L.get(**subfields)
    def sort(self,by=None,reverse=False):
        if not by:
            self.L.sort(reverse=reverse)
        else:
            key=operator.attrgetter(*by)
            self.L.sort(key=key,reverse=reverse)

class CloneableModList(ModList):
    def clone(self,**options):
        R=CloneableModList([])
        for item in self:
            R.append(item.clone(**options))
        return R
