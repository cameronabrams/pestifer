"""

.. module:: mods
   :synopsis: Manages all modifications
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
# from pestifer.mutation import Mutation
# from pestifer.ssbond import SSBond
# from pestifer.graft import Graft
# from pestifer.crot import Crot
# from pestifer.attach import Attach
# from pestifer.link import Link
# from pestifer.deletion import Deletion
# from pestifer.cleavage import Cleavage
# from pestifer.missing import Missing
import yaml
from .seqadv import Seqadv
class BaseMod:
    req_attr=[]
    opt_attr=[]
    alt_attr=[]
    def __init__(self,input_dict):
        if 'shortcode' in input_dict:
            self.parse_shortcode(input_dict)
        else:
            for k,v in input_dict.items():
                if k in self.req_attr or k in self.opt_attr:
                    self.__dict__[k]=v
        assert all([att in self.__dict__.keys() for att in self.req_attr])
        assert all([(x[0] in self.__dict__.keys() or x[1] in self.__dict__.keys()) for x in self.alt_attr])
    def __eq__(self,other):
        test_list=[self.__dict__[k]==other.__dict__[k] for k in self.req_attr]
        for k in self.opt_attr:
            if k in self.__dict__ and k in other.__dict__:
                test_list.append(self.__dict__[k]==other.__dict__[k])
        return all(test_list)
    def dump(self):
        return yaml.dump(self.td)
    def psfgen_lines(self):
        lines=[]
        return lines

class Mutation(BaseMod):
    req_attr=['chainID','origresname','newresname']
    opt_attr=['resseqnum','resseqnumi','insertion']
    alt_attr=[('resseqnum','insertion')]
    @classmethod
    def from_seqdav(cls,sq:Seqadv):
        input_dict={}
        input_dict['chainID']=sq.chainID
        input_dict['origresname']=sq.resName
        input_dict['newresname']=sq.dbRes
        input_dict['resseqnum']=sq.seqNum
        input_dict['insertion']=sq.iCode
        input_dict['resseqnumi']=f'{sq.seqNum}{sq.iCode}'
        inst=cls(input_dict)
        return inst
    def parse_shortcode(self,input_dict):
        ### c:nnn,rrr,nnn
        shortcode=input_dict['shortcode']
        s1=shortcode.split(':')
        self.chainID=s1[0]
        s2=s1[1].split(',')
        self.origresname=s2[0]
        self.resseqnumi=s2[1]
        self.newresname=s2[2]
        if not self.resseqnumi[-1].isdigit():
            self.insertion=self.resseqnumi[-1]
            self.resseqnum=self.resseqnumi[:-1]
        else:
            self.insertion=' '
            self.resseqnum=self.resseqnumi
    
    def clone(self,newchainID):
        input_dict={k:v for k,v in self.__dict__.items() if k in self.req_attr}
        input_dict.update({k:v for k,v in self.__dict__.items() if k in self.opt_attr})
        newmut=Mutation(input_dict)
        newmut.chainID=newchainID
        return newmut
    
    def psfgen_lines(self):
         lines=['    mutate {self.resseqnumi} {self.newresname}']
         return lines

ModTypes={
    'Mutations':Mutation,
    # 'Grafts':Graft,
    # 'Deletions':Deletion,
    # 'Crotations':Crot,
    # 'Attachments':Attach,
    # 'Links':Link,
    # 'SSbonds':SSBond,
    # 'Cleavages':Cleavage,
    # 'Missing':Missing
}

class ModsContainer:
    def __init__(self,input_dict):
        self.M=[]
        for k,v in input_dict.items():
            if k in ModTypes:
                self.M.append(ModTypes[k](v))
    def update(self,input_dict):
        for k,v in input_dict.items():
            if k in ModTypes:
                self.M.append(ModTypes[k](v))
    def dump(self):
        retstr=''
        for m in self.M:
            retstr+=yaml.dump(m)
        return retstr
