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
    req_attr={}
    opt_attr={}
    alt_attr=[]
    def __init__(self,input_dict):
        assert all([att in input_dict for att in self.req_attr])
        assert all([(x[0] in input_dict or x[1] in input_dict) for x in self.alt_attr])
        self.__dict__.update(input_dict)
        self.td=input_dict
    def dump(self):
        return yaml.dump(self.td)
    def psfgen_lines(self):
        lines=[]
        return lines
    # def clone_to_chain(self,chainID):
    #     new=BaseMod(self.td)
    #     new.chainID=chainID
    #     return new
    
class Mutation(BaseMod):
    req_attr={'chainID':'','origresname':'','newresname':''}
    opt_attr={'resseqnum':0,'resseqnumi':'','insertion':''}
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
    
    def clone(self,newchainID):
        newmut=Mutation(self.td)
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
