# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import networkx as nx

from functools import singledispatchmethod

from .psftopoelement import PSFTopoElement, PSFTopoElementList, LineList


class PSFBond(PSFTopoElement):

    def __eq__(self,other):
        return [self.serial1,self.serial2]==[other.serial1,other.serial2] or [self.serial1,self.serial2]==[other.serial2,other.serial1]
    
class PSFBondList(PSFTopoElementList):
    @singledispatchmethod
    def __init__(self,input_data,**kwargs):
        super().__init__(input_data)

    @__init__.register(PSFTopoElementList)
    def _from_tel(self,input_obj,**kwargs):
        super().__init__(input_obj)

    @__init__.register(LineList)
    def _from_list(self,lines,include_serials=[]):
        B=[]
        for line in lines:
            tokens=[int(x) for x in line.split()]
            assert len(tokens)%2==0,f'Poorly formatted BOND line in psffile?'
            if not include_serials:
                for l,r in zip(tokens[:-1:2],tokens[1::2]):
                    B.append(PSFBond([l,r]))
            else:
                for l,r in zip(tokens[:-1:2],tokens[1::2]):
                    if include_serials[l-1] and include_serials[r-1]:
                        B.append(PSFBond([l,r]))
        super().__init__(B)

    def to_graph(self):
        g=nx.Graph()
        for b in self:
            g.add_edge(b.serial1,b.serial2)
        return g
