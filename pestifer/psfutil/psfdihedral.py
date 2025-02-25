# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging

from functools import singledispatchmethod

from .psftopoelement import PSFTopoElement, PSFTopoElementList, LineList

logger=logging.getLogger(__name__)

class PSFDihedral(PSFTopoElement):
    def __eq__(self,other):
        return [self.serial1,self.serial2,self.serial3,self.serial4]==[other.serial1,other.serial2,other.serial3,other.serial4] or [self.serial1,self.serial2,self.serial3,self.serial4]==[other.serial4,other.serial3,other.serial2,other.serial1] 

class PSFDihedralList(PSFTopoElementList):
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
            assert len(tokens)%4==0,f'Poorly formatted PHI line in psffile?'
            if not include_serials:
                for li,i,j,rj in zip(tokens[:-3:4],tokens[1:-2:4],tokens[2:-1:4],tokens[3::4]):
                    B.append(PSFDihedral([li,i,j,rj]))
            else:
                for li,i,j,rj in zip(tokens[:-3:4],tokens[1:-2:4],tokens[2:-1:4],tokens[3::4]):
                    if include_serials[li-1] and include_serials[i-1] and include_serials[j-1] and include_serials[rj-1]:
                        B.append(PSFDihedral([li,i,j,rj]))
        super().__init__(B)
