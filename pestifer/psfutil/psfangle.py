# Author: Cameron F. Abrams, <cfa22@drexel.edu>

from functools import singledispatchmethod

from .psftopoelement import PSFTopoElement, PSFTopoElementList, LineList
from .psfbond import PSFBond

class PSFAngle(PSFTopoElement):

    def __eq__(self,other):
        return [self.serial1,self.serial2,self.serial3]==[other.serial1,other.serial2,other.serial3] or [self.serial1,self.serial2,self.serial3]==[other.serial3,other.serial2,other.serial1] 

class PSFAngleList(PSFTopoElementList):
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
            assert len(tokens)%3==0,f'Poorly formatted ANGLE line in psffile?'
            if not include_serials:
                for l,m,r in zip(tokens[:-2:3],tokens[1:-1:3],tokens[2::3]):
                    B.append(PSFAngle([l,m,r]))
            else:
                for l,m,r in zip(tokens[:-2:3],tokens[1:-1:3],tokens[2::3]):
                    if include_serials[l-1] and include_serials[m-1] and include_serials[r-1]:
                        B.append(PSFBond([l,m,r]))
        super().__init__(B)
