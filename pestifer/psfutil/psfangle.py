# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`PSFAngle` and :class:`PSFAngleList` classes for handling angles in PSF topology files.
These classes are used to represent angles formed by three atoms in a molecular structure, as defined in PSF files.
The :class:`PSFAngle` class represents a single angle, while the :class:`PSFAngleList` class is a collection of such angles.
These are descendants of the :class:`PSFTopoElement <.psftopoelement.PSFTopoElement>` and :class:`PSFTopoElementList <.psftopoelement.PSFTopoElementList>` classes, respectively, which provide a framework for representing and manipulating generic PSF topology elements.
"""

from functools import singledispatchmethod

from .psftopoelement import PSFTopoElement, PSFTopoElementList, LineList
from .psfbond import PSFBond

class PSFAngle(PSFTopoElement):
    """
    A class representing an angle in a PSF topology file.
    An angle is defined by three atom serial numbers, which represent the atoms forming the angle.
    """
    
    def __eq__(self,other):
        """
        Check if two PSFAngle objects are equal by comparing their serial numbers.
        Two angles are considered equal if they have the same serial numbers in either forward or reverse order.
        """
        return [self.serial1,self.serial2,self.serial3]==[other.serial1,other.serial2,other.serial3] or [self.serial1,self.serial2,self.serial3]==[other.serial3,other.serial2,other.serial1] 

class PSFAngleList(PSFTopoElementList):
    """
    A class representing a list of :class:`PSFAngle` objects.
    This class inherits from :class:`PSFTopoElementList <.psftopoelement.PSFTopoElementList>` and provides methods for managing a collection of angles in a PSF topology file.
    It can be initialized from a list of lines or from a :class:`PSFTopoElementList <.psftopoelement.PSFTopoElementList>` object.
    """
    
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
