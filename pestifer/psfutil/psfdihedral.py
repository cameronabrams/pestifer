# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
# Definition of the :class:`PSFDihedral` and :class:`PSFDihedralList` classes for handling dihedrals in PSF topology files.
# These classes are used to represent dihedrals formed by four atoms in a molecular structure, as defined in PSF files.
# The :class:`PSFDihedral` class represents a single dihedral, while the :class:`PSFDihedralList` class is a collection of such dihedrals.
# These are descendants of the :class:`PSFTopoElement <.psftopoelement.PSFTopoElement>` and :class:`PSFTopoElementList <.psftopoelement.PSFTopoElementList>` classes, respectively,
# which provide a framework for representing and manipulating generic PSF topology elements.
"""

import logging

from functools import singledispatchmethod

from .psftopoelement import PSFTopoElement, PSFTopoElementList, LineList

logger=logging.getLogger(__name__)

class PSFDihedral(PSFTopoElement):
    """
    A class representing a dihedral in a PSF topology file.
    """
    
    def __eq__(self,other):
        """
        Check if two PSFDihedral objects are equal by comparing their serial numbers.
        Two dihedrals are considered equal if they have the same serial numbers in either forward or reverse order.
        
        Parameters
        ----------
        other : PSFDihedral
            The other dihedral to compare against.
        
        Returns
        -------
        bool
            True if the dihedrals are equal, False otherwise.
        """
        return [self.serial1,self.serial2,self.serial3,self.serial4]==[other.serial1,other.serial2,other.serial3,other.serial4] or [self.serial1,self.serial2,self.serial3,self.serial4]==[other.serial4,other.serial3,other.serial2,other.serial1] 

class PSFDihedralList(PSFTopoElementList):
    """
    A class representing a list of :class:`PSFDihedral` objects.  This class can be initialized from a list of lines or from a :class:`PSFTopoElementList <.psftopoelement.PSFTopoElementList>` object.
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
            assert len(tokens)%4==0,f'Poorly formatted PHI line in psffile?'
            if not include_serials:
                for li,i,j,rj in zip(tokens[:-3:4],tokens[1:-2:4],tokens[2:-1:4],tokens[3::4]):
                    B.append(PSFDihedral([li,i,j,rj]))
            else:
                for li,i,j,rj in zip(tokens[:-3:4],tokens[1:-2:4],tokens[2:-1:4],tokens[3::4]):
                    if include_serials[li-1] and include_serials[i-1] and include_serials[j-1] and include_serials[rj-1]:
                        B.append(PSFDihedral([li,i,j,rj]))
        super().__init__(B)
