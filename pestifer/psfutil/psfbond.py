# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Definition of the :class:`PSFBond` and :class:`PSFBondList` classes for handling bonds in PSF topology files.
These classes are used to represent bonds formed between pairs of atoms in a molecular structure, as defined in PSF files.
The :class:`PSFBond` class represents a single bond, while the :class:`PSFBondList` class is a collection of such bonds.
These are descendants of the :class:`PSFTopoElement <.psftopoelement.PSFTopoElement>` and :class:`PSFTopoElementList <.psftopoelement.PSFTopoElementList>` classes, respectively,
which provide a framework for representing and manipulating generic PSF topology elements"""

import networkx as nx

from functools import singledispatchmethod

from .psftopoelement import PSFTopoElement, PSFTopoElementList, LineList


class PSFBond(PSFTopoElement):
    """
    A class representing a bond in a PSF topology file.
    """

    def __eq__(self,other):
        """
        Check if two PSFBond objects are equal by comparing their serial numbers.
        Two bonds are considered equal if they have the same serial numbers in either forward or reverse order.
        
        Parameters
        ----------
        other : PSFBond
            The other bond to compare against.
        """
        return [self.serial1,self.serial2]==[other.serial1,other.serial2] or [self.serial1,self.serial2]==[other.serial2,other.serial1]
    
class PSFBondList(PSFTopoElementList):
    """
    A class representing a list of :class:`PSFBond` objects.
    This class inherits from :class:`PSFTopoElementList <.psftopoelement.PSFTopoElementList>` and provides methods for managing a collection of bonds in a PSF topology file.
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
        """
        Create a :class:`networkx.Graph` representation of the bond list, where bonds are edges between atom serial numbers.

        Returns
        -------
        networkx.Graph
            A networkx graph representing the bond list, with atom serial numbers as nodes and bonds as edges.
        """
        g=nx.Graph()
        for b in self:
            g.add_edge(b.serial1,b.serial2)
        return g
