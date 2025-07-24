# Author: Cameron F. Abrams, <cfa22@drexel.edu>

from .stateinterval import StateIntervalList, StateInterval
from ..core.baseobj_new import BaseObj, BaseObjList
from pydantic import Field
from typing import List
import logging

logger = logging.getLogger(__name__) 

class Residue(BaseObj):
    """
    A class to represent a residue in a molecule.
    This class extends StateInterval to include additional properties specific to residues.
    """
    _required_fields = ['resname', 'resseqnum', 'insertion', 'chainID']
    resname: str = Field(..., description="Name of the residue")
    resseqnum: int = Field(..., description="Sequence number of the residue")
    insertion: str = Field(..., description="Insertion code for the residue")
    chainID: str = Field(..., description="Chain ID for the residue")

    def describe(self):
        return f"Residue(resname={self.resname}, resseqnum={self.resseqnum}, insertion={self.insertion}, chainID={self.chainID})"

class ResidueList(BaseObjList):
    """
    A list of StateInterval objects representing residues.
    This class extends BaseObjList to manage a collection of StateInterval instances.
    """
    def describe(self) -> str:
        return f"ResidueList with {len(self)} residues"

    def _validate_item(self, item: Residue) -> None:
        """
        Validate that the item is an instance of StateInterval.
        
        Parameters
        ----------
        item : StateInterval
            The item to validate.
        
        Raises
        ------
        TypeError
            If the item is not an instance of StateInterval.
        """
        if not isinstance(item, Residue):
            raise TypeError(f"Item must be an instance of StateInterval, got {type(item)}")
        
    def state_bounds(self,state_func):
        """
        Reduces calling instance (a list of residues) to a list of state intervals.  A state interval is defined as a contiguous set of residues with the same state, defined by a left and right bound.  The left bound is the index of the first residue in the contiguous set, and the right bound is the index of the last residue in the contiguous set. The parameter ``state_func`` is a method that returns the state value of a single element of the calling instance.
        
        Parameters
        ----------
        state_func : method
            A function that returns the state value of a single 
            element of the calling instance

        Returns
        -------
        StateIntervalList
            A list of StateInterval objects representing the state intervals
            of the calling instance.
        
        """
        slices=StateIntervalList([])
        if len(self)==0: # empty list means empty state interval list
            return slices
        first_slice=StateInterval({'state':state_func(self[0]),'bounds':[0,0]})
        slices.append(first_slice)
        for i,item in enumerate(self[1:],start=1):
            state=state_func(item)
            last_state=slices[-1].state
            if state==last_state:
                slices[-1].increment_rightbound()
            else:
                slices.append(StateInterval({'state':state_func(item),'bounds':[i,i]}))
        return slices