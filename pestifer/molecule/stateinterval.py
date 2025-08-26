# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A "state interval" is a representation of a run of frames within a finitely resolved 
trajectory, all of which are in the same state.  Each trajectory frame is represented
as a counting integer index, and the state interval is defined by a pair of indices
that represent the start and end of the interval.  The state is a string that describes
the state of the system over that interval, such as "RESOLVED" or "MISSING".  The
state interval can also include additional metadata such as a build flag, selection 
name, and PDB file name.
"""

import logging

from pydantic import Field
from typing import Callable, TYPE_CHECKING

from ..core.baseobj import BaseObj, BaseObjList
if TYPE_CHECKING:
    from ..molecule.residue import Residue

logger = logging.getLogger(__name__)

class StateInterval(BaseObj):
    """
    A simple class to represent a state interval with a start and end time.
    This is useful for tracking the state of an object over a period of time.
    """
    _required_fields = {'state', 'bounds'}
    _optional_fields = {'build', 'selname', 'pdb', 'sacres'}
    _attr_choices = {
        'state': {'RESOLVED', 'MISSING', 'TESTING-A', 'TESTING-B', 'TESTING-C'},
        'build': {True, False}
    }

    state: str = Field(..., description="State of the interval")
    bounds: list[int] = Field(..., description="Bounds of the interval")
    build: bool = Field(False, description="build flag, if applicable")
    selname: str = Field(None, description="Selection name, if applicable")
    pdb: str = Field(None, description="Temporary PDB file name, if applicable")
    sacres: 'Residue' = Field(None, description="Saccharide residue, if applicable")

    def declare_buildable(self) -> None:
        """
        Declare this object as buildable.
        This method sets the build attribute to True.
        """
        self.build = True

    def increment_rightbound(self, increment: int = 1) -> None:
        """
        Increment the right bound of the interval by a specified amount.
        
        Parameters
        ----------
        increment : int, optional
            The amount to increment the right bound by (default is 1).
        """
        if self.bounds and isinstance(self.bounds[-1], int):
            self.bounds[-1] += increment
        else:
            raise ValueError("Bounds must be a list of integers with at least one element.")

    def num_items(self) -> int:
        """
        Get the number of items in the interval.

        Returns
        -------
        int
            The number of items in the interval.
        """
        return self.bounds[1] - self.bounds[0] + 1

class StateIntervalList(BaseObjList[StateInterval]):
    """
    A list of StateInterval objects.
    This class extends BaseObjList to manage a collection of StateInterval instances.
    """
    def describe(self) -> str:
        return f"StateIntervalList with {len(self)} intervals"

    @classmethod
    def process_itemlist(cls, itemlist: list[object] = None, state_func: Callable = None):
        """
        Process an item to extract StateInterval instances.

        Parameters
        ----------
        item : Any, optional
            The item to process (default is None).
        state_func : function, optional
            A function to determine the state of the item (default is None).

        Returns
        -------
        StateIntervalList :
            A list of StateInterval instances.
        """
        if itemlist is None or state_func is None:
            return None

        if not callable(state_func):
            raise ValueError("state_func must be a callable function.")

        new_list = cls()
        for i, item in enumerate(itemlist):     
            state = state_func(item)
            state_change = len(new_list) == 0 or new_list[-1].state != state
            if state_change:
                bounds = [i, i]
                interval = StateInterval(state=state, bounds=bounds)
                new_list.append(interval)
            else:
                new_list[-1].bounds[1] = i
        return new_list

    def add_interval(self, interval: StateInterval) -> None:
        """
        Inserts a StateInterval into the calling instance.

        Parameters
        ----------
        interval : StateInterval
            the StateInterval to be inserted into the calling instance

        """
        # find a member whose bounds contain the incoming interval
        for member in self.data:
            if interval.bounds[0]>member.bounds[0] and interval.bounds[1]<member.bounds[1]:
                # the new interval is completely contained within an existing interval
                break
        else:
            member = None
        if member:
            if member.state != interval.state:
                # insert the new interval into the existing member
                ins_idx = self.index(member) + 1
                self.insert(ins_idx, interval)
                baby_interval=StateInterval(
                    state=member.state,
                    bounds=[interval.bounds[1]+1, member.bounds[1]],
                    build=member.build
                )
                self.insert(ins_idx + 1, baby_interval)
                member.bounds[1] = interval.bounds[0] - 1
            else:
                # do nothing -- new interval does not change the interval list
                return
        else:
            # look for adjacent members that might permit insertion
            for l, r in zip(self.data[:-1], self.data[1:]):
                assert l.state != r.state, "Adjacent intervals must have different states"
                if l.bounds[1] > interval.bounds[0] and r.bounds[0] < interval.bounds[1]:
                    if l.state == interval.state: # grow the left member by absorbing the new interval
                        l.bounds[1] = interval.bounds[1]
                        r.bounds[0] = interval.bounds[1] + 1
                    elif r.state == interval.state: # grow the right member by absorbing the new interval
                        l.bounds[1] = interval.bounds[0] - 1
                        r.bounds[0] = interval.bounds[0]
                    else:
                        l.bounds[1] = interval.bounds[0] - 1
                        r.bounds[0] = interval.bounds[1] + 1
                        self.insert(self.index(r), interval)
                    return
            raise ValueError("No suitable member found to insert the interval into.")

    def num_items(self) -> int:
        """
        Returns the number of items in the bounds list based on index calculation
        
        Returns
        -------
        int :
            The number of items in the bounds list.
        """
        for interval in self.data:
            item_count += interval.num_items()
        return item_count
    