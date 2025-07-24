
from ..core.baseobj_new import BaseObj, BaseObjList
from pydantic import Field
from typing import List
import logging
logger=logging.getLogger(__name__)

class StateInterval(BaseObj):
    """
    A simple class to represent a state interval with a start and end time.
    This is useful for tracking the state of an object over a period of time.
    """
    _required_fields = ['state', 'bounds', 'build']
    state: str = Field(..., description="State of the interval")
    bounds: List[int] = Field(..., description="Bounds of the interval")
    build: bool = Field(False, description="build flag, if applicable")

    def describe(self):
        return f"StateInterval(state={self.state}, bounds={self.bounds}, build={self.build})"

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
        
class StateIntervalList(BaseObjList):
    """
    A list of StateInterval objects.
    This class extends BaseObjList to manage a collection of StateInterval instances.
    """
    def describe(self) -> str:
        return f"StateIntervalList with {len(self)} intervals"
    
    def _validate_item(self, item: StateInterval) -> None:
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
        if not isinstance(item, StateInterval):
            raise TypeError(f"Item must be an instance of StateInterval, got {type(item)}")
        
    def add_interval(self, interval: StateInterval) -> None:
        """
        Inserts a StateInterval into the calling instance.

        Parameters
        ----------
        interval : StateInterval
            the StateInterval to be inserted into the calling instance

        """
        # find a member whose bounds contain the incoming interval
        for member in self:
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
            for l,r in zip(self[:-1], self[1:]):
                assert l.state!= r.state, "Adjacent intervals must have different states"
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
        return self.bounds[-1] - self.bounds[0] + 1 if self.bounds else 0

    def describe(self) -> str:
        return f"StateInterval(state={self.state}, bounds={self.bounds}, build={self.build})"