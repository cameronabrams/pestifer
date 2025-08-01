# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
This class allows the user to orient a coordinate set along a specified axis
and optionally align it with a reference atom.  The orientation is performed
using the VMD ``orient`` command, which computes the principal axes of the
coordinate set and aligns it with the specified axis.    
"""
import logging
logger=logging.getLogger(__name__)
from pydantic import Field
from typing import ClassVar, Any
from functools import singledispatchmethod

from ..core.baseobj_new import BaseObj, BaseObjList
from ..core.scripters import VMDScripter

class Orient(BaseObj):
    """
    A class for handling orientation of a coordinate set in VMD.
    """

    _required_fields = {'axis'}
    """

    Required attributes for an Orient object.
    These attributes must be provided when creating an Orient object.

    - ``axis``: The axis along which the coordinate set will be oriented.
    """
    _optional_fields = {'refatom'}
    """

    Optional attributes for an Orient object.
    These attributes may be present but are not required.

    - ``refatom``: The name of the reference atom to align the coordinate set with. If not specified,
      the coordinate set is centered at the origin.
    """

    axis: str = Field(..., description="Axis along which the coordinate set will be oriented (e.g., 'x', 'y', 'z')")
    refatom: str | None = Field(None, description="Name of the reference atom to align the coordinate set with. If not specified, the coordinate set is centered at the origin.")

    _yaml_header: ClassVar[str] = 'orient'
    """
    YAML header for Orient objects.
    This header is used to identify Orient objects in YAML files.
    """

    _objcat: ClassVar[str] = 'coord'
    """
    Category of the Orient object.
    This categorization is used to group Orient objects in the object manager.
    """

    @staticmethod
    def _adapt(*args) -> dict:
        """
        Adapts the input to a dictionary format suitable for Orient instantiation.
        This method is used to convert various input formats into a dictionary that can be used to create an Orient object.
        """
        if isinstance(args[0], str):
            parts = args[0].split(',')
            if len(parts) == 1:
                return dict(axis=parts[0])
            elif len(parts) == 2:
                return dict(axis=parts[0], refatom=parts[1])
            else:
                raise ValueError(f"Invalid format for Orient: {args[0]}")

    def write_TcL(self, W:VMDScripter):
        """
        Write the Tcl commands to orient the coordinate set in VMD.
        This method generates the Tcl commands to orient the coordinate set along the specified axis
        and optionally align it with a reference atom. The commands are written to the provided VMD
        script writer object.
        
        Parameters
        ----------
        W: VMD
            The VMD script writer object to which the Tcl commands will be written.
        """
        W.addline('set a [atomselect top all]')
        W.addline('set I [draw principalaxes $a]')
        adict=dict(x=r'{1 0 0}',y=r'{0 1 0}',z=r'{0 0 1}')
        W.addline(f'set A [orient $a [lindex $I 2] {adict[self.axis]}]')
        W.addline(r'$a move $A')
        if hasattr(self,'refatom'):
            W.addline(r'set com [measure center $a]')
            W.addline(r'$a moveby [vecscale $com -1]')
            W.addline(f'set z [[atomselect top "name {self.refatom}"] get z]')
            W.addline(r'if { $z < 0.0 } {')
            W.addline(r'   $a move [transaxis x 180 degrees]')
            W.addline(r'}')
            
class OrientList(BaseObjList[Orient]):
    """
    A class for handling a list of Orient objects.
    This class inherits from BaseObjList and provides methods to write Tcl commands
    for each Orient object in the list.
    """

    def describe(self):
        return f'<OrientList with {len(self)} orientations>'

    def write_TcL(self, W:VMDScripter):
        """
        Write the Tcl commands for each Orient object in the list.
        This method iterates over each Orient object in the list and calls its `write_TcL` method   
        to generate the Tcl commands. The commands are written to the provided VMD script writer object.
        
        Parameters
        ----------
        W: VMD
            The VMD script writer object to which the Tcl commands will be written.
        """
        for c in self:
            c.write_TcL(W)
