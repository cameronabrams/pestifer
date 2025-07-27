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
    refatom: str = Field(None, description="Name of the reference atom to align the coordinate set with. If not specified, the coordinate set is centered at the origin.")

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

    def describe(self):
        """
        Describe the Orient object.
        
        Returns
        -------
        str
            A string description of the Orient object, including axis and reference atom.
        """
        return f"Orient(axis={self.axis}, refatom={self.refatom})"
    
    class Adapter:
        """
        A class to represent the shortcode format for Orient, so that we can register to BaseObj.from_input rather than defining a local from_input.

        The shortcode format is axis,refatom
        where:
        - axis is the axis along which the coordinate set will be oriented (e.g., 'x', 'y', 'z')
        - refatom is the name of the reference atom to align the coordinate set with (optional)
        """
        def __init__(self, axis: str, refatom: str = None):
            self.axis = axis
            self.refatom = refatom

        @classmethod
        def from_string(cls, raw: str):
            parts = raw.split(',')
            if len(parts) == 1:
                return cls(parts[0])
            elif len(parts) == 2:
                return cls(parts[0], parts[1])
            else:
                raise ValueError(f"Invalid format for Orient: {raw}")
   
    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_adapter(cls, adapter: Adapter):
        """
        Create an Orient object from an Adapter instance, registered by BaseObj.from_input.
        
        Parameters
        ----------
        adapter : Adapter
            The Adapter instance containing the axis and optional reference atom.

        Returns
        -------
        Orient
            A new Orient instance created from the Adapter.
        """
        return cls(axis=adapter.axis, refatom=adapter.refatom)

    @singledispatchmethod
    @classmethod
    def new(cls, raw: Any) -> "Orient":
        pass

    @new.register(str)
    @classmethod
    def _from_str(cls, raw: str) -> "Orient":
        """
        Create a new Orient instance from a shortcode string.
        
        Parameters
        ----------
        raw : str
            The shortcode string in the format axis,refatom.
        
        Returns
        -------
        Orient
            A new instance of Orient.
        """
        adapter = cls.Adapter.from_string(raw)
        return cls._from_shortcode(adapter)
    
    def write_TcL(self,W:VMDScripter):
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
        return f'OrientList with {len(self)} orientations'
    
    def _validate_item(self, item: Orient):
        if not isinstance(item, Orient):
            raise TypeError(f"Expected Orient instance, got {type(item)}")

    def write_TcL(self,W:VMDScripter):
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
