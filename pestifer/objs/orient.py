# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
This class allows the user to orient a coordinate set along a specified axis
and optionally align it with a reference atom.  The orientation is performed
using the VMD ``orient`` command, which computes the principal axes of the
coordinate set and aligns it with the specified axis.    
"""
import logging
logger = logging.getLogger(__name__)
from pydantic import Field
from typing import ClassVar

from ..core.baseobj import BaseObj, BaseObjList

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

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Override the _adapt classmethod to handle initialization from a shortcode.
        """
        if args and isinstance(args[0], str):
            parts = args[0].split(',')
            if len(parts) == 1:
                return dict(axis=parts[0])
            elif len(parts) == 2:
                return dict(axis=parts[0], refatom=parts[1])
            else:
                raise ValueError(f"Invalid format for Orient: {args[0]}")
        return super()._adapt(*args, **kwargs)
            
class OrientList(BaseObjList[Orient]):
    """
    A class for handling a list of Orient objects.
    This class inherits from BaseObjList and provides methods to write Tcl commands
    for each Orient object in the list.
    """

    def describe(self):
        return f'<OrientList with {len(self)} orientations>'
