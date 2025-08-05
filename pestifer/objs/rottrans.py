# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for handling translation and rotation of segments in a molecular structure.
This class represents a move operation that can either be a translation or a rotation.
"""
import logging
logger = logging.getLogger(__name__)
from typing import ClassVar
from pydantic import Field
from ..core.baseobj import BaseObj, BaseObjList

class RotTrans(BaseObj):
    """
    A class for handling translation and rotation of segments in a molecular structure.
    """

    _required_fields = {'movetype'}
    """
    Required attributes for a RotTrans object.
    These attributes must be provided when creating a RotTrans object.
    
    - ``movetype``: The type of move operation, either ``trans`` for translation or ``rot`` for rotation.
    """

    _optional_fields = {'x', 'y', 'z', 'axis', 'angle'}
    """
    Optional attributes for a RotTrans object.
    These attributes may be required for certain move types but are not mandatory for all.
    - ``x``, ``y``, ``z``: The translation vector components for a translation operation.
    - ``axis``: The axis of rotation for a rotation operation, specified as a string (e.g., 'x', 'y', 'z').
    - ``angle``: The angle of rotation in degrees for a rotation operation.
    """

    _attr_choices = {'movetype': {'TRANS', 'ROT'},
                     'axis': {'x', 'y', 'z'}}

    _attr_dependencies = {'movetype': {
                            'TRANS': {'x', 'y', 'z'},
                            'ROT'  : {'axis', 'angle'}}}
    """
    This dictionary defines the dependencies between attributes in RotTrans objects.
    - For a translation operation (``movetype`` is ``trans`` or ``TRANS``), the attributes ``x``, ``y``, and ``z`` are required.
    - For a rotation operation (``movetype`` is ``rot`` or ``ROT``), the attributes ``axis`` and ``angle`` are required.
    """

    movetype: str = Field(..., description="Type of move operation, either 'trans' for translation or 'rot' for rotation.")
    x: float | None = Field(None, description="Translation vector component in the x direction.")
    y: float | None = Field(None, description="Translation vector component in the y direction.")
    z: float | None = Field(None, description="Translation vector component in the z direction.")
    axis: str | None = Field(None, description="Axis of rotation for a rotation operation (e.g., 'x', 'y', 'z').")
    angle: float | None = Field(None, description="Angle of rotation in degrees for a rotation operation.")

    _yaml_header: ClassVar[str] = 'transrot'
    """
    YAML header for RotTrans objects.
    This header is used to identify RotTrans objects in YAML files.
    """

    _objcat: ClassVar[str] = 'coord'
    """
    Category of the RotTrans object.
    This categorization is used to group RotTrans objects in the object manager.
    """

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Adapts the input to a dictionary format suitable for RotTrans instantiation.
        This method is used to convert various input types into a dictionary of parameters.
        """
        if args and isinstance(args[0], str):
            raw = args[0]
            parts = raw.split(',')
            movetype = parts[0].upper()
            if movetype == 'TRANS':
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                return dict(movetype=movetype, x=x, y=y, z=z)
            elif movetype == 'ROT':
                axis = str(parts[1])
                angle = float(parts[2])
                return dict(movetype=movetype, axis=axis, angle=angle)
            else:
                raise ValueError(f'move type {movetype} not recognized')
        return super()._adapt(*args, **kwargs)
        
    def shortcode(self) -> str:
        if self.movetype == 'TRANS':
            return f"{self.movetype},{self.x},{self.y},{self.z}"
        elif self.movetype == 'ROT':
            return f"{self.movetype},{self.axis},{self.angle}"

    def __str__(self):
        """
        Convert the RotTrans object to a shortcode representation.

        Returns
        -------
        str
            A string representation of the RotTrans object in shortcode format.
        """
        return self.shortcode()

class RotTransList(BaseObjList[RotTrans]):
    """
    A class for handling lists of RotTrans objects.
    This class inherits from BaseObjList and provides methods to manage
    a list of RotTrans objects.
    """

    def describe(self):
        """
        Describe the RotTransList object.
        
        Returns
        -------
        str
            A string description of the RotTransList object, including the number of items in the list.
        """
        return f'<RotTransList: {len(self)} items>'
