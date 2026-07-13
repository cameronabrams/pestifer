# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for handling translation and rotation of segments in a molecular structure.
This class represents a move operation that can either be a translation or a rotation.
"""
import logging
logger = logging.getLogger(__name__)
from typing import ClassVar
from pydantic import Field, model_validator
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

    _optional_fields = {'x', 'y', 'z', 'axis', 'angle', 'sel', 'source', 'target', 'axis_atoms'}
    """
    Optional attributes for a RotTrans object.
    These attributes may be required for certain move types but are not mandatory for all.
    - ``x``, ``y``, ``z``: The translation vector components for a translation operation.
    - ``axis``: The axis of rotation for a rotation operation, specified as a string (e.g., 'x', 'y', 'z').
    - ``angle``: The angle of rotation in degrees for a rotation operation.
    - ``sel``: A VMD atomselection string naming the fragment to transform (default ``all``).  The
      fragment must be fully disconnected from the rest of the system; a disconnection guard in the
      generated script hard-errors if any bond crosses the selection boundary.
    - ``source``, ``target``: For an ``ALIGN`` operation, the source and target vectors.  Each is
      either a literal 3-vector ``[x, y, z]`` or a pair of VMD atomselections ``[selA, selB]`` whose
      mass-weighted centers define the vector (from ``selA`` to ``selB``).  The fragment is rotated
      about its own center of mass by the minimal (roll-free) rotation carrying ``source`` onto
      ``target``.
    - ``axis_atoms``: For an ``AXISANGLE`` operation, a list of three VMD atomselections
      ``[selI, selJ, selK]``.  The fragment is rotated by ``angle`` degrees about the axis normal to
      the I-J-K angle (``(rI - rJ) x (rJ - rK)``), pivoting at the center of ``selJ``.  (This is the
      rigid-body rotation formerly offered as the ``ANGLEIJK`` irotation.)
    """

    _attr_choices = {'movetype': {'TRANS', 'ROT', 'ALIGN', 'AXISANGLE'},
                     'axis': {'x', 'y', 'z'}}

    _attr_dependencies = {'movetype': {
                            'TRANS': {'x', 'y', 'z'},
                            'ROT'  : {'axis', 'angle'},
                            'ALIGN': {'source', 'target'},
                            'AXISANGLE': {'axis_atoms', 'angle'}}}
    """
    This dictionary defines the dependencies between attributes in RotTrans objects.
    - For a translation operation (``movetype`` is ``trans`` or ``TRANS``), the attributes ``x``, ``y``, and ``z`` are required.
    - For a rotation operation (``movetype`` is ``rot`` or ``ROT``), the attributes ``axis`` and ``angle`` are required.
    - For a vector-alignment operation (``movetype`` is ``ALIGN``), the attributes ``source`` and ``target`` are required.
    - For an axis-angle rotation (``movetype`` is ``AXISANGLE``), the attributes ``axis_atoms`` and ``angle`` are required.
    """

    movetype: str = Field(..., description="Type of move operation: 'TRANS' (translate), 'ROT' (rotate about a principal axis), 'ALIGN' (rotate to carry a source vector onto a target vector), or 'AXISANGLE' (rotate about the axis normal to a three-atom angle).")
    x: float | None = Field(None, description="Translation vector component in the x direction.")
    y: float | None = Field(None, description="Translation vector component in the y direction.")
    z: float | None = Field(None, description="Translation vector component in the z direction.")
    axis: str | None = Field(None, description="Axis of rotation for a rotation operation (e.g., 'x', 'y', 'z').")
    angle: float | None = Field(None, description="Angle of rotation in degrees for a rotation operation.")
    sel: str | None = Field(None, description="VMD atomselection naming the fragment to transform (default 'all').")
    source: list | None = Field(None, description="ALIGN source vector: a literal [x,y,z] or a pair of atomselections [selA, selB].")
    target: list | None = Field(None, description="ALIGN target vector: a literal [x,y,z] or a pair of atomselections [selA, selB].")
    axis_atoms: list | None = Field(None, description="AXISANGLE: three VMD atomselections [selI, selJ, selK]; the fragment rotates about the axis normal to the I-J-K angle, pivoting at selJ.")

    _yaml_header: ClassVar[str] = 'transrot'
    """
    YAML header for RotTrans objects.
    This header is used to identify RotTrans objects in YAML files.
    """

    _yaml_aliases: ClassVar[list] = ['rottrans']
    """
    Alternate YAML headers accepted for RotTrans objects.  ``rottrans`` is accepted as a synonym for
    ``transrot`` so either spelling may be used in configuration files.
    """

    _objcat: ClassVar[str] = 'coord'
    """
    Category of the RotTrans object.
    This categorization is used to group RotTrans objects in the object manager.
    """

    @staticmethod
    def _validate_vecspec(name: str, spec: list) -> None:
        """Validate an ALIGN vector spec: a literal 3-vector of numbers, or a 2-element pair of
        atomselection strings.  Raises ValueError with a clear message otherwise."""
        if not isinstance(spec, (list, tuple)):
            raise ValueError(f"ALIGN '{name}' must be a list, not {type(spec).__name__}")
        if len(spec) == 3 and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in spec):
            return
        if len(spec) == 2 and all(isinstance(v, str) for v in spec):
            return
        raise ValueError(f"ALIGN '{name}' must be a literal 3-vector [x,y,z] of numbers or a pair "
                         f"of atomselections [selA, selB]; got {spec!r}")

    @model_validator(mode='after')
    def _check_align_vectors(self):
        if self.movetype == 'ALIGN':
            self._validate_vecspec('source', self.source)
            self._validate_vecspec('target', self.target)
        if self.movetype == 'AXISANGLE':
            if (not isinstance(self.axis_atoms, (list, tuple)) or len(self.axis_atoms) != 3
                    or not all(isinstance(v, str) for v in self.axis_atoms)):
                raise ValueError(f"AXISANGLE 'axis_atoms' must be a list of three atomselection "
                                 f"strings [selI, selJ, selK]; got {self.axis_atoms!r}")
        return self

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
        elif self.movetype == 'ALIGN':
            return f"{self.movetype},source={self.source},target={self.target}"
        elif self.movetype == 'AXISANGLE':
            return f"{self.movetype},axis_atoms={self.axis_atoms},angle={self.angle}"

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
