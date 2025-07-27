# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for handling translation and rotation of segments in a molecular structure.
This class represents a move operation that can either be a translation or a rotation.
"""
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod
from typing import ClassVar, Any
from pydantic import Field
from ..core.baseobj_new import BaseObj, BaseObjList
from ..core.scripters import VMDScripter

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

    _attr_dependencies = {'movetype': {'TRANS': {'x', 'y', 'z'},
                            'ROT': {'axis', 'angle'}}}
    """
    This dictionary defines the dependencies between attributes in RotTrans objects.
    - For a translation operation (``movetype`` is ``trans`` or ``TRANS``), the attributes ``x``, ``y``, and ``z`` are required.
    - For a rotation operation (``movetype`` is ``rot`` or ``ROT``), the attributes ``axis`` and ``angle`` are required.
    """

    movetype: str = Field(..., description="Type of move operation, either 'trans' for translation or 'rot' for rotation.")
    x: float = Field(None, description="Translation vector component in the x direction.")
    y: float = Field(None, description="Translation vector component in the y direction.")
    z: float = Field(None, description="Translation vector component in the z direction.")
    axis: str = Field(None, description="Axis of rotation for a rotation operation (e.g., 'x', 'y', 'z').")
    angle: float = Field(None, description="Angle of rotation in degrees for a rotation operation.")

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

    def describe(self):
        """
        Describe the RotTrans object.
        
        Returns
        -------
        str
            A string description of the RotTrans object, including move type and parameters.
        """
        if self.movetype=='TRANS':
            return f"RotTrans(movetype={self.movetype}, x={self.x}, y={self.y}, z={self.z})"
        elif self.movetype=='ROT':
            return f"RotTrans(movetype={self.movetype}, axis={self.axis}, angle={self.angle})"
        else:
            return f"RotTrans(movetype={self.movetype})"
        
    class Adapter:
        """
        A class to represent the shortcode format for RotTrans, so that we can register to BaseObj.from_input rather than defining a local from_input.
        """
        def __init__(self, movetype: str, x: float = None, y: float = None, z: float = None, axis: str = None, angle: float = None):
            self.movetype = movetype.upper()
            self.x = x
            self.y = y
            self.z = z
            self.axis = axis
            self.angle = angle

        @classmethod
        def from_string(cls, raw: str):
            # Parse the raw string to create a RotTrans object.
            parts = raw.split(',')
            movetype = parts[0].upper()
            if movetype == 'TRANS':
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                return cls(movetype=movetype, x=x, y=y, z=z)
            elif movetype == 'ROT':
                axis = str(parts[1])
                angle = float(parts[2])
                return cls(movetype=movetype, axis=axis, angle=angle)
            else:
                logger.warning(f'move type {movetype} not recognized')
                return None
        
        def to_dict(self) -> dict:
            """
            Convert the Adapter instance to a dictionary representation.
            
            Returns
            -------
            dict
                A dictionary representation of the RotTrans object.
            """
            return {k: v for k, v in self.__dict__.items() if v is not None}    
        
        def to_string(self) -> str:
            """
            Convert the Adapter instance to a string representation.
            
            Returns
            -------
            str
                A string representation of the RotTrans object in shortcode format.
            """
            if self.movetype == 'TRANS':
                return f"{self.movetype},{self.x},{self.y},{self.z}"
            elif self.movetype == 'ROT':
                return f"{self.movetype},{self.axis},{self.angle}"
            else:
                logger.warning(f'move type {self.movetype} not recognized')
                return None

    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_adapter(cls, adapter: Adapter):
        """
        Create a RotTrans object from an Adapter instance, registered by BaseObj.from_input.
        
        Parameters
        ----------
        adapter : Adapter
            An instance of the Adapter class containing the necessary attributes.
        
        Returns
        -------
        RotTrans
            A new RotTrans object created from the attributes of the adapter.
        """
        input_dict = adapter.to_dict()
        return cls(**input_dict)
    
    @singledispatchmethod
    @classmethod
    def new(cls, raw: Any) -> "RotTrans":
        raise TypeError(f"Cannot create RotTrans from {type(raw)}: {raw}")

    @new.register(str)
    @classmethod
    def _from_string(cls, raw: str) -> "RotTrans":
        """
        Create a new RotTrans instance from a shortcode string.

        Parameters
        ----------
        raw : str
            The shortcode string in the format movetype,x,y,z or movetype,axis,angle.

        Returns
        -------
        RotTrans
            A new instance of RotTrans.
        """
        adapter = cls.Adapter.from_string(raw)
        instance = cls._from_adapter(adapter)
        return instance

    def to_input_string(self) -> str:
        """
        Converts the RotTrans object to a string representation for input.
        
        Returns
        -------
        str
            A string representation of the RotTrans object in shortcode format.
        """
        return self.Adapter(**(self.model_dump())).to_string()

    def __str__(self):
        """
        Convert the RotTrans object to a shortcode representation.

        Returns
        -------
        str
            A string representation of the RotTrans object in shortcode format.
        """
        return self.to_input_string()

    def write_TcL(self,W:VMDScripter,**kwargs):
        """
        Write the Tcl commands to perform the move operation in VMD.
        This method generates the Tcl commands to either translate or rotate a segment in VMD.
        
        Parameters
        ----------
        W : VMD
            The VMD script writer object to which the Tcl commands will be written.
        **kwargs : dict
            Additional keyword arguments (not used in this method).
        """
        molid_varname=W.molid_varname
        molid=f'${molid_varname}'
        W.addline(f'set mover [atomselect {molid} all]')
        if self.movetype=='TRANS':
            W.addline(f'$mover moveby [list {self.x} {self.y} {self.z}]')
        elif self.movetype=='ROT':
            W.addline(f'set COM [measure center $mover weight mass]')
            W.addline(f'$mover move [trans origin $COM axis {self.axis} {self.angle}]')

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

    def write_TcL(self,W:VMDScripter,**kwargs):
        """
        Write the Tcl commands for each RotTrans object in the list.
        This method iterates over each RotTrans object in the list and calls its `write_TcL` method
        to generate the Tcl commands. The commands are written to the provided VMD script writer object.
        
        Parameters
        ----------
        W : VMD
            The VMD script writer object to which the Tcl commands will be written.
        **kwargs : dict
            Additional keyword arguments (not used in this method).
        """
        for c in self:
            c.write_TcL(W,**kwargs)