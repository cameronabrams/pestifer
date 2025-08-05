# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A C-rotation is a transformation in which atoms are rotated around a given bond by a given amount.  The "C" 
designation means that only the "downstream" atoms of the bond are moved; upstream atoms, along with the atoms
of the bond itself, are not moved.  The set of upstream atoms and the set of downstream atoms must naturally
have no topological connection *other* than the bond itself.  Typically, this can be used to execute rotation
of a backbone loop in a C-terminal loop, a side-chain angle, usually in the service
of reducing steric clashes.  The primary job of this class is to translate the C-rotation shortcodes
specified by the user into TcL commands to be incorporated in a psfgen script.
"""
import logging
logger=logging.getLogger(__name__)

from pydantic import Field
from typing import ClassVar
from ..core.baseobj import BaseObj, BaseObjList
from .resid import ResID

class Crot(BaseObj):
    """
    A class for managing so-called "C-rotations" in a molecular structure.
    """
    _required_fields = {'angle'}
    _optional_fields = {'chainID', 'resid1', 'resid2', 'resid3', 'segname',
                        'atom1', 'atom2',
                        'segname1', 'segname2', 'segnamei', 'residi', 'atomi',
                        'segnamejk', 'residj', 'atomj', 'residk', 'atomk', 'degrees'}
    """
    Optional attributes for a Crot object.
    These attributes are used to specify additional parameters for the C-rotation:
    
    - ``chainID``: The chain ID of the segment to which the rotation applies.
    - ``resid1``: The residue resid of the first residue involved in the rotation.
    - ``resid2``: The residue resid of the second residue involved in the rotation.
    - ``resid3``: The residue resid of the third residue involved in the rotation (optional, used for ALPHA rotations).
    - ``segname``: The segment name of the segment to which the rotation applies.
    - ``atom1``: The name of the first atom involved in the rotation.
    - ``atom2``: The name of the second atom involved in the rotation.
    - ``segname1``: The segment name of the first atom involved in the rotation.
    - ``segname2``: The segment name of the second atom involved in the rotation.
    - ``segnamei``: The segment name of the first atom in the ANGLEIJK rotation.
    - ``residi``: The residue resid of the first atom in the ANGLEIJK rotation.
    - ``atomi``: The name of the first atom in the ANGLEIJK rotation.
    - ``segnamejk``: The segment name of the second atom in the ANGLEIJK rotation.
    - ``residj``: The residue resid of the second atom in the ANGLEIJK rotation.
    - ``atomj``: The name of the second atom in the ANGLEIJK rotation.
    - ``residk``: The residue resid of the third atom in the ANGLEIJK rotation.
    - ``atomk``: The name of the third atom in the ANGLEIJK rotation.
    - ``degrees``: The angle in degrees by which the atoms will be rotated.
    """
    _attr_choices = {'angle': ['PHI', 'PSI', 'OMEGA', 'CHI1', 'CHI2', 'ANGLEIJK', 'ALPHA']}
    """
    Optional attribute dependencies for the Crot object.
    This dictionary defines the dependencies for optional attributes based on the value of the ``angle`` attribute.
    The keys are the possible values of the ``angle`` attribute, and the values are lists
    of required optional attributes for that angle type.
    The dependencies are as follows:

    - ``PHI``: Requires ``chainID``, ``resid1``, and ``resid2``.
    - ``PSI``: Requires ``chainID``, ``resid1``, and ``resid2``.
    - ``OMEGA``: Requires ``chainID``, ``resid1``, and ``resid2``.
    - ``CHI1``: Requires ``chainID`` and ``resid1``.
    - ``CHI2``: Requires ``chainID`` and ``resid1``.
    - ``ANGLEIJK``: Requires ``segnamei``, ``residi``, ``atomi``, ``segnamejk``, ``residj``, ``atomj``, ``residk``, and ``atomk``.
    - ``ALPHA``: Requires ``chainID``, ``resid1``, ``resid2``, and ``resid3``.
    """
    _attr_dependencies = {'angle': {
                            'PHI':      {'chainID', 'resid1', 'resid2', 'degrees'},
                            'PSI':      {'chainID', 'resid1', 'resid2', 'degrees'},
                            'OMEGA':    {'chainID', 'resid1', 'resid2', 'degrees'},
                            'CHI1':     {'chainID', 'resid1', 'degrees'},
                            'CHI2':     {'chainID', 'resid1', 'degrees'},
                            'ANGLEIJK': {'segnamei', 'residi', 'atomi',
                                         'segnamejk', 'residj', 'atomj', 'residk', 'atomk', 'degrees'},
                            'ALPHA':    {'chainID', 'resid1', 'resid2', 'resid3'}}
                        }
    
    angle: str = Field(..., description="Type of angle to be rotated (e.g., PHI, PSI, OMEGA, CHI1, CHI2, ANGLEIJK, ALPHA)")
    chainID: str | None = Field(None, description="Chain ID of the segment to which the rotation applies")
    resid1: ResID | None = Field(None, description="Residue information of the first residue involved in the rotation")
    resid2: ResID | None = Field(None, description="Residue information of the second residue involved in the rotation")
    resid3: ResID | None = Field(None, description="Residue information of the third residue involved in the rotation")
    segname: str | None = Field(None, description="Segment name of the segment to which the rotation applies")
    atom1: str | None = Field(None, description="Name of the first atom involved in the rotation")
    atom2: str | None = Field(None, description="Name of the second atom involved in the rotation")
    segname1: str | None = Field(None, description="Segment name of the first atom involved in the rotation")
    segname2: str | None = Field(None, description="Segment name of the second atom involved in the rotation")
    segnamei: str | None = Field(None, description="Segment name of the first atom involved in the rotation")
    residi: ResID | None = Field(None, description="Residue information of the first atom in the ANGLEIJK rotation")
    atomi: str | None = Field(None, description="Name of the first atom in the ANGLEIJK rotation")
    segnamejk: str | None = Field(None, description="Segment name of the second atom in the ANGLEIJK rotation")
    residj: ResID | None = Field(None, description="Residue information of the second atom in the ANGLEIJK rotation")
    atomj: str | None = Field(None, description="Name of the second atom in the ANGLEIJK rotation")
    residk: ResID | None = Field(None, description="Residue information of the third atom in the ANGLEIJK rotation")
    atomk: str | None = Field(None, description="Name of the third atom in the ANGLEIJK rotation")
    degrees: float | None = Field(None, description="Angle in degrees by which the atoms will be rotated")

    _yaml_header: ClassVar[str] = 'crotations'
    """
    YAML header for Crot objects.
    This header is used to identify Crot objects in YAML files.
    """

    _objcat: ClassVar[str] = 'coord'
    """
    Category of the Crot object.
    This categorization is used to group Crot objects in the object manager. 
    """
    
    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Override the _adapt classmethod to handle initialization from a shortcode.
        """
        if args and isinstance(args[0], str):
            input_dict = Crot._from_shortcode(args[0])
            return input_dict
        return super()._adapt(*args, **kwargs)

    @staticmethod
    def _from_shortcode(raw: str) -> dict:
        """
        Converts the shortcode format for Crot to a dict.

        The shortcode format is:
        - For PHI, PSI, OMEGA: angle,chainID,resid1,resid2,degrees
        - For CHI1, CHI2: angle,chainID,resid1,degrees
        - For ANGLEIJK: angle,segnamei,resid1,atomi,segnamejk,resid2,atomj,resid3,atomk,degrees
        - For ALPHA: angle,chainID,resid1,resid2,[resid3]
        """
        data = raw.split(',')
        angle = data[0].upper()
        match angle:
            case 'PHI' | 'PSI' | 'OMEGA':
                chainID = data[1]
                resid1 = ResID(data[2])
                resid2 = ResID(data[3])
                degrees = float(data[4])
                return dict(angle=angle, chainID=chainID, resid1=resid1, resid2=resid2, degrees=degrees)
            case 'CHI1' | 'CHI2' :
                chainID = data[1]
                resid1 = ResID(data[2])
                degrees = float(data[3])
                return dict(angle=angle, chainID=chainID, resid1=resid1, degrees=degrees)
            case 'ANGLEIJK':
                segnamei = data[1]
                residi = ResID(data[2])
                atomi = data[3]
                segnamejk = data[4]
                residj = ResID(data[5])
                atomj = data[6]
                residk = ResID(data[7])
                atomk = data[8]
                degrees = float(data[9])
                return dict(angle=angle, segnamei=segnamei, residi=residi, atomi=atomi,
                            segnamejk=segnamejk, residj=residj, atomj=atomj,
                            residk=residk, atomk=atomk, degrees=degrees)
            case 'ALPHA':
                chainID = data[1]
                resid1 = ResID(data[2])
                resid2 = ResID(data[3])
                if len(data) < 5:
                    resid3 = resid2
                else:
                    resid3 = ResID(data[4])
                return dict(angle=angle, chainID=chainID, resid1=resid1, resid2=resid2, resid3=resid3)
            case _:
                raise ValueError(f"Unknown angle type: {angle}")

    def shortcode(self) -> str:
        match self.angle:
            case 'PHI' | 'phi' | 'PSI' | 'psi' | 'OMEGA' | 'omega':
                return f'{self.angle.upper()},{self.chainID},{self.resid1.resid},{self.resid2.resid},{self.degrees:.4f}'
            case 'CHI1' | 'chi1' | 'CHI2' | 'chi2':
                return f'{self.angle.upper()},{self.chainID},{self.resid1.resid},{self.degrees:.4f}'
            case 'ANGLEIJK' | 'angleijk':
                return f'{self.angle},{self.segnamei},{self.resid1.resid},{self.atomi},{self.segnamejk},{self.residj.resid},{self.atomj},{self.residk.resid},{self.atomk},{self.degrees:.4f}'
            case 'ALPHA' | 'alpha':
                return f'{self.angle},{self.chainID},{self.resid1.resid},{self.resid2.resid},{self.resid3.resid}'
            
class CrotList(BaseObjList[Crot]):
    """
    A class for managing lists of Crot objects.
    This class inherits from BaseObjList and provides methods to handle multiple Crot objects.
    It allows for writing Tcl commands for all Crot objects in the list.
    """

    def describe(self):
        """
        Returns a string description of the CrotList, including the number of Crot objects it contains.
        
        Returns
        -------
        str
            A string describing the CrotList.
        """
        return f"CrotList with {len(self)} Crot objects"

    # def write_TcL(self,W:PsfgenScripter,chainIDmap={},**kwargs):
    #     """
    #     Write the Tcl commands for all Crot objects in the list.
        
    #     Parameters
    #     ----------
    #     W : PsfgenScripter
    #         The Psfgen script writer object to which the Tcl commands will be written.
    #     chainIDmap : dict, optional
    #         A dictionary mapping original chain IDs to new chain IDs. If not provided, the original chain ID will be used.
    #     **kwargs : dict, optional
    #         Additional keyword arguments that may be used in the future.
    #     """
    #     for c in self:
    #         c.write_TcL(W,chainIDmap=chainIDmap,**kwargs)
