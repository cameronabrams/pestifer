# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
CFusions are user-specified modifications which fuse residues from a
named residue range of a named protein chain of a named input coordinate
file to the C-terminus of a named segment of base molecule.
"""
import logging
logger=logging.getLogger(__name__)

from pydantic import Field
from typing import ClassVar

from ..core.baseobj_new import BaseObj, BaseObjList
from ..core.scripters import PsfgenScripter
from ..core.stringthings import split_ri, join_ri

class Cfusion(BaseObj):
    """
    A class for handling fusions of residues represented by an existing 
    coordinate file to the C-termini of base-molecule segments
    """

    _required_fields = ['sourcefile', 'sourceseg', 'resseqnum1', 'insertion1', 'resseqnum2', 'insertion2', 'chainID', 'id']
    _optional_fields = ['yaml_header', 'objcat']

    sourcefile: str = Field(..., description="Path to the source coordinate file containing the residues to be fused")
    sourceseg: str = Field(..., description="Segment in the source file from which residues are taken")
    resseqnum1: int = Field(..., description="N-terminal residue number of the fusion sequence")
    insertion1: str = Field(..., description="Insertion code of the N-terminal residue")
    resseqnum2: int = Field(..., description="C-terminal residue number of the fusion sequence")
    insertion2: str = Field(..., description="Insertion code of the C-terminal residue")
    chainID: str = Field(..., description="Chain ID of the segment in the base molecule to which the fusion is applied")
    id: int = Field(..., description="Unique identifier for the Cfusion object")

    _yaml_header: ClassVar[str] = 'Cfusions'
    _objcat: ClassVar[str] = 'seq'
    _counter: ClassVar[int] = 0  # Class variable to keep track of Cfusion instances

    def describe(self):
        return f"Cfusion(sourcefile={self.sourcefile}, sourceseg={self.sourceseg}, resseqnum1={self.resseqnum1}, insertion1={self.insertion1}, resseqnum2={self.resseqnum2}, insertion2={self.insertion2}, chainID={self.chainID}, id={self.id})"
    
    class Adapter:
        """
        A class to represent the shortcode format for Cfusion, so that we can register to BaseObj.from_input rather than defining a local from_input.

        The shortcode format is filename:C:nnn-ccc,S
        where:
        - filename is the fusion coordinate filename
        - C is the source segment
        - nnn is the N-terminal residue of the fusion sequence
        - ccc is the C-terminal residue of the fusion sequence
        - S is the chainID, segment in base-molecule the fusion is fused to
        """
        def __init__(self, sourcefile: str, sourceseg: str, resseqnum1: int, insertion1: str, resseqnum2: int, insertion2: str, chainID: str, id: int = 0):
            self.sourcefile = sourcefile
            self.sourceseg = sourceseg
            self.resseqnum1 = resseqnum1
            self.insertion1 = insertion1
            self.resseqnum2 = resseqnum2
            self.insertion2 = insertion2
            self.chainID = chainID
            self.id = id

        @classmethod
        def from_string(cls, raw: str):
            sourcefile, sourceseg, seq_range_chainID = raw.split(":")
            seq_range, chainID = seq_range_chainID.split(",")
            resrange = seq_range.split("-")
            resseqnum1, insertion1 = split_ri(resrange[0])
            resseqnum2, insertion2 = split_ri(resrange[1])
            return cls(sourcefile, sourceseg, resseqnum1, insertion1, resseqnum2, insertion2, chainID)

        def to_string(self) -> str:
            return f"{self.sourcefile}:{self.sourceseg}:{join_ri(self.resseqnum1,self.insertion1)}-{join_ri(self.resseqnum2,self.insertion2)},{self.chainID}"
        
        def to_dict(self) -> dict:
            return {k: v for k, v in self.__dict__.items() if not v is None}

    # if registering to BaseObj.from_input, the type must be unique to this class
    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_str(cls, shortcode: Adapter):
        input_dict=shortcode.to_dict()
        input_dict['id'] = cls._counter  # Use the class variable to set the id
        # Increment the counter for the next instance
        cls._counter+=1
        return cls(**input_dict)

    def to_input_string(self) -> str:
        """
        Converts the Cfusion object to a string representation for input.
        
        Returns
        -------
        str
            A string representation of the Cfusion object in the format:
            filename:C:nnn-ccc,S
        """
        # return self.Adapter(self.sourcefile, self.sourceseg, self.resseqnum1, self.insertion1, self.resseqnum2, self.insertion2, self.chainID).to_string()
        return self.Adapter(**(self.model_dump())).to_string()

    def write_pre_segment(self, W: PsfgenScripter):
        """
        Writes the Tcl commands to create a fusion segment in the Psfgen script.

        Parameters
        ----------
        W : Psfgen
            The Psfgen script writer object to which the Tcl commands will be written.
            This method generates the Tcl commands to create a new segment in the Psfgen
            script for the fusion of residues from the source coordinate file.
        """

        W.addline(f'set topid [molinfo top get id]')
        W.addline(f'mol new {self.sourcefile}')
        W.addline(f'set cfusid [molinfo top get id]')
        W.addline(f'mol top $topid')
        W.addline(f'set fusres [atomselect $cfusid "protein and chain {self.sourceseg} and resid {self.resseqnum1}{self.insertion1} to {self.resseqnum2}{self.insertion2}"]')
        self.segfile=f'Cfusion{self.id}_{self.sourceseg}_{self.resseqnum1}{self.insertion1}_to_{self.resseqnum2}{self.insertion2}.pdb'
        W.addline(f'$fusres writepdb {self.segfile}')
        W.addline(f'delete $cfusid')

    def write_in_segment(self,W:PsfgenScripter):
        """
        Writes the Tcl commands to add the fusion into the active segment of the base molecule.

        Parameters
        ----------
        W : Psfgen
            The Psfgen script writer object to which the Tcl commands will be written.
            This method generates the Tcl commands to add the fusion segment to the active segment
            of the base molecule in the Psfgen script.
        """
        W.addline (f'    pdb {self.segfile}')

    def write_post_segment(self,W:PsfgenScripter):
        """
        Writes the Tcl commands to finalize the fusion segment in the Psfgen script.

        Parameters
        ----------
        W : Psfgen
            The Psfgen script writer object to which the Tcl commands will be written.
            This method generates the Tcl commands to finalize the fusion segment in the Psfgen script.
        """
        W.addline(f'coordpdb {self.segfile} {self.chainID}')

class CfusionList(BaseObjList[Cfusion]):
    def describe(self) -> str:
        return f"CfusionList with {len(self)} cfusions"
    def _validate_item(self, item: Cfusion) -> None:
        """
        Validate that the item is an instance of Cfusion.
        """
        if not isinstance(item, Cfusion):
            raise TypeError(f"Expected Cfusion instance, got {type(item)}")
