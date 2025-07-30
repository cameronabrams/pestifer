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

    _required_fields = {'sourcefile', 'sourceseg', 'resseqnum1', 'resseqnum2', 'chainID'}
    _optional_fields = {'insertion1', 'insertion2', 'obj_id'}

    sourcefile: str = Field(..., description="Path to the source coordinate file containing the residues to be fused")
    sourceseg: str = Field(..., description="Segment in the source file from which residues are taken")
    resseqnum1: int = Field(..., description="N-terminal residue number of the fusion sequence")
    insertion1: str = Field('', description="Insertion code of the N-terminal residue")
    resseqnum2: int = Field(..., description="C-terminal residue number of the fusion sequence")
    insertion2: str = Field('', description="Insertion code of the C-terminal residue")
    chainID: str = Field(..., description="Chain ID of the segment in the base molecule to which the fusion is applied")
    obj_id: int = Field(0, description="Unique identifier for the Cfusion object")

    _yaml_header: ClassVar[str] = 'Cfusions'
    _objcat: ClassVar[str] = 'seq'
    _counter: ClassVar[int] = 0  # Class variable to keep track of Cfusion instances

    @staticmethod
    def _adapt(*args) -> dict:
        """
        Adapts the input to a dictionary format suitable for Cfusion instantiation.
        This method is used to convert various input types into a dictionary of parameters.
        """
        if isinstance(args[0], str):
            input_dict=Cfusion._from_shortcode(args[0])
            input_dict['obj_id'] = Cfusion._counter
            Cfusion._counter += 1
            return input_dict
        raise TypeError(f"Cannot convert {type(args[0])} to Cfusion")

    @staticmethod
    def _from_shortcode(raw: str) -> dict:
        """
        Converts a shortcode string to a dictionary of parameters for Cfusion.
        The shortcode format is: filename:C:nnn-ccc,S
        where:
        - filename is the source coordinate file
        - C is the source segment
        - nnn is the N-terminal residue number of the fusion sequence
        - ccc is the C-terminal residue number of the fusion sequence
        - S is the chain ID of the segment in the base molecule to which the fusion is applied
        """
        sourcefile, sourceseg, seq_range_chainID = raw.split(":")
        seq_range, chainID = seq_range_chainID.split(",")
        resrange = seq_range.split("-")
        resseqnum1, insertion1 = split_ri(resrange[0])
        resseqnum2, insertion2 = split_ri(resrange[1])
        return dict(
            sourcefile=sourcefile,
            sourceseg=sourceseg,
            resseqnum1=resseqnum1,
            resseqnum2=resseqnum2,
            insertion1=insertion1,
            insertion2=insertion2,
            chainID=chainID)

    def shortcode(self) -> str:
        return f"{self.sourcefile}:{self.sourceseg}:{join_ri(self.resseqnum1, self.insertion1)}-{join_ri(self.resseqnum2, self.insertion2)},{self.chainID}"

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
        return f"<CfusionList with {len(self)} items>"
