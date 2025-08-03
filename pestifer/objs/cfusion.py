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

from ..core.baseobj import BaseObj, BaseObjList
from ..core.scripters import PsfgenScripter
from .resid import ResID

class Cfusion(BaseObj):
    """
    A class for handling fusions of residues represented by an existing 
    coordinate file to the C-termini of base-molecule segments
    """

    _required_fields = {'sourcefile', 'sourceseg', 'resid1', 'resid2', 'chainID'}
    _optional_fields = {'insertion1', 'insertion2', 'obj_id'}

    sourcefile: str = Field(..., description="Path to the source coordinate file containing the residues to be fused")
    sourceseg: str = Field(..., description="Segment in the source file from which residues are taken")
    resid1: ResID = Field(..., description="Residue information for the N-terminal residue")
    resid2: ResID = Field(..., description="Residue information for the C-terminal residue")
    chainID: str = Field(..., description="Chain ID of the segment in the base molecule to which the fusion is applied")
    obj_id: int = Field(0, description="Unique identifier for the Cfusion object")

    _yaml_header: ClassVar[str] = 'Cfusions'
    _objcat: ClassVar[str] = 'seq'
    _counter: ClassVar[int] = 0  # Class variable to keep track of Cfusion instances

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Override the _adapt classmethod to handle initialization from a shortcode.
        """
        if args and isinstance(args[0], str):
            input_dict = Cfusion._from_shortcode(args[0])
            input_dict['obj_id'] = Cfusion._counter
            Cfusion._counter += 1
            return input_dict
        return super()._adapt(*args, **kwargs)

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
        resid1 = ResID(resrange[0])
        resid2 = ResID(resrange[1])
        base_dict = {
            'sourcefile': sourcefile,
            'sourceseg': sourceseg,
            'resid1': resid1,
            'resid2': resid2,
            'chainID': chainID
        }
        return base_dict

    def shortcode(self) -> str:
        return f"{self.sourcefile}:{self.sourceseg}:{self.resid1.resid}-{self.resid2.resid},{self.chainID}"

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
        W.addline(f'set fusres [atomselect $cfusid "protein and chain {self.sourceseg} and resid {self.resid1.resid} to {self.resid2.resid}"]')
        self.segfile=f'Cfusion{self.id}_{self.sourceseg}_{self.resid1.resid}_to_{self.resid2.resid}.pdb'
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
