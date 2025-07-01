# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
CFusions are user-specified modifications which fuse residues from a
named residue range of a named protein chain of a named input coordinate
file to the C-terminus of a named segment of base molecule.
"""
import logging
logger=logging.getLogger(__name__)
import os

from functools import singledispatchmethod

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..core.scripters import PsfgenScripter
from ..core.stringthings import split_ri

class Cfusion(AncestorAwareObj):
    """
    A class for handling fusions of residues represented by an existing 
    coordinate file to the C-termini of base-molecule segments
    """

    req_attr=AncestorAwareObj.req_attr+['sourcefile','sourceseg','resseqnum1','insertion1','resseqnum2','insertion2','chainID','id']
    """
    Required attributes for a Cfusion object.
    These attributes must be provided when creating a Cfusion object.

    - ``sourcefile``: The path to the source coordinate file containing the residues to be fused.
    - ``sourceseg``: The segment identifier in the source file from which residues are taken.
    - ``resseqnum1``: The N-terminal residue number of the fusion sequence.
    - ``insertion1``: The insertion code of the N-terminal residue.
    - ``resseqnum2``: The C-terminal residue number of the fusion sequence.
    - ``insertion2``: The insertion code of the C-terminal residue.
    - ``chainID``: The chain ID of the segment in the base molecule to which the fusion is applied.
    - ``id``: A unique identifier for the Cfusion object.
    """
    
    yaml_header='Cfusions'
    """
    YAML header for Cfusion objects.
    This header is used to identify Cfusion objects in YAML files.
    """

    objcat='seq'
    """
    Category of the Cfusion object.
    This categorization is used to group Cfusion objects in the object manager.
    """

    _Cfusion_counter=0
    
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # shortcode format: filename:C:nnn-ccc,S
        # filename -- fusion coordinate filename
        # C -- source segment
        # nnn -- N-terminal resid of fusion sequence
        # ccc -- C-terminal resid of fusion sequence
        # S -- chainID, segment in base-molecule the fusion is fused to
        p1=shortcode.split(':')
        sourcefile=p1[0]
        assert os.path.exists(sourcefile),f'Cfusion sourcefile {sourcefile} not found.'
        sourceseg=p1[1]
        p2=p1[2].split(',')
        seqrange=p2[0]
        chainID=p2[1]
        seq=seqrange.split('-')
        r1,i1=split_ri(seq[0])
        r2,i2=split_ri(seq[1])
        input_dict={
            'sourcefile':sourcefile,
            'sourceseg':sourceseg,
            'resseqnum1':int(r1),
            'insertion1':i1,
            'resseqnum2':int(r2),
            'insertion2':i2,
            'chainID':chainID,
            'id':Cfusion._Cfusion_counter
        }
        Cfusion._Cfusion_counter+=1
        super().__init__(input_dict)    

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

class CfusionList(AncestorAwareObjList):
    pass
