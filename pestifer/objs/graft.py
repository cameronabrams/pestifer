# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A graft refers to a set of residues that are added to the base molecule and sourced from a separate pdb. 
The graft's "target" is a residue that is congruent to the "reference" residue in the graft.  The graft is positioned by aligning all of
its atoms using a triple on the reference and the same triplet on the target as an alignment generator. The target residue's atoms
are replaced by those in the reference and the rest of the reference is incorporated.
"""
from __future__ import annotations

import logging

from pydantic import Field
logger=logging.getLogger(__name__)
import os

from pydantic import Field
from typing import ClassVar
from ..core.baseobj_new import BaseObj, BaseObjList
from ..core.scripters import PsfgenScripter
from .link import Link, LinkList
from ..molecule.residue import Residue, ResidueList
from .resid import ResID, ResIDList

class Graft(BaseObj):
    """
    A class for handling grafts.
    """
    _required_fields = {'target_chainID', 'target_root',
                        'source_pdbid', 'source_chainID', 'source_root'}
    _optional_fields = {'source_end', 'target_partner', 'source_partner', 'source_end', 'obj_id'}

    target_chainID: str = Field(..., description="Chain ID of the target segment in the base molecule")
    target_root: ResID = Field(..., description="Root residue ID of the target segment")
    source_pdbid: str = Field(..., description="Basename of the source PDB file or PDB ID from which the graft is sourced")
    source_chainID: str = Field(..., description="Chain ID in the source PDB file")
    source_root: ResID = Field(..., description="Root residue ID of the source segment")
    """
    Required attributes for a Graft object.
    These attributes must be provided when creating a Graft object.
    
    - ``target_chainID``: Chain ID of the target segment in the base molecule.
    - ``target_root``: Resid of residue in the structure that is targeted for grafting.
    - ``source_pdbid``: Basename of the source PDB file or PDB ID from which the graft is sourced.
    - ``source_chainID``: Chain ID of the graft in the source PDB file.
    - ``source_root``: Resid of the first residue of the graft.
    """

    source_end: ResID | None = Field(None, description="End residue of the source segment")
    target_partner: ResID | None = Field(None, description="Additional residue of target segment")
    source_partner: ResID | None = Field(None, description="Additional residue of the source segment")
    obj_id: int | None = Field(0, description="Unique identifier for the Graft object")
    """
    Optional attributes for a Graft object.
    These attributes are not required when creating a Graft object.

    - ``source_end``: Last residue of the graft, which may be the same as the first if the entire source is used.
    - ``target_partner``: Additional residue in the target structure used for aligning the graft (``source_partner`` is aligned to this).
    - ``source_partner``: Additional residue in the source structure used for aligning the graft.
    - ``obj_id``: Unique identifier/tag for the graft; probably an integer.
    """
    
    _yaml_header: ClassVar[str] = 'grafts'
    _objcat: ClassVar[str] = 'seq'
    _counter: ClassVar[int] = 0  # Class variable to keep track of Graft instances

    _residues: ResidueList = None
    _source_molecule: object = None
    """
    The Molecule object from which this graft was derived.
    """

    _source_seg: object = None
    """
    The Segment object in the source molecule from which this graft was derived."""
    _mover_residues: ResidueList = None
    _index_residues: ResidueList = None
    _my_links: LinkList = None

    @staticmethod
    def _adapt(*args):
        if isinstance(args[0], str):
            input_dict=Graft._from_shortcode(args[0])
            input_dict['obj_id'] = Graft._counter
            Graft._counter += 1
            return input_dict
        raise TypeError(f"Cannot convert {type(args[0])} to Graft")

    @staticmethod
    def _from_shortcode(raw: str) -> dict:
        """
        Parse a shortcode string into an Adapter instance.
        The shortcode format is ``<target:source>``, where:
        - ``target`` has the format ``C_RRR[#AAA]``
            - C is the chainID
            - RRR is the resid of the first residue in the target alignment basis
            - AAA is the resid of the _optional_ additional residue in the target alignment basis
        - ``source`` has the format ``pdbid,C_RRR[#AAA][-NNN]``
            - pdbid is the basename of the PDB file or PDB ID (with or without the .pdb extension)
            - C is the chainID in the source PDB
            - RRR is the resid of the first residue in the source alignment basis
            - AAA is the resid of the optional additional residue in the source alignment basis
            - NNN is the optional resid of the last residue in the source alignment basis; if unset, all residues in the graft are used

        The number resids that must be specified in the source is context-dependent:

        - if the target specifies one resid, and
          - the source also specifies one resid, then it is assumed the entire graft is used. target_partner, source_partner, and source_end are all None; or
          - the source specifies two resids, then the second resid is the end resid of the graft.  target_partner and source_partner are None.
        - if the target specifies two resids, and:
          - the source specifies one, this is an error; the source must also specify an addition resid congruent to the second one in the target for alignment.
          - the source specifies two, then the second resid of the source must be the congruent resid to the second one in the target for alignment, and we assume the entire graft is used; source_end is None.
          - the source specifies three, then the second resid is congruent to the second resid of the target, and the third is the terminal resid of the source. residues to the second and third ones in the target for alignment.
        """
        try:
            target, source = raw.split(':')
        except ValueError:
            raise ValueError(f'Malformed graft shortcode {raw}')
        target_chainID, target_resrange = target.split('_')
        trr = ResIDList(target_resrange)
        nresid_in_target = len(trr)
        assert nresid_in_target in [1, 2], f'Graft target must have 1 or 2 resids, not {nresid_in_target}'
        target_root = trr[0]
        if len(trr) > 1:
            target_partner = trr[1]
        try:
            source_pdbid, source_chainID_and_resrange = source.split(',')
        except:
            raise ValueError(f'Malformed graft shortcode {raw}')
        try:
            source_chainID, source_resrange = source_chainID_and_resrange.split('_')
        except ValueError:
            raise ValueError(f'Malformed graft shortcode {raw}')
        
        srr = ResIDList(source_resrange)
        source_root = srr[0]
        nresid_in_source = len(srr)
        assert nresid_in_source in [1, 2, 3], f'Graft source must have 1, 2, or 3 resids, not {nresid_in_source}'
        base_input_dict = dict(
            target_chainID=target_chainID,
            target_root=target_root,
            source_pdbid=source_pdbid,
            source_chainID=source_chainID,
            source_root=source_root
        )
        if nresid_in_target == 1 and nresid_in_source == 1:
            return base_input_dict
        elif nresid_in_target == 1 and nresid_in_source == 2:
            base_input_dict['source_end'] = srr[1]
            return base_input_dict
        elif nresid_in_target == 1 and nresid_in_source > 2:
            raise ValueError(f'With only one resid specified in the target, you may not align to more than two residues in the source')
        elif nresid_in_target == 2 and nresid_in_source == 1:
            raise ValueError(f'With only one resid specified in the source, you may not align to two residues in the target')
        elif nresid_in_target == 2 and nresid_in_source == 2:
            base_input_dict['target_partner'] = target_partner
            base_input_dict['source_partner'] = srr[1]
            return base_input_dict
        elif nresid_in_target == 2 and nresid_in_source == 3:
            base_input_dict['target_partner'] = target_partner
            base_input_dict['source_partner'] = srr[1]
            base_input_dict['source_end'] = srr[2]
            return base_input_dict
        else:
            raise ValueError(f'Unexpected number of residues in target ({nresid_in_target}) and source ({nresid_in_source})')

    def shortcode(self) -> str:
        """
        Convert the Adapter instance to a string representation for input.
        
        Returns
        -------
        str
            A string representation of the Graft object in the format:
            target:source
        """
        target = f"{self.target_chainID}_{self.target_root.resid}"
        if self.target_partner is not None:
            target += f"#{self.target_partner.resid}"
        source = f"{self.source_pdbid},{self.source_chainID}_{self.source_root.resid}"
        if self.source_partner is not None:
            source += f"#{self.source_partner.resid}"
        if self.source_end is not None:
            source += f"-{self.source_end.resid}"
        return f"{target}:{source}"

    def activate(self, mol):
        """
        Activate the graft by linking it to a source molecule.

        Parameters
        ----------
        mol: Molecule
            The source molecule from which the graft will be sourced.
        """
        self._source_molecule = mol
        self._source_seg = self._source_molecule.asymmetric_unit.segments.get(segname=self.source_chainID)
        self._mover_residues = ResidueList([])
        self._index_residues = ResidueList([])
        # split the source residues into the index, mover, and any leftovers (which we ignore)
        for residue in self._source_seg.residues:
            if residue == self.source_root.resid:
                self._index_residues.append(residue)
            elif self.source_partner is not None and residue == self.source_partner.resid:
                self._index_residues.append(residue)
        for residue in self._source_seg.residues:
            if residue in self._index_residues:
                continue
            if residue == self.source_root.resid:
                self._mover_residues.append(residue)
            elif self.source_end is not None and residue <= self.source_end.resid:
                self._mover_residues.append(residue)
            elif self.source_end is None and self.source_partner is not None and residue == self.source_partner.resid:
                self._mover_residues.append(residue)

        # only ingest links that are internal to this set of residues or that link to index residues
        self._my_links = LinkList([])
        for l in self._source_molecule.objmanager.get('topol', {}).get('links', LinkList([])):
            if l.residue1 in self._mover_residues and l.residue2 in self._mover_residues:
                self._my_links.append(l)
            elif l.residue1 in self._index_residues and l.residue2 in self._mover_residues:
                self._my_links.append(l)
            elif l.residue2 in self._index_residues and l.residue1 in self._mover_residues:
                self._my_links.append(l)
        logger.debug(f'Activated graft {self.id}')

    def set_links(self, next_resid: ResID) -> ResID:
        """
        Set the links for the graft based on the resolved receiver residues.

        Parameters
        ----------
        next_resid: ResID
            The next available residue ID to be assigned to the graft residues.

        Returns
        -------
        next_resid: ResID
            The updated next available residue ID after setting the links.
        """
        assert self.residues is not None
        logger.debug(f'Setting links for graft {self.id}')
        logger.debug(f'-> resolved receiver residues {[repr(x) for x in self.residues]}')
        self._index_residues.set(chainID=self.residues[0].chainID)
        self._mover_residues.set(chainID=self.residues[0].chainID)
        for r in self._mover_residues:
            r.set(resseqnum=next_resid)
            next_resid += 1
        for l in self._my_links:
            if l.residue1 in self._index_residues and l.residue2 in self._mover_residues:
                l.residue1=self.residues[self._index_residues.index(l.residue1)]
            elif l.residue2 in self._index_residues and l.residue1 in self._mover_residues:
                l.residue2=self.residues[self._index_residues.index(l.residue2)]
        return next_resid

    def __str__(self):
        res=f'Graft {self.id}: index residues {[str(x) for x in self._index_residues]} delivers {len(self._mover_residues)} residues to \n'
        res+=f'  target {[str(x) for x in self.residues]} along with {len(self._my_links)} internal links:\n'
        for l in self._my_links:
            res+=f'  -> link {str(l)}\n'
        return res

    def write_pre_segment(self, W:PsfgenScripter):
        """
        Write the pre-segment Tcl commands for the graft operation.
        This method generates the Tcl commands to prepare the graft operation in VMD.
        It sets up the molecule information, selects the appropriate atoms for the graft,
        and applies a transformation to align the graft with the target segment in the base molecule.
        It also writes the grafted segment to a PDB file.

        Parameters
        ----------
        W : PsfgenScripter
            The Psfgen script writer object to which the Tcl commands will be written.
        """
        W.comment(f'{str(self)}')
        W.addline(f'set topid [molinfo ${W.molid_varname} get id]')
        if os.path.exists(f'{self.source_pdbid}.pdb'):
            W.addline(f'mol new {self.source_pdbid}.pdb')
        else:
            W.addline(f'mol new {self.source_pdbid}')
        W.addline(f'set graftid [molinfo top get id]')
        W.addline(f'mol top $topid')

        alignment_target_resid_logic = f'(resid {self.target_root.resid}'
        if self.target_partner is not None:
            alignment_target_resid_logic += f' or resid {self.target_partner.resid}'
        alignment_target_resid_logic += ')'
        
        W.addline(f'set target_sel [atomselect $topid "chain {self.target_chainID} and {alignment_target_resid_logic}"]')

        alignment_source_resid_logic = f'(resid {self.source_root.resid}'
        if self.source_partner is not None:
            alignment_source_resid_logic += f' or resid {self.source_partner.resid}'
        if self.source_end is not None:
            alignment_source_resid_logic += f' or resid {self.source_end.resid}'
        alignment_source_resid_logic += ')'
        
        W.addline(f'set source_sel [atomselect $graftid "chain {self.source_chainID} and {alignment_source_resid_logic}"]')
        W.addline(f'vmdcon -info "[$source_sel num] atoms in source, [$target_sel num] atoms in target"')
        
        # Now we need to select the residues that will be moved
        # The mover residues are those that are in the source segment *excluding* those used in the alignment

        movers_source_resid_logic = f'(not (resid {self.source_root.resid}'
        if self.source_partner is not None:
            movers_source_resid_logic += f' or resid {self.source_partner.resid})'
        if self.source_end is not None:
            movers_source_resid_logic += f' and resid <= {self.source_end.resid}'
        movers_source_resid_logic += ')'

        W.addline(f'set mover_sel [atomselect $graftid "chain {self.source_chainID} and {movers_source_resid_logic}"]')
        W.addline(f'vmdcon -info "[$mover_sel num] atoms will be moved"')
        
        # Now we need to align the mover residues to the target residues
        # We will use the transidentity command to get the transformation matrix
        # and then apply it to the mover residues

        # W.addline(f'set TT [transidentity]')
        W.addline(f'set TT [measure fit $source_sel $target_sel]')
        W.addline(f'vmdcon -info "Homog. trans. matrix: $TT"')
        W.addline(f'$mover_sel move $TT')
        self.segfile=f'graft{self.obj_id}.pdb'
        new_residlist=[]
        for y in self._mover_residues:
            new_residlist.extend([f'{y.resid.resid}' for x in y.atoms])  # r
        W.addline(f'$mover_sel set resid [list {" ".join(new_residlist)}]')
        W.addline(f'$mover_sel set chain {self._mover_residues[0].chainID}')
        W.addline(f'$mover_sel writepdb {self.segfile}')
        W.addline(f'$mover_sel delete')

    def write_in_segment(self, W:PsfgenScripter):
        """
        Write the Tcl commands to add the graft into the active segment of the base molecule.
        This method generates the Tcl commands to add the graft segment to the active segment
        of the base molecule in the Psfgen script.
        
        Parameters
        ----------
        W : PsfgenScripter
            The Psfgen script writer object to which the Tcl commands will be written."""
        W.addline (f'    pdb {self.segfile}')

    def write_post_segment(self, W:PsfgenScripter):
        """
        Write the Tcl commands to finalize the graft segment in the Psfgen script.
        This method generates the Tcl commands to finalize the graft segment in the Psfgen script.

        Parameters
        ----------
        W : PsfgenScripter
            The Psfgen script writer object to which the Tcl commands will be written.
        """
        W.addline(f'coordpdb {self.segfile} {self.chainID}')

    def assign_receiver_residues(self, Residues):
        """
        Assign the residues objects for the graft based on the original chain ID and residue sequence numbers.
        This method scans the residues in the original chain and assigns those that fall within the specified range
        of residue sequence numbers and insertions to the graft's residues attribute.

        Parameters
        ----------
        Residues: ResidueList
            The list of residues from which the graft will be assigned.
        """
        assert self._residues == None
        logger.debug(f'Assigning receiver residues to graft {self.obj_id} ({self.target_root.resid})')

        target_residue = Residues.get(chainID=self.target_chainID, resseqnum=self.target_root.resid) # change when residue definition is updated for ResID
        if target_residue is not None:
            self.residues=type(Residues)([])
            self._residues.append(target_residue)
            if self.target_partner is not None:
                target_addl_residue = Residues.get(chainID=self.target_chainID, resseqnum=self.target_partner.resid)
                if target_addl_residue is not None:
                    self._residues.append(target_addl_residue)

class GraftList(BaseObjList[Graft]):
    """
    A class for handling lists of grafts.
    This class inherits from BaseObjList and provides methods to manage a list of Graft objects.
    It allows for the assignment of residue objects to each graft and the removal of grafts that do not have any assigned residues.
    It also handles the destruction of down-links from terminal residues in the grafts.
    """
    def describe(self):
        return f"<GraftList with {len(self)} items>"

    def assign_residues(self, Residues, Links):
        """
        Assign residue objects to each graft in the list.
        This method iterates through each graft in the list, assigns residues from the provided Residues
        object, and handles the removal of down-links from terminal residues.

        Parameters
        ----------
        Residues: ResidueList
            The list of residues from which the grafts will be assigned.
        Links: LinkList
            The list of links that may need to be modified based on the assigned residues.
        """
        logger.debug(f'Assigning residue objects to {len(self)} grafts')
        delete_us=[]
        down_group=[]
        down_links=[]
        for graft in self:
            logger.debug(f' -> graft {s.id}')
            graft.assign_receiver_residues(Residues)
            # a graft that received no residue assignments means that the builder was instructed to delete
            # residues that comprise it prior to the graft operation.  We can remove it from the list.
            if graft._residues is None or len(graft._residues) == 0:
                logger.debug(f'-> removing graft {graft.obj_id} because it has no residues')
                delete_us.append(graft)
                continue
            # the graft specifications allow us to keep only the first residue in the target segment, AND if specified, the additional residue.  If there are downlinks in the first residue that do not link to the additional residue, we remove them.  If there are downlinks in the additional residue, we remove them.
            root_res = graft._residues[0]
            downs = root_res.down
            adowns = []
            addl_res = None
            if graft.target_partner is not None:
                addl_res = graft._residues[1]
                assert addl_res in downs, f'Graft {graft.obj_id} target_partner {graft.target_partner.resid} not found in down-links of root residue {root_res.chainID}_{root_res.resname}{root_res.resseqnum}{root_res.insertion}'
                adowns = addl_res.down
                # remove residues down-linked to the root residue that are not the additional residue
            for down_res in downs+adowns:
                logger.debug(f'-> removing down-link from {root_res.chainID}_{root_res.resname}{root_res.resseqnum}{root_res.insertion} to {down_res.chainID}_{down_res.resname}{down_res.resseqnum}{down_res.insertion}')
                down_group.append(Residues.remove(down_res))
                down_links.extend(Links.remove_links_to(down_res))
            # now we can remove the down-links from the root residue
        for s in delete_us:
            self.remove(s)
        return delete_us, down_group, down_links

