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
logger = logging.getLogger(__name__)

from pydantic import Field
from typing import ClassVar, TYPE_CHECKING
from ..core.baseobj import BaseObj, BaseObjList
from .link import LinkList
from ..molecule.residue import ResidueList
from .resid import ResID, ResIDList

if TYPE_CHECKING:
    from ..molecule.molecule import Molecule
    from ..molecule.segment import Segment

class Graft(BaseObj):
    """
    A class for handling grafts.
    """
    _required_fields = {'chainID', 'target_root',
                        'source_pdbid', 'source_chainID', 'source_root'}
    _optional_fields = {'source_end', 'target_partner', 'source_partner', 'source_end', 'obj_id',
                        'residues', 'source_molecule', 'source_seg', 'donor_residues', 
                        'index_residues', 'graft_residues', 'donor_internal_links', 'donor_external_links', 'segfile'}

    chainID: str = Field(..., description="Chain ID of the target segment in the base molecule")
    target_root: ResID = Field(..., description="Root residue ID of the target segment")
    source_pdbid: str = Field(..., description="Basename of the source PDB file or PDB ID from which the graft is sourced")
    source_chainID: str = Field(..., description="Chain ID in the source PDB file")
    source_root: ResID = Field(..., description="Root residue ID of the source segment")
    """
    Required attributes for a Graft object.
    These attributes must be provided when creating a Graft object.
    
    - ``chainID``: Chain ID of the target segment in the base molecule.
    - ``target_root``: Resid of residue in the structure that is targeted for grafting.
    - ``source_pdbid``: Basename of the source PDB file or PDB ID from which the graft is sourced.
    - ``source_chainID``: Chain ID of the graft in the source PDB file.
    - ``source_root``: Resid of the first residue of the graft.
    """

    source_end: ResID | None = Field(None, description="End residue of the source segment")
    target_partner: ResID | None = Field(None, description="Additional residue of target segment")
    source_partner: ResID | None = Field(None, description="Additional residue of the source segment")
    obj_id: int | None = Field(0, description="Unique identifier for the Graft object")
    residues: ResidueList | None = Field(None, description="List of residues in the graft")
    source_molecule: 'Molecule' = Field(None, description="Molecule object from which the graft was derived")
    source_seg: 'Segment' = Field(None, description="Segment object from which the graft was derived")
    donor_residues: ResidueList | None = Field(None, description="List of residues to be moved during the grafting process")
    index_residues: ResidueList | None = Field(None, description="List of residues used for indexing the graft")
    graft_residues: ResidueList | None = Field(None, description="List of residues in the graft")
    donor_internal_links: LinkList | None = Field(None, description="List of internal links associated with the graft")
    donor_external_links: LinkList | None = Field(None, description="List of external links associated with the graft")
    segfile: str | None = Field(None, description="Path to the segment file for the graft")
    """
    Optional attributes for a Graft object.

    - ``source_end``: Last residue of the graft, which may be the same as the first if the entire source is used.
    - ``target_partner``: Additional residue in the target structure used for aligning the graft (``source_partner`` is aligned to this).
    - ``source_partner``: Additional residue in the source structure used for aligning the graft.
    - ``obj_id``: Unique identifier/tag for the graft; probably an integer.

    These optional attributes are set by "activation" of the graft:

    - ``residues``: List of residue objects in the graft.
    - ``source_molecule``: Molecule object from which the graft was derived.
    - ``donor_residues``: List of residues to be donated during the grafting process.
    - ``index_residues``: List of residues used for indexing the graft.
    - ``donor_internal_links``: List of internal links associated with the graft.
    - ``donor_external_links``: List of external links associated with the graft.
    - ``segfile``: Path to the segment file for the graft.
    """
    
    _yaml_header: ClassVar[str] = 'grafts'
    _objcat: ClassVar[str] = 'seq'
    _counter: ClassVar[int] = 0  # Class variable to keep track of Graft instances

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Override the _adapt classmethod to handle initialization from a shortcode.
        This method is used to convert various input types into a dictionary of parameters.
        """
        if args and isinstance(args[0], str):
            input_dict = Graft._from_shortcode(args[0])
            input_dict['obj_id'] = Graft._counter
            Graft._counter += 1
            return input_dict
        return super()._adapt(*args, **kwargs)

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
        chainID, target_resrange = target.split('_')
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
            chainID=chainID,
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
        target = f"{self.chainID}_{self.target_root.resid}"
        if self.target_partner is not None:
            target += f"#{self.target_partner.resid}"
        source = f"{self.source_pdbid},{self.source_chainID}_{self.source_root.resid}"
        if self.source_partner is not None:
            source += f"#{self.source_partner.resid}"
        if self.source_end is not None:
            source += f"-{self.source_end.resid}"
        return f"{target}:{source}"

    def activate(self, mol: 'Molecule'):
        """
        Activate the graft by linking it to a source molecule and populating
        its native residue object lists

        Parameters
        ----------
        mol: Molecule
            The source molecule from which the graft will be sourced.
        """
        self.source_molecule = mol
        self.source_seg: 'Segment' = self.source_molecule.asymmetric_unit.segments.get(lambda x: x.segname == self.source_chainID)
        self.donor_residues = ResidueList([])
        self.index_residues = ResidueList([])
        # split the source residues into the index, mover, and any leftovers (which we ignore)
        for residue in self.source_seg.residues.data:
            if residue.resid == self.source_root:
                self.index_residues.append(residue)
            elif self.source_partner is not None and residue.resid == self.source_partner:
                self.index_residues.append(residue)
            else:
                self.donor_residues.append(residue)
        self.graft_residues = self.source_seg.residues

        self.donor_internal_links = LinkList([])
        self.donor_external_links = LinkList([])
        source_links: LinkList = self.source_molecule.objmanager.get('topol', {}).get('links', LinkList([]))
        for l in source_links.data:
            if l.residue1 in self.donor_residues.data and l.residue2 in self.donor_residues.data:
                self.donor_internal_links.append(l)
            elif l.residue1 in self.index_residues.data and l.residue2 in self.donor_residues.data:
                self.donor_external_links.append(l)
            elif l.residue2 in self.index_residues.data and l.residue1 in self.donor_residues.data:
                self.donor_external_links.append(l)
        logger.debug(f'Activated graft {self.obj_id}: {len(self.donor_residues)} mover residues and {len(self.index_residues)} index residues')

    def set_internal_resids(self, next_resid: ResID) -> ResID:
        """
        Set the internal residue IDs for the graft based on the resolved receiver residues.

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
        logger.debug(f'Graft {self.obj_id}:')
        logger.debug(f'  -> donated residues adopt chainID {self.residues[0].chainID}')
        # self.index_residues.set(chainID=self.residues[0].chainID)
        self.donor_residues.set(chainID=self.residues[0].chainID)
        for r in self.donor_residues:
            r.set(resid=next_resid)
            next_resid += 1
        for l in self.donor_external_links.data:
            if l.residue1 in self.index_residues and l.residue2 in self.donor_residues:
                l.residue1 = self.residues[self.index_residues.index(l.residue1)]
                logger.debug(f'  -> external link {str(l)}')
            elif l.residue2 in self.index_residues and l.residue1 in self.donor_residues:
                l.residue2 = self.residues[self.index_residues.index(l.residue2)]
                logger.debug(f'  -> external link {str(l)}')
        return next_resid

    def __str__(self):
        res = f'Graft {self.obj_id}:\n'
        for rr, ir in zip(self.residues.data, self.index_residues.data):
            res += f'   {str(rr):>10s} <== {str(ir):<10s}\n'
        for dr in self.donor_residues.data:
            res += f'              <== {str(dr):<10s}\n'
        res += f'{len(self.donor_internal_links)} internal links:\n'
        for l in self.donor_internal_links:
            res += f'  -> link {str(l)}\n'
        res += f'{len(self.donor_external_links)} external links:\n'
        for l in self.donor_external_links:
            res += f'  -> link {str(l)}\n'
        return res

    def assign_receiver_residues(self, Residues: ResidueList):
        """
        Populate the list of residues from the base molecule ("receiver residues") onto the graft.

        Parameters
        ----------
        Residues: ResidueList
            The list of residues from which the graft will be assigned.
        """
        assert self.residues == None
        # logger.debug(f'Assigning receiver residues to graft {self.obj_id} ({self.chainID} {self.target_root.resid})')
        # logger.debug(f'from list of {len(Residues)} residues')
        # for r in Residues.data:
        #     logger.debug(f' -> {r.chainID} {r.resid}')
        target_residue = Residues.get(lambda x: x.chainID == self.chainID and x.resid == self.target_root)
        # target_residue = Residues.get(chainID=self.chainID, resid=self.target_root)
        # logger.debug(f'-> found target residue {target_residue.resid}')
        if target_residue is not None:
            self.residues = type(Residues)([])
            self.residues.append(target_residue)
            if self.target_partner is not None:
                target_addl_residue = Residues.get(lambda x: x.chainID == self.chainID and x.resid == self.target_partner)
                if target_addl_residue is not None:
                    self.residues.append(target_addl_residue)
        else:
            raise ValueError(f'Graft {self.obj_id} target residue ({self.chainID} {self.target_root}) not found in receiver residues ({len(Residues)} total)')
        if len(self.residues) != len(self.index_residues):
            raise ValueError(f'Graft {self.obj_id} ({self.chainID} {self.target_root}) has {len(self.index_residues)} index residues but only {len(self.residues)} receiver residues assigned')
        logger.debug(f'Graft {self.obj_id} assigned {len(self.residues)} receiver residues')

class GraftList(BaseObjList[Graft]):
    """
    A class for handling lists of grafts.
    This class inherits from BaseObjList and provides methods to manage a list of Graft objects.
    It allows for the assignment of residue objects to each graft and the removal of grafts that do not have any assigned residues.
    It also handles the destruction of down-links from terminal residues in the grafts.
    """
    def describe(self):
        return f"<GraftList with {len(self)} items>"

    def assign_residues(self, Residues: ResidueList, Links: LinkList):
        """
        Assign residue objects to each graft in the list.
        Parameters
        ----------
        Residues: ResidueList
            The list of residues from which the grafts will be assigned.
        Links: LinkList
            The list of links that may need to be modified based on the assigned residues.
        """
        logger.debug(f'Assigning residue objects to {len(self)} grafts')
        delete_us = GraftList([])
        # down_group = ResidueList([])
        # down_links = LinkList([])
        for graft in self.data:
            logger.debug(f'  Graft {graft.obj_id}:')
            graft.assign_receiver_residues(Residues)
            # a graft that received no residue assignments means that the builder was instructed to delete
            # residues that comprise it prior to the graft operation.  We can remove it from the list.
            if graft.residues is None or len(graft.residues) == 0:
                logger.debug(f'-> removing graft {graft.obj_id} because it has no receivers')
                delete_us.append(graft)
                continue
            logger.debug(str(graft))
            # # the graft specifications allow us to keep only the first residue in the target segment, AND if specified, the additional residue.
            # root_res = graft.residues[0]
            # downs = root_res.down
            # adowns = []
            # addl_res = None
            # if graft.target_partner is not None:
            #     addl_res = graft.residues[1]
            #     assert addl_res in downs, f'Graft {graft.obj_id} target_partner {graft.target_partner.resid} not found in down-links of root residue {root_res.chainID}_{root_res.resname}{root_res.resid.resid}'
            #     adowns = addl_res.down
            #     # remove residues down-linked to the root residue that are not the additional residue
            # for down_res in downs + adowns:
            #     logger.debug(f'-> removing residue {down_res.chainID}_{down_res.resname}{down_res.resid.resid} from base list and graft list')
            #     down_group.append(Residues.remove_and_return(down_res))
            #     graft.residues.remove(down_res)
            #     logger.debug(f'-> removing down-link from {root_res.chainID}_{root_res.resname}{root_res.resid.resid} to {down_res.chainID}_{down_res.resname}{down_res.resid.resid}')
            #     down_links.extend(Links.remove_links_to(down_res))
            # now we can remove the down-links from the root residue
        for s in delete_us:
            self.remove(s)
        return delete_us #, down_group, down_links

    # def sequential_renumber(self, r: ResID):
    #     """
    #     Renumber the residues in each graft sequentially.
    #     This method iterates through each graft in the list and renumbers the residues
    #     in each graft sequentially, starting from 1.
    #     """
    #     logger.debug(f'Renumbering {len(self)} grafts')
    #     for g in self.data:
    #         logger.debug(f' -> renumbering graft {g.obj_id}')
    #         g.residues.renumber(start=1)
    #         g.donor_residues.renumber(start=1)
    #         g.index_residues.renumber(start=1)
    #         g.my_links.renumber(start=1)
