#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the Segment class for generating and managing segments declared for psfgen
"""
from typing import ClassVar, List, TYPE_CHECKING
from pydantic import Field
import logging
logger=logging.getLogger(__name__)

from .chainidmanager import ChainIDManager
from .residue import Residue, ResidueList
from .stateinterval import StateIntervalList

from ..core.baseobj import BaseObj, BaseObjList
from ..core.objmanager import ObjManager

from ..objs.cleavagesite import CleavageSite
from ..objs.deletion import DeletionList
from ..objs.graft import GraftList
from ..objs.link import Link, LinkList
from ..objs.mutation import Mutation, MutationList
from ..objs.patch import PatchList
from ..objs.ssbond import SSBond, SSBondList

if TYPE_CHECKING:
    from ..molecule.molecule import Molecule

class Segment(BaseObj):
    """
    A class for handling segments in a molecular structure.
    This class represents a segment defined by a list of residue indices and provides methods to check if a bond intersects the segment.
    It also provides methods to yield treadmilled versions of the segment and to check for equality with another segment.
    """
    _required_fields = {'segtype', 'segname', 'chainID', 
                        'residues', 'subsegments', 'parent_chain', 
                        'specs'}
    _optional_fields = {'mutations', 'deletions', 'grafts', 'patches', 
                        'attachments', 'psfgen_segname', 'objmanager', 'parent_molecule'}
    
    segtype: str = Field(..., description="The type of the segment (e.g., 'protein', 'nucleicacid', 'glycan').")
    segname: str = Field(..., description="The name of the segment as it is understood by psfgen and in a PSF file.")
    chainID: str = Field(..., description="The chain ID associated with the segment.")
    residues: ResidueList = Field(default_factory=ResidueList, description="A list of residues that make up the segment.")
    subsegments: StateIntervalList = Field(..., description="A list of subsegments within the segment, each with a state indicating whether it is resolved or missing.")
    parent_chain: str = Field(..., description="The chain ID of the parent chain to which this segment belongs.")
    specs: dict = Field(..., description="Specifications for the segment, including details like loops and terminal modifications.")

    mutations: MutationList | None = Field(default=None, description="A list of mutations applied to the segment.")
    deletions: DeletionList | None = Field(default=None, description="A list of deletions applied to the segment.")
    grafts: GraftList | None = Field(default=None, description="A list of grafts applied to the segment.")
    patches: PatchList | None = Field(default=None, description="A list of patches applied to the segment.")
    attachments: List | None = Field(default=None, description="A list of attachments applied to the segment (Attachments are currently unimplemented).")
    psfgen_segname: str | None = Field(default=None, description="The segment name as understood by psfgen.")
    objmanager: ObjManager = Field(default_factory=ObjManager, description="The object manager for the segment, used to track and manage objects associated with the segment.")
    parent_molecule: 'Molecule' = Field(default=None, description="The parent molecule to which this segment belongs. This is set after the segment is constructed.")

    _inheritable_objs: ClassVar[set] = {'mutations', 'patches', 'Cfusions', 'grafts'}
    """
    List of inheritable modifications for the segment.
    These modifications include mutations, patches, C-terminal fusions, and grafts.
    """

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        if args: 
            segname_override = kwargs.pop('segname', 'UNSET')
            specs = kwargs.pop('specs', {})
            if isinstance(args[0], Residue):
                residue = args[0]
                apparent_chainID = residue.chainID
                apparent_segtype = residue.segtype
                apparent_segname = residue.chainID if segname_override == 'UNSET' else segname_override
                myRes = ResidueList([])
                myRes.append(residue)
                subsegments = myRes.state_bounds(lambda x: 'RESOLVED' if len(x.atoms) > 0 else 'MISSING')
                input_dict = {
                    'specs': specs,
                    'segtype': apparent_segtype,
                    'segname': apparent_segname,
                    'chainID': apparent_chainID,
                    'residues': myRes,
                    'subsegments': subsegments,
                    'parent_chain': apparent_chainID
                }
                return input_dict
            elif isinstance(args[0], ResidueList):
                residues: ResidueList = args[0]
                apparent_chainID = residues[0].chainID
                apparent_segtype = residues[0].segtype
                apparent_segname = residues[0].chainID if segname_override == 'UNSET' else segname_override
                if apparent_segtype in ['protein', 'nucleicacid']:
                    # a protein segment must have unique residue numbers
                    assert residues.puniq(['resid']), f'ChainID {apparent_chainID} has duplicate resid!'
                    # a protein segment may not have more than one protein chain
                    assert all([x.chainID==residues[0].chainID for x in residues.data])
                    residues.sort()
                    # for r in residues:
                    #     logger.debug(f'Residue {r.resname}{r.resid.resid} with {len(r.atoms)} atoms (resolved: {r.resolved})')
                    subsegments = residues.state_bounds(lambda x: 'RESOLVED' if x.resolved else 'MISSING')
                    # logger.debug(f'Segment {apparent_segname} has {len(residues)} residues across {len(subsegments)} subsegments')
                else:
                    logger.debug(f'Calling puniqify on residues of non-protein segment {apparent_segname}')
                    residues.puniquify(fields=['resid'],make_common=['chainID'])
                    count = sum([1 for x in residues if len(x.atoms)>0 and 'resid' in x.atoms[0].ORIGINAL_ATTRIBUTES])
                    if count > 0:
                        logger.debug(f'{count} residue(s) were affected by puniquify:')
                        for x in residues:
                            if len(x.atoms) > 0 and len(x.atoms[0].ORIGINAL_ATTRIBUTES) > 0:
                                logger.debug(f'    {x.chainID} {x.resname} {x.resid.resid} was resid {x.atoms[0].ORIGINAL_ATTRIBUTES["resid"].resid}')
                    # this assumes residues are in a linear sequence?  not really..
                    subsegments = residues.state_bounds(lambda x: 'RESOLVED' if len(x.atoms)>0 else 'MISSING')
                logger.debug(f'Segment {apparent_segname} has {len(residues)} residues across {len(subsegments)} subsegments')
                input_dict = {
                    'specs': specs,
                    'segtype': apparent_segtype,
                    'segname': apparent_segname,
                    'chainID': apparent_chainID,
                    'residues': residues,
                    'subsegments': subsegments,
                    'parent_chain': apparent_chainID
                }
                return input_dict
        else:
            # generate an input_dict for an empty, unconstructed segment
            # built from kwargs
            input_dict = {
                'specs': kwargs.get('specs', {}),
                'segtype': kwargs.get('segtype', 'unconstructed'),
                'segname': kwargs.get('segname', 'unconstructed'),
                'chainID': kwargs.get('chainID', 'unconstructed'),
                'residues': kwargs.get('residues', ResidueList([])),
                'subsegments': kwargs.get('subsegments', StateIntervalList([])),
                'parent_chain': kwargs.get('parent_chain', 'unconstructed')
            }
            return input_dict
        return {}

    def __str__(self):
        return f'{self.segname}: type {self.segtype} chain {self.chainID} with {len(self.residues)} residues'

    def set_parent_molecule(self, parent_molecule):
        """
        Set the parent molecule for this segment.
        
        Parameters
        ----------
        parent_molecule : object
            The parent molecule to which this segment belongs.
        """
        self.parent_molecule = parent_molecule

    def cleave(self, clv:CleavageSite, daughter_chainID):
        """
        Cleave the segment at a specified cut location, creating a daughter segment.
        
        Parameters
        ----------
        clv : Cleavage
            The cut location specifying where to cleave the segment.
        daughter_chainID : str
            The chain ID to assign to the daughter segment.
        """
        assert clv.chainID == self.segname
        assert self.segtype == 'protein'
        # r2i = self.residues.iget(resid=clv.resid2, chainID=clv.chainID)
        r2i = self.residues.iget(lambda x: x.resid == clv.resid2 and x.chainID == clv.chainID)
        # These two slice operations create *copies* of lists of residues
        # The original list of residues contains elements that are referenced
        # by other structures, like links.
        assert isinstance(self.residues, ResidueList)
        daughter_residues = self.residues[r2i:]
        assert isinstance(daughter_residues, ResidueList)
        assert daughter_residues[0].resid == clv.resid2
        parent_residues = self.residues[:r2i]
        assert daughter_residues[0] is self.residues[r2i]
        self.residues = parent_residues
        # this set_chainID call must act DEEPLY on structures in the list items
        # that have chainID attributes AND we have to revisit all links!
        daughter_residues.set(chainID=daughter_chainID)
        Dseg = Segment(daughter_residues, segname=daughter_chainID)
        assert isinstance(Dseg, Segment)
        Dseg.objmanager = self.objmanager.expel(daughter_residues)
        Dseg.set_parent_molecule(self.parent_molecule)
        # Dseg.ancestor_obj=self.ancestor_obj
        return Dseg

class SegmentList(BaseObjList[Segment]):
    """
    A list of segments in a molecular structure.
    This class represents a collection of segments and provides methods to manage and manipulate them.
    It inherits from `AncestorAwareObjList` to maintain the context of the parent molecule.
    """

    _segnames: list[str] | None = None
    _counters_by_segtype: dict[str, int] | None = None
    _daughters: dict[str, list[str]] | None = None
    _segtype_of_segname: dict[str, str] | None = None
    _segtypes_ordered: list[str] | None = None

    def describe(self) -> str:
        return f'<SegmentList with {len(self)} segments>'

    @classmethod
    def generate(cls, seq_spec: dict = {}, residues: ResidueList = [], chainIDmanager: ChainIDManager = None):
        """
        Generate a SegmentList from a sequence specification and a list of residues.
        
        Parameters
        ----------
        seq_spec : dict
            A dictionary containing sequence specifications.
        residues : ResidueList
            A list of residues to be processed into segments.
        chainIDmanager : ChainIDManager
            An object managing chain IDs to ensure uniqueness.
        
        Returns
        -------
        SegmentList
            A new SegmentList instance containing the generated segments.
        """
        self = cls([])
        self._counters_by_segtype = {}
        self._segnames = []
        self._daughters = {}
        self._segtype_of_segname = {}
        assert all([x.segtype != 'UNSET' for x in residues]), f'There are residues with UNSET segtype: {[(x.resname, x.chainID, x.resid.resid) for x in residues if x.segtype == "UNSET"]}'
        self._segtypes_ordered = []
        for r in residues:
            if not r.chainID in self._segtype_of_segname:
                self._segtype_of_segname[r.chainID] = r.segtype
            if not r.segtype in self._segtypes_ordered:
                self._segtypes_ordered.append(r.segtype)
        logger.debug(f'Generating segments from list of {len(residues)} residues. segtypes detected: {self._segtypes_ordered}')
        initial_chainIDs = list(self._segtype_of_segname.keys())
        logger.debug(f'ChainIDs detected: {initial_chainIDs}')
        chainIDmanager.sandbag(initial_chainIDs)
        for stype in self._segtypes_ordered:
            self._counters_by_segtype[stype] = 0
            # res = residues.filter(segtype=stype)
            res = residues.filter(lambda x: x.segtype == stype)
            orig_chainIDs = res.uniqattrs(['chainID'])['chainID']
            logger.debug(f'Processing {len(res)} residues of segtype {stype} in {len(orig_chainIDs)} unique chainIDs')
            # orig_res_groups = {chainID: res.filter(chainID=chainID) for chainID in orig_chainIDs}
            orig_res_groups = {chainID: res.filter(lambda x: x.chainID == chainID) for chainID in orig_chainIDs}
            logger.debug(f'Found {len(orig_res_groups)} original chainIDs for segtype {stype}: {list(orig_res_groups.keys())} {list([len(x) for x in orig_res_groups.values()])}')
            for chainID, c_res in orig_res_groups.items():
                assert isinstance(c_res, ResidueList), f'ChainID {chainID} residues are not a ResidueList: {type(c_res)}'
                this_chainID = chainID
                # c_res=res.filter(chainID=this_chainID)
                logger.debug(f'-> original chainID {chainID} in {len(c_res)} residues')
                this_chainID = chainIDmanager.check(this_chainID)
                if this_chainID != chainID:
                    logger.debug(f'{len(c_res)} residues with original chainID {chainID} and segtype {stype} are assigned new daughter chainID {this_chainID}')
                    if not chainID in self._daughters:
                        self._daughters[chainID] = []
                    self._daughters[chainID].append(this_chainID)
                    c_res.set_chainIDs(this_chainID)
                    self._segtype_of_segname[this_chainID] = stype
                    if stype == 'protein':
                        if chainID in seq_spec['build_zero_occupancy_C_termini']:
                            logger.debug(f'-> appending new chainID {this_chainID} to build_zero_occupancy_C_termini due to member {chainID}')
                            seq_spec['build_zero_occupancy_C_termini'].insert(seq_spec['build_zero_occupancy_C_termini'].index(chainID), this_chainID)
                        if chainID in seq_spec['build_zero_occupancy_N_termini']:
                            logger.debug(f'-> appending new chainID {this_chainID} to build_zero_occupancy_N_termini due to member {chainID}')
                            seq_spec['build_zero_occupancy_N_termini'].insert(seq_spec['build_zero_occupancy_N_termini'].index(chainID), this_chainID)
                        # seq_spec['build_zero_occupancy_C_termini'].remove(chainID)
                num_mis = sum([1 for x in c_res if len(x.atoms) == 0])
                thisSeg = Segment(c_res, segname=this_chainID, specs=seq_spec)
                logger.debug(f'Made segment: stype {stype} chainID {this_chainID} segname {thisSeg.segname} ({num_mis} missing) (seq_spec {seq_spec})')
                self.append(thisSeg)
                self._segnames.append(thisSeg.segname)
                self._counters_by_segtype[stype] += 1
        return self
    
    def collect_residues(self):
        """
        Collect all residues from the segments in the list.

        Returns
        -------
        ResidueList
            A list of all residues associated with the segments."""
        residues = ResidueList([])
        for seg in self.data:
            residues.extend(seg.residues)
        return residues

    def collect_working_files(self):
        """
        Collect all working files (PDB files) from the segments in the list.

        Returns
        -------
        list
            A list of PDB files associated with the segments."""
        working_files = []
        for seg in self.data:
            for subseg in seg.subsegments:
                if hasattr(subseg, 'pdb'):
                    logger.debug(f'wf: {subseg.pdb}')
                    working_files.append(subseg.pdb)
        return working_files

    def get_segment_of_residue(self, residue: Residue) -> Segment | None:
        """
        Get the segment that contains a specified residue.

        Parameters
        ----------
        residue : Residue
            The residue for which to find the containing segment.

        Returns
        -------
        Segment or None
            The segment containing the specified residue, or None if not found.
        """
        logger.debug(f'Looking for {residue.chainID}_{residue.resid.resid} in all segments {", ".join([S.segname for S in self.data])}')
        for S in self.data:
            logger.debug(f'   -> {residue.chainID}_{residue.resid.resid} in segment {S.segname}')
            # for r in S.residues:
            #     logger.debug(f'     {r.chainID}_{r.resid.resid}')
            if residue in S.residues:
                logger.debug(f' found it!')
                return S
        raise ValueError(f'Residue {residue.chainID}_{residue.resid.resid} not found in any segment!')

    def remove(self, item):
        """
        Remove a segment from the list and update associated attributes.
        
        Parameters
        ----------
        item : Segment
            The segment to remove from the list.
        """
        self._segnames.remove(item.segname)
        self._counters_by_segtype[self._segtype_of_segname[item.segname]] -= 1
        return super().remove(item)

    def prune_topology(self, mutations: MutationList, links: LinkList, ssbonds: SSBondList):
        """ 
        Prunes links, ssbonds, and even whole segments based on mutations 
        """
        pruned_objects = {'residues': ResidueList([]), 'links': LinkList([]), 'ssbonds': SSBondList([]), 'segments': SegmentList([])}
        for mutation in mutations.data:
            bad_ssbonds: SSBondList = ssbonds.filter(lambda x: (x.chainID1 == mutation.chainID and x.resid1 == mutation.resid) or (x.chainID2 == mutation.chainID and x.resid2 == mutation.resid))
            logger.debug(f'Pruning ssbonds for mutation {mutation}: {len(bad_ssbonds)} ssbonds found: {[str(s) for s in bad_ssbonds]}')
            if bad_ssbonds:
                for bond in bad_ssbonds:
                    ssbonds.remove(bond)
                    pruned_objects['ssbonds'].append(bond)
            bad_links: LinkList = links.filter(lambda x: (x.chainID1 == mutation.chainID and x.resid1 == mutation.resid) or (x.chainID2 == mutation.chainID and x.resid2 == mutation.resid))
            logger.debug(f'Pruning links for mutation {mutation}: {len(bad_links)} links found: {[str(l) for l in bad_links]}')
            if bad_links:
                Allres, Alllin = ResidueList([]), LinkList([])  
                for link in bad_links.data:
                    links.remove(link)
                    pruned_objects['links'].append(link)
                    res, lin = link.residue2.get_down_group()
                    Allres.extend(res)
                    Allres.insert(0, link.residue2)
                    Alllin.extend(lin)
                logger.debug(f'Found downstream residues for mutation {mutation}: {[f"{r.chainID}:{r.resname}_{r.resid.resid}" for r in Allres]}, links: {[str(l) for l in Alllin]}')
                if Allres:
                    logger.debug(f'Deleting residues down from and including {str(Allres[0])} due to a mutation')
                    S = self.get_segment_of_residue(Allres[0])
                    logger.debug(f'Segment {S.segname} contains residues that must be deleted because they are downstream (right) of a deleted link.')
                    for r in Allres:
                        logger.debug(f'...{str(r)}')
                        S.residues.remove(r)
                        pruned_objects['residues'].append(r)
                    if len(S.residues) == 0:
                        logger.debug(f'All residues of {S.segname} are deleted; {S.segname} is deleted')
                        self.remove(S)
                        pruned_objects['segments'].append(S)
                if Alllin:
                    for l in Alllin:
                        links.remove(l)
                        pruned_objects['links'].append(l)
        return pruned_objects
            # llist, rlist = [], []

            # left = ssbonds.get(chainID1=mutation.chainID, resid1=mutation.resid) # ssbonds for which partner 1 is the mutation
            # right = ssbonds.get(chainID2=mutation.chainID, resid2=mutation.resid) # ssbonds for which partner 2 is the mutation
            # if left:
            #     ssbonds.remove(left)
            #     pruned_objects['ssbonds'].append(left)
            # if right:
            #     ssbonds.remove(right)
            #     pruned_objects['ssbonds'].append(right)
            # left = links.get(chainID1=mutation.chainID, resid1=mutation.resid) # links for which partner 1 is the mutation
            # right = links.get(chainID2=mutation.chainID, resid2=mutation.resid) # links for which partner 2 is the mutation
            # ### it is assumed that a residue can be represented as a partner in only two shared links
            # assert isinstance(left, Link) or left is None, f'Left link for mutation {mutation} is not a Link: {type(left)}'
            # assert isinstance(right, Link) or right is None, f'Right link for mutation {mutation} is not a Link: {type(right)}'
            # if left:  # this is a link in which this mutation is partner 1
            #     links.remove(left) # get rid of this link
            #     # we need to remove residue2 and everything downstream
            #     # remove downstream residues!
            #     rlist, llist = left.residue2.get_down_group()
            #     rlist.insert(0, left.residue2)
            # elif right: # this is a link in which this mutation is the right member (should be very rare)
            #     links.remove(right)
            #     rlist, llist = right.residue2.get_down_group()
            #     rlist.insert(0, right.residue2)
            # logger.debug(f'Found upstream residues for mutation {mutation}: {[repr(r) for r in llist]}')
        #     logger.debug(f'Found downstream residues for mutation {mutation}: {[f"{r.chainID}:{r.resname}_{r.resid.resid}" for r in rlist]}')
        #     if rlist and llist:
        #         logger.debug(f'Deleting residues down from and including {str(rlist[0])} due to a mutation')
        #         S = self.get_segment_of_residue(rlist[0])
        #         logger.debug(f'Segment {S.segname} contains residues that must be deleted because they are downstream (right) of a deleted link.')
        #         for r in rlist:
        #             logger.debug(f'...{str(r)}')
        #             S.residues.remove(r)
        #             pruned_objects['residues'].append(r)
        #         if len(S.residues) == 0:
        #             logger.debug(f'All residues of {S.segname} are deleted; {S.segname} is deleted')
        #             self.remove(S)
        #             pruned_objects['segments'].append(S)
        #         for l in llist:
        #             links.remove(l)
        #             pruned_objects['links'].append(l)
        # return pruned_objects

    def inherit_objs(self, objmanager: ObjManager):
        """
        Inherit objects from the object manager for each segment in the list.
        This method updates the object manager for each segment with inherited objects
        such as mutations, grafts, and patches.
        
        Parameters
        ----------
        objmanager : ObjectManager
            The object manager from which to inherit objects.
        """
        for item in self.data:
            item.objmanager = objmanager.filter_copy(lambda x: x.chainID == item.segname, objnames=item._inheritable_objs)
            counts = item.objmanager.counts_by_header()
            logger.debug(f'Inherit objs: seg {item.segname} has {counts}')
