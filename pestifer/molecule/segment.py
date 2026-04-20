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
from ..objs.link import LinkList
from ..objs.mutation import MutationList
from ..objs.patch import PatchList
from ..objs.resid import ResID
from ..objs.ssbond import SSBondList

from ..psfutil.psfcontents import PSFSegmentList

from ..util.stringthings import my_logger

if TYPE_CHECKING:
    from ..molecule.molecule import Molecule

class Segment(BaseObj):
    """
    This class represents a segment defined by a list of residue indices and provides methods to check if a bond intersects the segment.
    It also provides methods to yield treadmilled versions of the segment and to check for equality with another segment.

    A segment represents a set of residues in which there are not repeated resids.
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
                apparent_chainID = residues.data[0].chainID
                apparent_segtype = residues.data[0].segtype
                apparent_segname = residues.data[0].segname if segname_override == 'UNSET' else segname_override
                if apparent_segtype in ['protein', 'nucleicacid']:
                    # a protein segment must have unique residue numbers
                    assert residues.puniq(['resid']), f'ChainID {apparent_chainID} has duplicate resid!'
                    # a protein segment may not have more than one protein chain
                    assert all([x.chainID==residues.data[0].chainID for x in residues.data])
                    residues.sort()
                    # for r in residues:
                    #     logger.debug(f'Residue {r.resname}{r.resid.resid} with {len(r.atoms)} atoms (resolved: {r.resolved})')
                    subsegments = residues.state_bounds(lambda x: 'RESOLVED' if x.resolved else 'MISSING')
                    # logger.debug(f'Segment {apparent_segname} has {len(residues)} residues across {len(subsegments)} subsegments')
                else:
                    skip_uniquify = kwargs.get('skip_uniquify', False)
                    if not skip_uniquify:
                        logger.debug(f'Calling puniquify on {len(residues)} residues of non-polymer segment {apparent_segname} by attribute \'resid\'')
                        residues.puniquify(attrs=['resid'])
                        logger.debug(f'counting affected residues...')
                        count = sum([1 for x in residues.data if len(x.atoms)>0 and 'resid' in x.atoms.data[0].ORIGINAL_ATTRIBUTES and x.resid != x.atoms.data[0].ORIGINAL_ATTRIBUTES['resid']])
                        if count > 0:
                            logger.debug(f'{count} residue(s) were affected by puniquify:')
                            for x in residues.data:
                                if len(x.atoms) > 0 and len(x.atoms.data[0].ORIGINAL_ATTRIBUTES) > 0:
                                    logger.debug(f'    {x.chainID} {x.resname} {x.resid.resid} was resid {x.atoms.data[0].ORIGINAL_ATTRIBUTES["resid"].resid}')
                        else:
                            logger.debug(f'No duplicate resids found in {len(residues)} residues of non-polymer segment {apparent_segname}')
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
        return Dseg

class SegmentList(BaseObjList[Segment]):
    """
    A list of segments in a molecular structure.
    This class represents a collection of segments and provides methods to manage and manipulate them.
    It inherits from `AncestorAwareObjList` to maintain the context of the parent molecule.
    """

    def __init__(self, data):
        super().__init__(data)
        self.segnames = []
        self.counters_by_segtype = {}
        self.daughters = {}
        self.glycan_segment_parents: dict[str, str] = {}
        self.segtype_of_segname = {}
        self.segtypes_ordered = []
        self.seq_spec = {}
        self.residues: ResidueList = None
        self.chainIDmanager: ChainIDManager = None
        self.psfcompanion: PSFSegmentList = None
        self.links = None

    def describe(self) -> str:
        return f'<SegmentList with {len(self)} segments>'

    def generate_from_residues(self, seq_spec: dict = {}, residues: ResidueList = [],
                               chainIDmanager: ChainIDManager = None, psfcompanion: PSFSegmentList = None,
                               links=None):
        """
        Generate a SegmentList from a sequence specification and a list of residues extracted from a structure file.

        Parameters
        ----------
        seq_spec : dict
            A dictionary containing sequence specifications.
        residues : ResidueList
            A list of residues to be processed into segments.
        chainIDmanager : ChainIDManager
            An object managing chain IDs to ensure uniqueness.
        psfcompanion: PSFSegmentList
            A companion PSF segment list to provide additional context.

        Returns
        -------
        SegmentList
            A new SegmentList instance containing the generated segments.
        """
        if not all([x.segtype != 'UNSET' for x in residues.data]):
            raise ValueError(f'There are residues with UNSET segtype: {[(x.resname, x.chainID, x.resid.resid) for x in residues.data if x.segtype == "UNSET"]}')
        self.residues = residues
        self.seq_spec = seq_spec
        self.chainIDmanager = chainIDmanager
        self.links = links
        if psfcompanion is None:
            return self.build_from_only_pdb_data()
        else:
            self.psfcompanion = psfcompanion
            self.build_from_psf_and_pdb_data()

    def _build_glycan_trees(self, all_glycan_residues: list, max_glycan_size: int = 30) -> list[dict]:
        """Find connected glycan trees using residue identity, then sort by protein attachment.

        Returns a list of dicts sorted by (parent_chainID, parent_resid):
            {'residues': [Residue, ...], 'parent_chainID': str, 'parent_resid': ResID}
        """
        glycan_ids: dict[int, object] = {id(r): r for r in all_glycan_residues}
        adj: dict[int, list[int]] = {id(r): [] for r in all_glycan_residues}
        protein_anchor: dict[int, tuple] = {}  # glycan root id → (protein_chainID, protein_resid)

        if self.links:
            for link in self.links:
                r1, r2 = link.residue1, link.residue2
                if r1 is None or r2 is None:
                    continue
                r1_glycan = id(r1) in glycan_ids
                r2_glycan = id(r2) in glycan_ids
                if r1_glycan and r2_glycan:
                    adj[id(r1)].append(id(r2))
                    adj[id(r2)].append(id(r1))
                elif r1_glycan and getattr(r2, 'segtype', '') == 'protein':
                    protein_anchor[id(r1)] = (r2.chainID, r2.resid)
                elif r2_glycan and getattr(r1, 'segtype', '') == 'protein':
                    protein_anchor[id(r2)] = (r1.chainID, r1.resid)

        visited: set[int] = set()
        trees: list[dict] = []
        for r in all_glycan_residues:
            if id(r) in visited:
                continue
            component: list = []
            queue = [id(r)]
            visited.add(id(r))
            while queue:
                curr_id = queue.pop(0)
                component.append(glycan_ids[curr_id])
                for nb_id in adj[curr_id]:
                    if nb_id not in visited:
                        visited.add(nb_id)
                        queue.append(nb_id)
            parent_chainID = None
            parent_resid = None
            for c in component:
                if id(c) in protein_anchor:
                    parent_chainID, parent_resid = protein_anchor[id(c)]
                    break
            if parent_chainID is None:
                parent_chainID = component[0].chainID
                parent_resid = component[0].resid
            trees.append({
                'residues': component,
                'parent_chainID': parent_chainID,
                'parent_resid': parent_resid,
            })

        trees.sort(key=lambda t: (t['parent_chainID'], t['parent_resid'].resseqnum))
        return trees

    def build_from_only_pdb_data(self):
        """ Build the segment list given residues that were populated only from PDB data """
        for r in self.residues.data:
            if not r.chainID in self.segtype_of_segname:
                self.segtype_of_segname[r.chainID] = r.segtype
            if not r.segtype in self.segtypes_ordered:
                self.segtypes_ordered.append(r.segtype)
        logger.debug(f'Generating segments from list of {len(self.residues)} residues. segtypes detected: {self.segtypes_ordered}')
        initial_chainIDs = list(self.segtype_of_segname.keys())
        logger.debug(f'ChainIDs detected: {initial_chainIDs}')
        self.chainIDmanager.sandbag(initial_chainIDs)
        for stype in self.segtypes_ordered:
            self.counters_by_segtype[stype] = 0
            res = self.residues.filter(lambda x: x.segtype == stype)
            logger.debug(f'Processing {len(res)} residues of segtype {stype} (out of {len(self.residues)})')
            if stype == 'glycan':
                max_glycan_size = self.seq_spec.get('glycans', {}).get('max_glycan_size', 30)
                trees = self._build_glycan_trees(list(res.data), max_glycan_size)
                logger.debug(f'Found {len(trees)} glycan trees')
                glycan_local_ctr: dict[str, int] = {}
                for tree in trees:
                    parent_chainID = tree['parent_chainID']
                    component = tree['residues']
                    glycan_local_ctr[parent_chainID] = glycan_local_ctr.get(parent_chainID, 0) + 1
                    k = glycan_local_ctr[parent_chainID]
                    segname = f'{parent_chainID}G{k:02d}'
                    base_resid = (k - 1) * max_glycan_size + 1
                    for i, r in enumerate(component):
                        new_resid = ResID(base_resid + i)
                        r.resid = new_resid
                        for a in r.atoms.data:
                            a.resid = new_resid
                        r.set_chainID(parent_chainID)
                    c_res = ResidueList(component)
                    self.segtype_of_segname[segname] = stype
                    self.glycan_segment_parents[segname] = parent_chainID
                    num_mis = sum(1 for x in component if len(x.atoms) == 0)
                    thisSeg = Segment(c_res, segname=segname, specs=self.seq_spec)
                    logger.debug(f'Made glycan segment: parent_chain {parent_chainID} segname {segname} k={k} resids {base_resid}-{base_resid+len(component)-1} ({num_mis} missing)')
                    self.append(thisSeg)
                    self.segnames.append(segname)
                    self.counters_by_segtype[stype] += 1
                continue
            orig_chainIDs = res.uniqattrs(['chainID'])['chainID']
            logger.debug(f'Found {len(orig_chainIDs)} unique chainIDs for segtype {stype}: {orig_chainIDs}')
            orig_res_groups = {chainID: res.filter(lambda x: x.chainID == chainID) for chainID in orig_chainIDs}
            for chainID, c_res in orig_res_groups.items():
                assert isinstance(c_res, ResidueList), f'ChainID {chainID} residues are not a ResidueList: {type(c_res)}'
                logger.debug(f'-> original chainID {chainID} in {len(c_res)} residues')
                this_chainID = chainID
                this_chainID = self.chainIDmanager.check(this_chainID)
                if this_chainID != chainID:
                    logger.debug(f'{len(c_res)} residues with original chainID {chainID} and segtype {stype} are assigned new daughter chainID {this_chainID}')
                    if not chainID in self.daughters:
                        self.daughters[chainID] = []
                    self.daughters[chainID].append(this_chainID)
                    c_res.set_chainIDs(this_chainID)
                    self.segtype_of_segname[this_chainID] = stype
                    if stype == 'protein':
                        if chainID in self.seq_spec['build_zero_occupancy_C_termini']:
                            logger.debug(f'-> appending new chainID {this_chainID} to build_zero_occupancy_C_termini due to member {chainID}')
                            self.seq_spec['build_zero_occupancy_C_termini'].insert(self.seq_spec['build_zero_occupancy_C_termini'].index(chainID), this_chainID)
                        if chainID in self.seq_spec['build_zero_occupancy_N_termini']:
                            logger.debug(f'-> appending new chainID {this_chainID} to build_zero_occupancy_N_termini due to member {chainID}')
                            self.seq_spec['build_zero_occupancy_N_termini'].insert(self.seq_spec['build_zero_occupancy_N_termini'].index(chainID), this_chainID)
                num_mis = sum([1 for x in c_res if len(x.atoms) == 0])
                thisSeg = Segment(c_res, segname=this_chainID, specs=self.seq_spec)
                logger.debug(f'Made segment: stype {stype} chainID {this_chainID} segname {thisSeg.segname} ({num_mis} missing)')
                my_logger(self.seq_spec, logger.debug)
                self.append(thisSeg)
                self.segnames.append(thisSeg.segname)
                self.counters_by_segtype[stype] += 1
        return self

    def build_from_psf_and_pdb_data(self):
        """ Build the segment list using the pre-populated PSF segment list and PDB-derived residues """
        if not self.psfcompanion.num_residues() == len(self.residues.data):
            logger.debug(f'Number of residues in PSF: {self.psfcompanion.num_residues()}')
            logger.debug(f'Number of residues in working molecule: {len(self.residues.data)}')
        if not self.psfcompanion.num_atoms() == sum(len(x.atoms) for x in self.residues.data):
            logger.debug(f'Number of atoms in PSF: {self.psfcompanion.num_atoms()}')
            logger.debug(f'Number of atoms in working molecule: {sum(len(x.atoms) for x in self.residues.data)}')
        # for each segment in the psfcompanion, build a list of actual residues extracted from self.residues for which atom serials match those in the PSFResidues

        for seg in self.psfcompanion.data:
            matching_residues = ResidueList(list(filter(lambda x: x.segname == seg.segname, self.residues.data)))
            if not matching_residues:
                raise ValueError(f'No matching residues found for segment {seg.segname}')
            logger.debug(f'For segment {seg.segname} in PSF, found {len(matching_residues)} matching residues in PDB data, from resid {matching_residues.data[0].resid.resid} to resid {matching_residues.data[-1].resid.resid}')
            chainIDs = list(set([x.chainID for x in matching_residues.data]))
            assert len(chainIDs) == 1, f'Multiple chainIDs found for segment {seg.segname}: {chainIDs}'
            my_chainID = chainIDs[0]
            """ Since we are building from a complete psf/pdb set, there should be no
            chainID collisions, but each segment should represent only one chainID. """
            self.chainIDmanager.touch(my_chainID)
            self.append(Segment(matching_residues, segname=seg.segname, specs=self.seq_spec, skip_uniquify=True))

    def get_residue_to_segname_map(self) -> dict:
        """Return a mapping from id(residue) → segname for all residues in all segments."""
        result = {}
        for S in self.data:
            for r in S.residues.data:
                result[id(r)] = S.segname
        return result

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
            if any(r is residue for r in S.residues):
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
        self.segnames.remove(item.segname)
        self.counters_by_segtype[self.segtype_of_segname[item.segname]] -= 1
        return super().remove(item)

    def prune_topology(self, mutations: MutationList, links: LinkList, ssbonds: SSBondList):
        """
        Prunes links, ssbonds, and even whole segments based on mutations
        """
        pruned_objects = {'residues': ResidueList([]), 'links': LinkList([]), 'ssbonds': SSBondList([]), 'segments': SegmentList([])}

        bad_ssbonds = SSBondList([])
        for ssbond in ssbonds.data:
            mutation1 = mutations.get(lambda x: x.chainID == ssbond.chainID1 and x.resid == ssbond.resid1)
            mutation2 = mutations.get(lambda x: x.chainID == ssbond.chainID2 and x.resid == ssbond.resid2)
            logger.debug(f'Checking disulfide bond {ssbond} for mutations...')
            logger.debug(f'   mutation1: {mutation1}, mutation2: {mutation2}')
            if mutation1 and mutation1.newresname != 'CYS':
                bad_ssbonds.append(ssbond)
            elif mutation2 and mutation2.newresname != 'CYS':
                bad_ssbonds.append(ssbond)

        for ssbond in bad_ssbonds.data:
            ssbonds.remove(ssbond)
            pruned_objects['ssbonds'].append(ssbond)

        bad_links = LinkList([])
        for link in links.data:
            mutation1 = mutations.get(lambda x: x.chainID == link.chainID1 and x.resid == link.resid1)
            mutation2 = mutations.get(lambda x: x.chainID == link.chainID2 and x.resid == link.resid2)
            if mutation1 or mutation2:
                bad_links.append(link)

        if len(bad_links) > 0:
            Alllin = LinkList([])
            for link in bad_links.data:
                links.remove(link)
                pruned_objects['links'].append(link)
                res, lin = link.residue2.get_down_group()
                to_delete = ResidueList([link.residue2] + list(res))
                Alllin.extend(lin)
                if to_delete:
                    logger.debug(f'Deleting residues down from and including {str(to_delete[0])} due to a mutation')
                    S = self.get_segment_of_residue(to_delete[0])
                    logger.debug(f'Segment {S.segname} contains residues that must be deleted because they are downstream (right) of a deleted link.')
                    for r in to_delete:
                        logger.debug(f'...{str(r)}')
                        S.residues.remove_instance(r)
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
            logger.debug(f'Segment {item.segname} inherits:')
            my_logger(counts, logger.debug, depth=1)
