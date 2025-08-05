# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for building the asymmetric unit from a PDB file
"""
import logging
logger = logging.getLogger(__name__)

from argparse import Namespace
from mmcif.api.PdbxContainers import DataContainer
from pathlib import Path
from pidibble.pdbrecord import PDBRecordDict
from typing import ClassVar

from .atom import AtomList, Atom, Hetatm
from .chainidmanager import ChainIDManager
from .residue import Residue, ResidueList, ResiduePlaceholderList
from .segment import SegmentList

from ..core.objmanager import ObjManager

from ..objs.graft import GraftList
from ..objs.link import LinkList
from ..objs.mutation import Mutation, MutationList
from ..objs.patch import PatchList
from ..objs.seqadv import SeqadvList
from ..objs.ssbond import SSBondList
from ..objs.ter import TerList

from ..psfutil.psfcontents import PSFContents

from ..util.util import write_residue_map

class AsymmetricUnit:
    """
    A class for building the asymmetric unit from a PDB/mmCIF file.

    Parameters
    ----------
    parsed : PDBRecordDict or DataContainer
        The parsed PDB or mmCIF data.  A PDBRecordDict is returned by `pidibble.pdbparse.PDBParser(spec).parse().parsed`, and a DataContainer is returned by `~pestifer.util.cifutil.CIFload(filepath)`.
    sourcespecs : dict, optional
        Specifications for the source of the data, including exclusions and sequence specifications.
    objmanager : ObjManager, optional
        An object manager to handle the objects in the asymmetric unit. Defaults to a new instance of `ObjManager`.
    chainIDmanager : ChainIDManager, optional
        A manager for chain IDs. Defaults to a new instance of `ChainIDManager`.
    psf : str or Path, optional
        Path to a PSF file to use for additional information about the structure. Defaults to None.  If set, it must be congruent to the PDB/mmCIF data.
    """

    _residue_excludables: ClassVar[dict] = {'resnames':'resname','chains':'chainID','segtypes':'segtype'}
    """
    Attributes that can be excluded from the asymmetric unit based on user specifications.
    Key is the YAML key used the input file and value is the attribute of the residue to check against.
    """

    _atom_excludables: ClassVar[dict] ={'altlocs':'altloc'}
    """
    Attributes that can be excluded from the asymmetric unit based on user specifications.
    Key is the YAML key used the input file and value is the attribute of the atom to check against.
    """

    _parent_molecule: ClassVar[object] = None
    """
    The parent molecule of the asymmetric unit, if applicable.
    This is used to link the asymmetric unit to a larger molecular structure.
    """

    def describe(self):
        return f'<AsymmetricUnit: {len(self.atoms)} atoms, {len(self.residues)} residues>'

    def __repr__(self):
        return self.describe()

    def __init__(self,
                 parsed: PDBRecordDict | DataContainer = None,
                 sourcespecs: dict = {},
                 objmanager: ObjManager = ObjManager(),
                 chainIDmanager: ChainIDManager = ChainIDManager(),
                 psf: str | Path = None):

        if not parsed:
            return

        self._parent_molecule = None
        missings = ResiduePlaceholderList([])
        ssbonds = SSBondList([])
        links = LinkList([])
        ters = TerList([])
        seqadvs = SeqadvList([])
        patches = PatchList([])
        grafts = GraftList([])

        if type(parsed) == PDBRecordDict: # PDB format
            # minimal parsed has ATOMS
            model_id = sourcespecs.get('model', None)
            atoms = AtomList.from_pdb(parsed, model_id=model_id)
            ters = TerList.from_pdb(parsed, model_id=model_id)
            if len(atoms) == 0:
                raise ValueError(f'No ATOM/HETATM records found in parsed data')
            if sourcespecs.get('reserialize', False):
                atoms.reserialize()
            else:
                atoms.sort(by=['serial'])
            if ters._has_serials: atoms.reserialize()
            objmanager.ingest(ters)
            missings = ResiduePlaceholderList.from_pdb(parsed)
            ssbonds = SSBondList.from_pdb(parsed)
            seqadvs = SeqadvList.from_pdb(parsed)
            logger.debug(f'Parsed {len(atoms)} atoms, {len(missings)} missing residues, {len(ssbonds)} ssbonds, {len(seqadvs)} seqadvs, {len(ters)} ters')
            links = LinkList.from_pdb(parsed)
            links.remove_duplicates(fields=['chainID1', 'resid1', 'chainID2', 'resid2']) # some pdb files list links multiple times (2ins, i'm looking at you)
        elif type(parsed) == DataContainer: # mmCIF format
            atoms = AtomList.from_cif(parsed)
            seqadvs = SeqadvList.from_cif(parsed)
            missings = ResiduePlaceholderList.from_cif(parsed)
            ssbonds = SSBondList.from_cif(parsed)
            links = LinkList.from_cif(parsed)
        else:
            raise TypeError(f'Unknown type for parsed: {type(parsed)}; expected dict or DataContainer')

        if psf:
            self.psf = PSFContents(psf)
            assert len(self.psf.atoms) == len(atoms), f'Error: psf file {psf} has wrong number of atoms {len(self.psf.atoms)}, expected {len(atoms)}'
            atoms.apply_psf_resnames(self.psf.atoms)
            # note we expect that there are NO ssbonds or links in the pdb file
            # if we are also using a psf file
            ssbonds.extend(self.psf.ssbonds)
            links.extend(self.psf.links)
            if len(self.psf.ssbonds) > 0:
                logger.debug(f'PSF file {psf} identifies {len(self.psf.ssbonds)} ssbonds; total ssbonds now {len(ssbonds)}')
            if len(self.psf.links) > 0:
                logger.debug(f'PSF file {psf} identifies {len(self.psf.links)} links; total links now {len(links)}')

        # at this point the objmanager is holding all mods that are in 
        # the input file but NOT in the PDB/mmCIF/psf file.
        seqmods = objmanager.get('seq', {})
        logger.debug(f'Seqmods: {seqmods}')
        topomods = objmanager.get('topol', {})
        grafts: GraftList = seqmods.get('grafts', GraftList([]))

        userlinks = topomods.get('links', LinkList([]))
        links.extend(userlinks)

        userssbonds = topomods.get('ssbonds', SSBondList([]))
        ssbonds.extend(userssbonds)

        userpatches = topomods.get('patches', PatchList([]))
        patches.extend(userpatches)

        # Build the list of residues
        excludes = sourcespecs.get('exclude', {})
        thru_dict = {resattr: excludes.get(yaml, []) for yaml, resattr in self._atom_excludables.items()}
        logger.debug(f'Atom exclusions: {thru_dict}')
        ignored_atoms = atoms.prune_exclusions(**thru_dict)
        logger.debug(f'{len(ignored_atoms)} atoms excluded by user-specified exclusions')
        fromAtoms = ResidueList.from_atomlist(atoms)
        fromResiduePlaceholders = ResidueList.from_ResiduePlaceholderlist(missings)
        residues: ResidueList = fromAtoms + fromResiduePlaceholders
        # for r in residues:
        #     logger.debug(f'Residue {r.resname}{r.resid.resid} with {len(r.atoms)} atoms (resolved: {r.resolved})')
        if sourcespecs.get('cif_residue_map_file', ''):
            write_residue_map(residues.cif_residue_map(), sourcespecs['cif_residue_map_file'])
        residues.apply_segtypes()
        # apply seqmods
        if 'deletions' in seqmods:
            residues.deletion(seqmods['deletions'])
        if 'substitutions' in seqmods:
            logger.debug(f'Applying {len(seqmods["substitutions"])} substitutions...')
            new_seqadv, wearegone = residues.substitutions(seqmods['substitutions'])
            seqadvs.extend(new_seqadv)
            logger.debug(f'There are now {len(seqadvs)} seqadv elements in seqadvs ({type(seqadvs)})')

        uniques = residues.uniqattrs(['segtype'], with_counts=True)
        nResolved = sum([1 for x in residues if x.resolved])
        nUnresolved = sum([1 for x in residues if not x.resolved])
        logger.debug(f'{len(residues)} total residues: {nResolved} resolved and {nUnresolved} unresolved')
        if 'deletions' in seqmods:
            nResolvedDelete = len(fromAtoms) - nResolved
            nUnresolvedDelete = len(fromResiduePlaceholders) - nUnresolved
            logger.debug(f'Deletions removed {nResolvedDelete} resolved and {nUnresolvedDelete} unresolved residues')
        if 'insertions' in seqmods:
            residues.apply_insertions(seqmods['insertions'])
        logger.debug(f'Renumbering residues')
        residues.renumber(links)

        logger.debug(f'Segtypes present: {uniques["segtype"]}')
        # Delete any residues dictated by user-specified exclusions
        thru_dict = {resattr: excludes.get(yaml, []) for yaml, resattr in self._residue_excludables.items()}
        logger.debug(f'Exclusions: {thru_dict}')
        # delete residues that are in user-specified exclusions
        ignored_residues = residues.prune_exclusions(**thru_dict)
        # populates a specific mapping of PDB chainID to CIF label_asym_id
        # This is only meaningful if mmCIF input is used
        residues.map_chainIDs_label_to_auth()

        ignored_seqadvs = seqadvs.assign_residues(residues)
        logger.debug(f'{len(seqadvs)} seqadvs after assign_residues')
        for s in seqadvs:
            logger.debug(f'{repr(s)}')
        ignored_ssbonds = ssbonds.assign_residues(residues)
        ignored_patches = patches.assign_residues(residues)
        logger.debug(f'{len(ssbonds)} ssbonds after assign_residues')
        for b in ssbonds:
            logger.debug(f'{b}')
        logger.debug(f'{len(patches)} patches after assign_residues')
        for p in patches:
            logger.debug(f'{p}')
        more_ignored_residues, ignored_links = links.assign_residues(residues)
        logger.debug(f'{len(links)} links after assign_residues')
        for l in links:
            logger.debug(f'{l}')
        # a deleted link may create a "free" glycan; in this case
        # we should also delete its residues; 
        ignored_residues.extend(more_ignored_residues)
        ignored_grafts, more_ignored_residues, new_ignored_links = grafts.assign_residues(residues, links)
        ignored_residues.extend(more_ignored_residues)
        ignored_links.extend(new_ignored_links)

        if excludes or len(ignored_residues + ignored_seqadvs + ignored_ssbonds + ignored_links + ignored_grafts) > 0:
            logger.debug(f'Exclusions result in deletion of:')
            logger.debug(f'    {len(ignored_residues)} residues, {len(residues)} remain')
            logger.debug(f'    {len(ignored_seqadvs)} seqadvs, {len(seqadvs)} remain')
            logger.debug(f'    {len(ignored_ssbonds)} ssbonds, {len(ssbonds)} remain')
            logger.debug(f'    {len(ignored_patches)} patches, {len(patches)} remain')
            logger.debug(f'    {len(ignored_links)} links, {len(links)} remain')
            logger.debug(f'    {len(ignored_grafts)} grafts, {len(grafts)} remain')

        logger.debug(f'{len(residues)} residues after assign_residues')

        # provide specifications of how to handle sequence issues
        # implied by PDB input
        seq_specs = sourcespecs.get('sequence', {})
        segments = SegmentList.generate(residues=residues, seq_spec=seq_specs, chainIDmanager=chainIDmanager)
        # this may have altered chainIDs for some residues.  So we must
        # be sure all mods that are chainID-specific but
        # not inheritable by segments are updated
        # for s in seqadvs:
        #     if type(s.residue) == ResidueList:
        #         logger.debug(f'{str(s)}')
        seqadvs.update_attr_from_obj_attr('chainID', 'residue', 'chainID')
        seqadvs.map_attr('dbRes', 'dbRes', {'HIS': 'HSD'})
        ssbonds.update_attr_from_obj_attr('chainID1', 'residue1', 'chainID')
        ssbonds.update_attr_from_obj_attr('chainID2', 'residue2', 'chainID')
        links.update_attr_from_obj_attr('chainID1', 'residue1', 'chainID')
        links.update_attr_from_obj_attr('chainID2', 'residue2', 'chainID')
        grafts.update_attr_from_objlist_elem_attr('chainID', 'residues', 0, 'chainID')
        logger.debug(f'Segnames in A.U.: {",".join(segments._segnames)}')
        if segments._daughters:
            logger.debug(f'Daughter chains generated: {segments._daughters}')
        logger.debug(f'Used chainIDs {chainIDmanager.Used}')
        # promote sequence numbers in any grafts to avoid collisions
        # and include the graft links in the overall list of links
        next_resid = max([x.resid for x in residues]) + 1
        for g in grafts:
            next_resid = g.set_internal_resids(next_resid)
            links.extend(g.my_links)

        # at this point, we have built the asymmetric unit according to the intention of the 
        # author of the structure AND the intention of the user in excluding certain parts
        # of that structure (ligands, ions, chains, etc).  At this point, we should apply
        # any user-defined modifications

        # First, scan all seqadv's for relevant mutations to apply
        mutations = MutationList([])
        if seq_specs.get('fix_engineered_mutations', False):
            mutations.extend(MutationList([Mutation(s) for s in seqadvs if s.typekey == 'engineered mutation']))
        if seq_specs.get('fix_conflicts', False):
            mutations.extend(MutationList([Mutation(s) for s in seqadvs if s.typekey == 'conflict']))
        mutations.extend(MutationList([Mutation(s) for s in seqadvs if s.typekey == 'user']))
        # Now append these to the objmanager's mutations
        mutations = objmanager.ingest(mutations)
        logger.debug(f'All mutations')
        mutations.sort(by=['typekey'])
        for m in mutations:
            logger.debug(str(m) + ' ' + m.typekey)
        pruned_ssbonds = ssbonds.prune_mutations(mutations)
        pruned = pruned_by_links = links.prune_mutations(mutations, segments)
        pruned_segments = pruned.get('segments', [])
        for s in pruned_segments:
            logger.debug(f'Unregistering chainID {s.segname} because this entire segment was pruned')
            chainIDmanager.unregister_chain(s.segname)

        # Now any explicitly deleted ssbonds
        topomods = objmanager.get('topol', {})
        if 'ssbondsdelete' in topomods:
            for s in ssbonds:
                if topomods['ssbondsdelete'].is_deleted(s):
                    ignored_ssbonds.append(ssbonds.remove(s))

        # finalize the objmanager
        ssbonds = objmanager.ingest(ssbonds, overwrite=True)
        links = objmanager.ingest(links, overwrite=True)
        grafts = objmanager.ingest(grafts, overwrite=True)
        patches = objmanager.ingest(patches, overwrite=True)

        segments.inherit_objs(objmanager)

        self.atoms = atoms
        self.residues = residues
        self.segments = segments
        self.objmanager = objmanager
        self.chainIDmanager = chainIDmanager
        self.psf = psf
        self.ignored = Namespace(residues=ignored_residues, links=ignored_links, ssbonds=ignored_ssbonds, seqadv=ignored_seqadvs)
        self.pruned = Namespace(ssbonds=pruned_ssbonds, residues=pruned_by_links['residues'], links=pruned_by_links['links'], segments=pruned_by_links['segments'])

    @property
    def parent_molecule(self):
        """
        The parent molecule of the asymmetric unit, if applicable.
        This is used to link the asymmetric unit to a larger molecular structure.
        """
        return self._parent_molecule
    
    def set_parent_molecule(self, parent_molecule):
        """
        Sets the parent molecule of the asymmetric unit.

        Parameters
        ----------
        parent_molecule : Any
            The parent molecule to set for the asymmetric unit.
        """
        if self._parent_molecule is not None:
            raise RuntimeError("Parent molecule is already set and cannot be changed.")
        self._parent_molecule = parent_molecule
        for segment in self.segments:
            segment.set_parent_molecule(parent_molecule)

    def add_segment(self, seg):
        """
        Adds a segment to the asymmetric unit.

        Parameters
        ----------
        seg : Segment
            The segment to add to the asymmetric unit.
        """
        self.segments.append(seg)

    def ingest_grafts(self, grafts: GraftList, residues: ResidueList, links: LinkList):
        """
        Ingests grafts into the asymmetric unit.

        Parameters
        ----------
        grafts : GraftList
            The list of grafts to ingest.
        residues : ResidueList
            The list of residues in the asymmetric unit.
        links : LinkList
            The list of links in the asymmetric unit.
        """
        for g in grafts:
            graft_chainID = g.chainID
            chain_residues = residues.filter(chainID=graft_chainID)
            # last_chain_residue_idx=residues.index(chain_residues[-1])
            next_available_resid = max([x.resid for x in chain_residues]) + 1
            g_topomods = g.source_molecule.objmanager.get('topol', {})
            g_links = g_topomods.get('links', LinkList([]))
            g.source_seg = g.source_molecule.asymmetric_unit.segments.get(segname=g.source_chainID)
            # build the list of new residues this graft contributes; the index residue is not included!
            g.my_residues = ResidueList([])
            for residue in g.source_seg.residues:
                if residue > f'{g.source_resid1}' and residue <= f'{g.source_resid2}':
                    residue.set_resid(next_available_resid)
                    residue.set_chainID(graft_chainID)
                    next_available_resid += 1
                    g.my_residues.append(residue)
            # only ingest links that are internal to this set of residues
            ingested_links = 0
            for l in g_links:
                if l.residue1 in g.my_residues and l.residue2 in g.my_residues:
                    links.append(l)
                    ingested_links += 1
            logger.debug(f'Chain {graft_chainID} of raw asymmetric unit ingests {len(g.my_residues)} residues and {ingested_links} links from graft {g.id}')
                # residues.insert(last_chain_residue_idx+1,r)
                # last_chain_residue_idx+=1

    def set_coords(self, altstruct: PDBRecordDict):
        """
        Sets the coordinates of the asymmetric unit from an alternative structure.
        
        Parameters
        ----------
        altstruct : dict
            The alternative structure containing atomic coordinates.
        """
        atoms = AtomList([Atom(p) for p in altstruct[Atom.PDB_keyword]])
        atoms.extend([Hetatm(p) for p in altstruct.get(Hetatm.PDB_keyword, [])])
        # ters = TerList([Ter(p) for p in altstruct.get(Ter.PDB_keyword, [])])
        atoms.reserialize()
        # atoms.adjustSerials(ters)
        altRes = ResidueList(atoms)
        overwrites = []
        for ar in altRes:
            tr: Residue = self.residues.get(resname=ar.resname, resid=ar.resid, chainID=ar.chainID)
            if tr:
                tr.atoms.overwrite_positions(ar.atoms)
                overwrites.append(ar)
        logger.debug(f'set_coords: {len(overwrites)} residues overwritten')

