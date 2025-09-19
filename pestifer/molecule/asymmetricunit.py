# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for building the asymmetric unit from a PDB file
"""
from __future__ import annotations

import logging

from argparse import Namespace
from mmcif.api.PdbxContainers import DataContainer
from pathlib import Path
from pidibble.pdbrecord import PDBRecordDict
from typing import ClassVar, TYPE_CHECKING

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
from ..util.stringthings import parse_filter_expression

if TYPE_CHECKING:
    from .molecule import Molecule

logger = logging.getLogger(__name__)

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

        self.parent_molecule: Molecule = None

        missings = ResiduePlaceholderList([])
        ssbonds = SSBondList([])
        links = LinkList([])
        ters = TerList([])
        seqadvs = SeqadvList([])
        patches = PatchList([])
        grafts = GraftList([])
    
        # These will be built below
        self.atoms = AtomList([])
        self.residues = ResidueList([])
        self.segments = SegmentList([])
        self.objmanager = objmanager
        self.chainIDmanager = chainIDmanager
        self.psfcontents = PSFContents(psf) if psf else None
        self.ignored = Namespace()
        self.pruned = Namespace()

        atom_include_logic = sourcespecs.get('include', [])
        atom_exclude_logic = sourcespecs.get('exclude', [])
        if len(atom_exclude_logic) > 0 and len(atom_include_logic) > 0:
            raise ValueError('Cannot specify both include and exclude logic for atoms')

        if not parsed:
            # return an empty unconstructed asymmetric unit
            return

        if isinstance(parsed, PDBRecordDict): # PDB format
            # minimal parsed has ATOMS
            model_id = sourcespecs.get('model', None)
            atoms = AtomList.from_pdb(parsed, model_id=model_id)
            ters = TerList.from_pdb(parsed, model_id=model_id)
            if len(atoms) == 0:
                raise ValueError(f'No ATOM/HETATM records found in parsed data')
            if sourcespecs.get('reserialize', False) or self.psfcontents is not None:
                atoms.reserialize()
            else:
                atoms.sort(by=['serial'])
            if ters.has_serials(): atoms.reserialize()
            objmanager.ingest(ters)
            missings = ResiduePlaceholderList.from_pdb(parsed)
            ssbonds = SSBondList.from_pdb(parsed)
            seqadvs = SeqadvList.from_pdb(parsed)
            logger.debug(f'Parsed {len(atoms)} atoms, {len(missings)} missing residues, {len(ssbonds)} ssbonds, {len(seqadvs)} seqadvs, {len(ters)} ters')
            links = LinkList.from_pdb(parsed)
            links.remove_duplicates(fields=['chainID1', 'resid1', 'chainID2', 'resid2']) # some pdb files list links multiple times (2ins, i'm looking at you)
        elif isinstance(parsed, DataContainer): # mmCIF format
            atoms = AtomList.from_cif(parsed)
            seqadvs = SeqadvList.from_cif(parsed)
            missings = ResiduePlaceholderList.from_cif(parsed)
            ssbonds = SSBondList.from_cif(parsed)
            links = LinkList.from_cif(parsed)
        else:
            raise TypeError(f'Cannot construct Asymmetric Unit from data of type {type(parsed)}; expected PDBRecordDict or DataContainer')

        uniq_attr = atoms.uniqattrs(['chainID'])
        logger.debug(f'ChainIDs in structure: {uniq_attr["chainID"]}')
        logger.debug(f'Samples of first atoms in each chainID:')
        for chainID in uniq_attr["chainID"]:
            first_atom = atoms.get(lambda x: x.chainID == chainID)[0]
            logger.debug(f'  {chainID}: {first_atom.resname} {first_atom.resid.resid}')

        if self.psfcontents is not None:
            assert len(self.psfcontents.atoms) == len(atoms), f'Error: psf file {psf} has wrong number of atoms {len(self.psfcontents.atoms)}, expected {len(atoms)}'
            logger.debug(f'PSF defines {len(self.psfcontents.atoms)} atoms')
            logger.debug(f'PSF defines {len(self.psfcontents.segments)} segments: {self.psfcontents.segnames}')
            logger.debug(f'Samples of first atom in each segment:')
            for segname in self.psfcontents.segnames:
                first_atom = self.psfcontents.atoms.get(lambda x: x.segname == segname)[0]
                logger.debug(f'  {segname}: {first_atom.resname} {first_atom.resid.resid}')
            atoms.apply_psf_attributes(self.psfcontents.atoms)
            # note we expect that there are NO ssbonds or links in the pdb file
            # because it is not a Structure File from the PDB (since it has a companion PSF)
            ssbonds.extend(self.psfcontents.ssbonds)
            links.extend(self.psfcontents.links)
            if len(self.psfcontents.ssbonds) > 0:
                logger.debug(f'PSF file {psf} identifies {len(self.psfcontents.ssbonds)} ssbonds; total ssbonds now {len(ssbonds)}')
            if len(self.psfcontents.links) > 0:
                logger.debug(f'PSF file {psf} identifies {len(self.psfcontents.links)} links; total links now {len(links)}')

        # at this point the objmanager is holding all mods that are in 
        # the input file but NOT in the PDB/mmCIF/psf file.
        seqmods = objmanager.get('seq', {})
        for k,v in seqmods.items():
            if k != 'grafts':
                if len(v) > 0:
                    logger.debug(f'Seqmod {k}')
                    for mod in v:
                        logger.debug(f'  {str(mod)}')
        topomods = objmanager.get('topol', {})
        grafts: GraftList = seqmods.get('grafts', GraftList([]))

        userlinks = topomods.get('links', LinkList([]))
        links.extend(userlinks)

        userssbonds = topomods.get('ssbonds', SSBondList([]))
        for u in userssbonds:
            logger.debug(f'Adding user-specified ssbond {u}')
        ssbonds.extend(userssbonds)

        userpatches = seqmods.get('patches', PatchList([]))
        for u in userpatches:
            logger.debug(f'Adding user-specified patch {u}')
        patches.extend(userpatches)

        # Build the list of residues
        ignored_atom_count = 0
        ignored_missing_residue_count = 0
        if atom_include_logic or atom_exclude_logic:
            logger.debug(f'Applying atom logics')
            if atom_include_logic:
                logger.debug(f'Atoms: Applying atom inclusion logic: {atom_include_logic}')
                ignored_atom_count += atoms.apply_inclusion_logics(atom_include_logic)
            elif atom_exclude_logic:
                logger.debug(f'Atoms: Applying atom exclusion logic: {atom_exclude_logic}')
                ignored_atom_count += atoms.apply_exclusion_logics(atom_exclude_logic)
            logger.debug(f'Ignored {ignored_atom_count} atoms from PDB by inclusion/exclusion logic')
            if atom_include_logic:
                logger.debug(f'Missing residues: Applying atom inclusion logic: {atom_include_logic}')
                ignored_missing_residue_count += missings.apply_inclusion_logics(atom_include_logic)
            elif atom_exclude_logic:
                logger.debug(f'Missing residues: Applying atom exclusion logic: {atom_exclude_logic}')
                ignored_missing_residue_count += missings.apply_exclusion_logics(atom_exclude_logic)
            logger.debug(f'Ignored {ignored_missing_residue_count} missing residues from PDB by inclusion/exclusion logic')
        if self.psfcontents is not None:
            ignored_psfatom_count = self.psfcontents.apply_atom_logics(atom_include_logic, atom_exclude_logic)
            assert len(atoms) == len(self.psfcontents.atoms), f'Atom logic is not consistent between PDB and PSF'
            assert ignored_psfatom_count == ignored_atom_count, f'Ignored atom count is not consistent between PDB and PSF'
            logger.debug(f'Ignored {ignored_psfatom_count} atoms from PSF by inclusion/exclusion logic')
        fromAtoms = ResidueList.from_residuegrouped_atomlist(atoms)
        fromResiduePlaceholders = ResidueList.from_ResiduePlaceholderlist(missings)
        if self.psfcontents is not None:
            assert len(fromResiduePlaceholders) == 0, f'Expected no missing residues because you specified a companion PSF, but found {len(fromResiduePlaceholders)}'
        residues: ResidueList = fromAtoms + fromResiduePlaceholders
            
        residues.apply_segtypes()
        # for r in residues:
        #     logger.debug(str(r))
        # apply seqmods
        if 'deletions' in seqmods:
            residues.deletion(seqmods['deletions'])
        if 'substitutions' in seqmods and len(seqmods['substitutions']) > 0:
            logger.debug(f'Applying {len(seqmods["substitutions"])} substitutions...')
            new_seqadv, wearegone = residues.substitutions(seqmods['substitutions'])
            seqadvs.extend(new_seqadv)
            logger.debug(f'There are now {len(seqadvs)} seqadv elements in seqadvs ({type(seqadvs)})')

        uniques = residues.uniqattrs(['segtype'], with_counts=True)
        nResolved = sum([1 for x in residues if x.resolved])
        nUnresolved = sum([1 for x in residues if not x.resolved])
        logger.debug(f'{len(residues)} total residues: {nResolved} resolved and {nUnresolved} unresolved')
        if 'deletions' in seqmods and len(seqmods['deletions']) > 0:
            nResolvedDelete = len(fromAtoms) - nResolved
            nUnresolvedDelete = len(fromResiduePlaceholders) - nUnresolved
            logger.debug(f'Deletions removed {nResolvedDelete} resolved and {nUnresolvedDelete} unresolved residues')
        if 'insertions' in seqmods and len(seqmods['insertions']) > 0:
            residues.apply_insertions(seqmods['insertions'])
        logger.debug(f'Renumbering residues')
        residues.renumber(links)

        logger.debug(f'Segtypes present: {uniques["segtype"]}')
        # Delete any residues dictated by user-specified exclusions
        # thru_dict = {resattr: excludes.get(yaml, []) for yaml, resattr in self._residue_excludables.items()}
        # logger.debug(f'Exclusions: {thru_dict}')
        # # delete residues that are in user-specified exclusions
        # ignored_residues = residues.prune_exclusions(residues.dict_to_condition(thru_dict))
        # populates a specific mapping of PDB chainID to CIF label_asym_id
        # This is only meaningful if mmCIF input is used
        residues.map_chainIDs_label_to_auth()

        ignored_seqadvs = seqadvs.assign_residues(residues)
        logger.debug(f'{len(seqadvs)} seqadvs after assign_residues')
        # for s in seqadvs:
        #     logger.debug(f'{str(s)}')
        ignored_ssbonds = ssbonds.assign_residues(residues)
        logger.debug(f'{len(patches)} patches before assign_residues')
        ignored_patches = patches.assign_residues(residues)
        logger.debug(f'{len(ssbonds)} ssbonds after assign_residues')
        # for b in ssbonds:
        #     logger.debug(f'{str(b)}')
        logger.debug(f'{len(patches)} patches after assign_residues')
        # for p in patches:
        #     logger.debug(f'{str(p)}')
        more_ignored_residues, ignored_links = links.assign_residues(residues)
        logger.debug(f'{len(links)} links after assign_residues')
        # for l in links:
        #     logger.debug(f'{str(l)}')
        ignored_residues = ResidueList([])
        ignored_residues.extend(more_ignored_residues)
        ignored_grafts = grafts.assign_residues(residues, links)
        ignored_residues.extend(more_ignored_residues)
        # ignored_links.extend(new_ignored_links)
        total_ignored_residue_count = len(ignored_residues) + ignored_missing_residue_count
        if (ignored_atom_count + total_ignored_residue_count + len(ignored_seqadvs) + len(ignored_ssbonds) + len(ignored_links) + len(ignored_grafts)) > 0:
            logger.debug(f'Inclusion/exclusion logic results in deletion of:')
            logger.debug(f'    {ignored_atom_count} atoms, {len(atoms)} remain')
            logger.debug(f'    {total_ignored_residue_count} residues, {len(residues)} remain')
            logger.debug(f'    {len(ignored_seqadvs)} seqadvs, {len(seqadvs)} remain')
            logger.debug(f'    {len(ignored_ssbonds)} ssbonds, {len(ssbonds)} remain')
            logger.debug(f'    {len(ignored_patches)} patches, {len(patches)} remain')
            logger.debug(f'    {len(ignored_links)} links, {len(links)} remain')
            logger.debug(f'    {len(ignored_grafts)} grafts, {len(grafts)} remain')

        logger.debug(f'{len(residues)} residues survived parsing')

        if self.psfcontents is not None:
            self.psfcontents.remove_ignored_residues(ignored_residues)

        if len(grafts)>0:
            logger.debug('Grafts:')
            for g in grafts:
                logger.debug(str(g))

        # provide specifications of how to handle sequence issues
        # implied by PDB input
        seq_specs = sourcespecs.get('sequence', {})
        psfcompanion = None
        if self.psfcontents:
            psfcompanion = self.psfcontents.segments
        self.segments.generate_from_residues(residues=residues, seq_spec=seq_specs, chainIDmanager=chainIDmanager, psfcompanion=psfcompanion)
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
        logger.debug(f'Segnames in A.U.: {", ".join(self.segments.segnames)}')
        if self.segments.daughters:
            logger.debug(f'Daughter chains generated: {self.segments.daughters}')
        logger.debug(f'Used chainIDs {chainIDmanager.Used}')
        # promote sequence numbers in any grafts to avoid collisions
        # and include the graft links in the overall list of links
        next_resid = max([x.resid for x in residues]) + 1
        for g in grafts.data:
            next_resid = g.set_internal_resids(next_resid)
            links.extend(g.donor_internal_links)
            links.extend(g.donor_external_links)

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
        objmanager.ingest(mutations)
        mutations = objmanager.get('seq', {}).get('mutations', MutationList([]))
        if len(mutations) > 0:
            logger.debug(f'All mutations')
            mutations.sort(by=['typekey'])
            for m in mutations:
                logger.debug(str(m) + ' ' + m.typekey)
        pruned_objects = self.segments.prune_topology(mutations, links, ssbonds)
                # pruned_ssbonds = ssbonds.prune_mutations(mutations)
                # pruned = pruned_by_links = links.prune_mutations(mutations, segments)
        # logger.debug(f'Pruned objects: {pruned_objects}')
        pruned_segments: SegmentList = pruned_objects.get('segments', [])
        for s in pruned_segments.data:
            logger.debug(f'Unregistering chainID {s.segname} because this entire segment was pruned')
            chainIDmanager.unregister_chain(s.segname)

        # Now any explicitly deleted ssbonds
        topomods = objmanager.get('topol', {})
        if 'ssbondsdelete' in topomods and len(topomods['ssbondsdelete']) > 0:
            logger.debug(f'SSBonds to delete: {topomods["ssbondsdelete"]} ')
            for s in ssbonds:
                logger.debug(f'Examining ssbond {s}')
                if topomods['ssbondsdelete'].is_deleted(s):
                    logger.debug(f'Deleting ssbond {s}')
                    ignored_ssbonds.append(ssbonds.remove_and_return(s))

        # finalize the objmanager
        for olist in [ssbonds, links, grafts, patches]:
            objmanager.ingest(olist, overwrite=True)

        self.segments.inherit_objs(objmanager)

        self.atoms = atoms
        self.residues = residues
        self.objmanager = objmanager
        self.chainIDmanager = chainIDmanager
        self.psf = psf
        self.ignored = Namespace(residues=ignored_residues, links=ignored_links, ssbonds=ignored_ssbonds, seqadv=ignored_seqadvs)
        self.pruned = Namespace(**pruned_objects)
    
    def set_parent_molecule(self, parent_molecule):
        """
        Sets the parent molecule of the asymmetric unit.

        Parameters
        ----------
        parent_molecule : Any
            The parent molecule to set for the asymmetric unit.
        """
        if self.parent_molecule is not None:
            raise RuntimeError("Parent molecule is already set and cannot be changed.")
        self.parent_molecule = parent_molecule
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
        for g in grafts.data:
            graft_chainID = g.chainID
            chain_residues = residues.newfilter(lambda x: x.chainID == graft_chainID)
            # last_chain_residue_idx=residues.index(chain_residues[-1])
            next_available_resid = max([x.resid for x in chain_residues]) + 1
            g_topomods = g.source_molecule.objmanager.get('topol', {})
            g_links = g_topomods.get('links', LinkList([]))
            g.source_seg = g.source_molecule.asymmetric_unit.segments.get(lambda x: x.segname == g.source_chainID)
            # build the list of new residues this graft contributes; the index residue is not included!
            g.my_residues = ResidueList([])
            for residue in g.source_seg.residues.data:
                if residue > f'{g.source_resid1}' and residue <= f'{g.source_resid2}':
                    residue.set_resid(next_available_resid)
                    residue.set_chainID(graft_chainID)
                    next_available_resid += 1
                    g.my_residues.append(residue)
            # only ingest links that are internal to this set of residues
            ingested_links = 0
            for l in g_links.data:
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
            tr: Residue = self.residues.get(lambda x: x.resname == ar.resname and x.resid == ar.resid and x.chainID == ar.chainID)
            if tr:
                tr.atoms.overwrite_positions(ar.atoms)
                overwrites.append(ar)
        logger.debug(f'set_coords: {len(overwrites)} residues overwritten')

