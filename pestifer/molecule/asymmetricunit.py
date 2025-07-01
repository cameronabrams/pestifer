# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for building the asymmetric unit from a PDB file
"""
import logging
logger=logging.getLogger(__name__)

from argparse import Namespace
from mmcif.api.PdbxContainers import DataContainer

from .atom import AtomList, Atom, Hetatm
from ..core.baseobj import AncestorAwareObj
from ..objs.graft import Graft, GraftList
from ..objs.link import Link, LinkList
from ..objs.mutation import Mutation, MutationList
from ..objs.seqadv import Seqadv, SeqadvList
from ..objs.ssbond import SSBond, SSBondList
from ..objs.patch import Patch, PatchList
from ..objs.ter import Ter, TerList
from ..core.objmanager import ObjManager, ObjCats
from ..psfutil.psfcontents import PSFContents
from .residue import ResidueList,EmptyResidue,EmptyResidueList
from .segment import SegmentList
from ..core.stringthings import my_logger
from ..util.util import write_residue_map
from ..util.cifutil import CIFdict

class AsymmetricUnit(AncestorAwareObj):
    """
    A class for building the asymmetric unit from a PDB file.
    """

    req_attr=AncestorAwareObj.req_attr+['atoms','residues','objmanager']
    """
    Required attributes:
    
    - atoms: AtomList of all atoms in the asymmetric unit
    - residues: ResidueList of all residues in the asymmetric unit
    - objmanager: ObjManager containing all modifications and additional objects
    """

    opt_attr=AncestorAwareObj.opt_attr+['segments','ignored','pruned']
    """
    Optional attributes:

    - segments: SegmentList of segments in the asymmetric unit
    - ignored: List of residues/atoms to ignore
    - pruned: List of residues/atoms that have been pruned
    """

    residue_excludables={'resnames':'resname','chains':'chainID','segtypes':'segtype'}
    """
    Attributes that can be excluded from the asymmetric unit based on user specifications.
    Key is the YAML key used the input file and value is the attribute of the residue to check against.
    """

    atom_excludables={'altlocs':'altloc'}
    """
    Attributes that can be excluded from the asymmetric unit based on user specifications.
    Key is the YAML key used the input file and value is the attribute of the atom to check against.
    """

    def __init__(self,**objs):
        if len(objs)==0:
            logger.debug('Generating an empty A.U.')
            input_dict={
                'atoms':AtomList([]),
                'residues':ResidueList([]),
                'objmanager':ObjManager(),
                'segments':SegmentList({},ResidueList([]))
            }
        else: #
            pr=objs.get('parsed',None)
            sourcespecs=objs.get('sourcespecs',{})
            excludes=sourcespecs.get('exclude',{})
            objmanager=objs.get('objmanager',ObjManager())
            chainIDmanager=objs.get('chainIDmanager',None)
            psf=objs.get('psf',None)

            missings=EmptyResidueList([])
            ssbonds=SSBondList([])
            links=LinkList([])
            ters=TerList([])
            seqadvs=SeqadvList([])
            patches=PatchList([])
            if type(pr)==dict: # PDB format
                # minimal pr has ATOMS
                if 'model' in sourcespecs:
                    atoms=AtomList([Atom(p) for p in pr[Atom.PDB_keyword] if p.model==sourcespecs['model']])
                    atoms.extend([Hetatm(p) for p in pr.get(Hetatm.PDB_keyword,[]) if p.model==sourcespecs['model']])
                    ters=TerList([Ter(p) for  p in pr.get(Ter.PDB_keyword,[]) if p.model==sourcespecs['model']])
                else:
                    atoms=AtomList([Atom(p) for p in pr[Atom.PDB_keyword]])
                    atoms.extend([Hetatm(p) for p in pr.get(Hetatm.PDB_keyword,[])])
                    ters=TerList([Ter(p) for  p in pr.get(Ter.PDB_keyword,[])])
                if sourcespecs.get('reserialize',False):
                    atoms.reserialize()
                else:
                    atoms.sort(by=['serial'])
                if len(ters)>0:
                    lnonemptyters=len([x for x in ters if x.serial!=None])
                    logger.debug(f'{lnonemptyters} TER records require adjusting atom serial numbers')
                    if lnonemptyters>0:
                        atoms.reserialize()
                # atoms.adjustSerials(ters)
                objmanager.ingest(ters)
                if 'REMARK.465' in pr:
                    missings=EmptyResidueList([EmptyResidue(p) for p in pr['REMARK.465'].tables['MISSING']])
                ssbonds=SSBondList([SSBond(p) for p in pr.get(SSBond.PDB_keyword,[])])
                logger.debug(f'PDB yields {len(ssbonds)} ssbonds')
                seqadvs=SeqadvList([Seqadv(p) for p in pr.get(Seqadv.PDB_keyword,[])])
                links=LinkList([Link(p) for p in pr.get(Link.PDB_keyword,[])])
                links.remove_duplicates(fields=['chainID1','resseqnum1','insertion1','chainID2','resseqnum2','insertion2']) # some pdb files list links multiple times (2ins, i'm looking at you)
            elif type(pr)==DataContainer: # mmCIF format
                obj=pr.getObj(Atom.mmCIF_name)
                atoms=AtomList([Atom(CIFdict(obj,i)) for i in range(len(obj))])
                obj=pr.getObj(Seqadv.mmCIF_name)
                seqadvs=SeqadvList([Seqadv(CIFdict(obj,i)) for i in range(len(obj))])
                logger.debug(f'{len(seqadvs)} {type(seqadvs)} seqadv detected; typekeys: {[s.typekey for s in seqadvs]}')
                obj=pr.getObj(EmptyResidue.mmCIF_name)
                missings=EmptyResidueList([EmptyResidue(CIFdict(obj,i)) for i in range(len(obj))])
                obj=pr.getObj('struct_conn')
                ssdicts=[]
                lndicts=[]
                for i in range(len(obj)):
                    conn_type_id=obj.getValue('conn_type_id',i)
                    if conn_type_id=='disulf':
                        ssdicts.append(CIFdict(obj,i))
                    elif conn_type_id in ['covale','metalc']:
                        lndicts.append(CIFdict(obj,i))
                ssbonds=SSBondList([SSBond(x) for x in ssdicts])
                links=LinkList([Link(x) for x in lndicts])

            if psf:
                self.psf=PSFContents(psf)
                assert len(self.psf.atoms)==len(atoms),f'Error: psf file {psf} has wrong number of atoms {len(self.psf.atoms)}, expected {len(atoms)}'
                atoms.apply_psf_resnames(self.psf.atoms)
                # note we expect that there are NO ssbonds or links in the pdb file
                # if we are also using a psf file
                ssbonds.extend(self.psf.ssbonds)
                links.extend(self.psf.links)
                if len(self.psf.ssbonds)>0:
                    logger.debug(f'PSF file {psf} identifies {len(self.psf.ssbonds)} ssbonds; total ssbonds now {len(ssbonds)}')
                if len(self.psf.links)>0:
                    logger.debug(f'PSF file {psf} identifies {len(self.psf.links)} links; total links now {len(links)}')
            
            # at this point the objmanager is holding all mods that are in 
            # the input file but NOT in the PDB/mmCIF/psf file.
            seqmods=objmanager.get('seq',{})
            topomods=objmanager.get('topol',{})
            grafts=seqmods.get('grafts',GraftList([]))

            userlinks=topomods.get('links',LinkList([]))
            links.extend(userlinks)
        
            userssbonds=topomods.get('ssbonds',SSBondList([]))
            ssbonds.extend(userssbonds)

            userpatches=topomods.get('patches',[])
            patches.extend(userpatches)

            # Build the list of residues
            thru_dict={resattr:excludes.get(yaml,[]) for yaml,resattr in self.atom_excludables.items()}
            logger.debug(f'Atom exclusions: {thru_dict}')
            ignored_atoms=atoms.prune_exclusions(**thru_dict)
            logger.debug(f'{len(ignored_atoms)} atoms excluded by user-specified exclusions')
            fromAtoms=ResidueList(atoms)
            fromEmptyResidues=ResidueList(missings)
            residues=fromAtoms+fromEmptyResidues
            # self.ingest_grafts(grafts,residues,links)
            if sourcespecs.get('cif_residue_map_file',''):
                write_residue_map(residues.cif_residue_map(),sourcespecs['cif_residue_map_file'])
            residues.apply_segtypes()
            # apply seqmods
            if 'deletions' in seqmods:
                residues.deletion(seqmods['deletions'])
            if 'substitutions' in seqmods:
                logger.debug(f'Applying {len(seqmods["substitutions"])} substitutions...')
                new_seqadv,wearegone=residues.substitutions(seqmods['substitutions'])
                seqadvs.extend(new_seqadv)
                logger.debug(f'There are now {len(seqadvs)} seqadv elements in seqadvs ({type(seqadvs)})')

            uniques=residues.uniqattrs(['segtype'],with_counts=True)
            nResolved=sum([1 for x in residues if x.resolved])
            nUnresolved=sum([1 for x in residues if not x.resolved])
            logger.debug(f'{len(residues)} total residues: {nResolved} resolved and {nUnresolved} unresolved')
            if 'deletions' in seqmods:
                nResolvedDelete=len(fromAtoms)-nResolved
                nUnresolvedDelete=len(fromEmptyResidues)-nUnresolved
                logger.debug(f'Deletions removed {nResolvedDelete} resolved and {nUnresolvedDelete} unresolved residues')
            if 'insertions' in seqmods:
                residues.apply_insertions(seqmods['insertions'])
            logger.debug(f'Renumbering residues')
            residues.renumber(links)

            logger.debug(f'Segtypes present: {uniques["segtype"]}')
            # Delete any residues dictated by user-specified exclusions
            thru_dict={resattr:excludes.get(yaml,[]) for yaml,resattr in self.residue_excludables.items()}
            logger.debug(f'Exclusions: {thru_dict}')
            # delete residues that are in user-specified exclusions
            ignored_residues=residues.prune_exclusions(**thru_dict)
            # populates a specific mapping of PDB chainID to CIF label_asym_id
            # This is only meaningful if mmCIF input is used
            residues.map_chainIDs_label_to_auth()
            logger.debug(f'{len(residues)} residues before assign_residues')
            for r in residues:
                if hasattr(r,'label_seq_id'):
                    logger.debug(f'{r.chainID}_{r.resname}{r.resseqnum}')
                else:
                    logger.debug(f'{r.chainID}_{r.resname}{r.resseqnum}{r.insertion}')

            # Give each Seqadv a residue identifier
            ignored_seqadvs=seqadvs.assign_residues(residues)
            ignored_ssbonds=ssbonds.assign_residues(residues)
            ignored_patches=patches.assign_residues(residues)
            logger.debug(f'{len(ssbonds)} ssbonds after assign_residues')
            for b in ssbonds:
                logger.debug(f'{b}')
            logger.debug(f'{len(patches)} patches after assign_residues')
            for p in patches:
                logger.debug(f'{p}')
            more_ignored_residues,ignored_links=links.assign_residues(residues)
            logger.debug(f'{len(links)} links after assign_residues')
            for l in links:
                logger.debug(f'{l}')
            # a deleted link may create a "free" glycan; in this case
            # we should also delete its residues; problem is that 
            ignored_residues.extend(more_ignored_residues)
            ignored_grafts,more_ignored_residues,new_ignored_links=grafts.assign_residues(residues,links)
            ignored_residues.extend(more_ignored_residues)
            ignored_links.extend(new_ignored_links)
            
            if excludes or len(ignored_residues+ignored_seqadvs+ignored_ssbonds+ignored_links+ignored_grafts)>0:
                logger.debug(f'Exclusions result in deletion of:')
                logger.debug(f'    {len(ignored_residues)} residues, {len(residues)} remain')
                logger.debug(f'    {len(ignored_seqadvs)} seqadvs, {len(seqadvs)} remain')
                logger.debug(f'    {len(ignored_ssbonds)} ssbonds, {len(ssbonds)} remain')
                logger.debug(f'    {len(ignored_patches)} patches, {len(patches)} remain')
                logger.debug(f'    {len(ignored_links)} links, {len(links)} remain')
                logger.debug(f'    {len(ignored_grafts)} grafts, {len(grafts)} remain')

            logger.debug(f'{len(residues)} residues after assign_residues')
            for r in residues:
                if hasattr(r,'label_seq_id'):
                    logger.debug(f'{r.chainID}_{r.resname}{r.resseqnum}')
                else:
                    logger.debug(f'{r.chainID}_{r.resname}{r.resseqnum}{r.insertion}')

            # provide specifications of how to handle sequence issues
            # implied by PDB input
            seq_specs=sourcespecs.get('sequence',{})
            segments=SegmentList(seq_specs,residues,chainIDmanager)
            # this may have altered chainIDs for some residues.  So we must
            # be sure all mods that are chainID-specific but
            # not inheritable by segments are updated
            for s in seqadvs:
                if type(s.residue)==ResidueList:
                    logger.debug(f'{str(s)}')
            seqadvs.update_attr_from_obj_attr('chainID','residue','chainID')
            seqadvs.map_attr('dbRes','dbRes',{'HIS':'HSD'})
            ssbonds.update_attr_from_obj_attr('chainID1','residue1','chainID')
            ssbonds.update_attr_from_obj_attr('chainID2','residue2','chainID')
            links.update_attr_from_obj_attr('chainID1','residue1','chainID')
            links.update_attr_from_obj_attr('chainID2','residue2','chainID')
            grafts.update_attr_from_objlist_elem_attr('chainID','residues',0,'chainID')
            logger.debug(f'Segnames in A.U.: {",".join(segments.segnames)}')
            if segments.daughters:
                logger.debug(f'Daughter chains generated: {segments.daughters}')
            logger.debug(f'Used chainIDs {chainIDmanager.Used}')
            # promote sequence numbers in any grafts to avoid collisions
            next_resseqnum=max([x.resseqnum for x in residues])+1
            for g in grafts:
                next_resseqnum=g.set_links(next_resseqnum)
                my_logger(str(g),logger.debug,fill=' ',just='<')
                links.extend(g.my_links)

            # at this point, we have built the asymmetric unit according to the intention of the 
            # author of the structure AND the intention of the user in excluding certain parts
            # of that structure (ligands, ions, chains, etc).  At this point, we should apply
            # any user-defined modifications

            # First, scan all seqadv's for relevant mutations to apply
            mutations=MutationList([])
            if seq_specs.get('fix_engineered_mutations',False):
                mutations.extend(MutationList([Mutation(s) for s in seqadvs if s.typekey=='engineered']))
            if seq_specs.get('fix_conflicts',False):
                mutations.extend(MutationList([Mutation(s) for s in seqadvs if s.typekey=='conflict']))
            mutations.extend(MutationList([Mutation(s) for s in seqadvs if s.typekey=='user']))
            # Now append these to the objmanager's mutations
            mutations=objmanager.ingest(mutations)
            logger.debug(f'All mutations')
            mutations.sort(by=['typekey'])
            for m in mutations:
                logger.debug(str(m)+' '+m.typekey)
            pruned_ssbonds=ssbonds.prune_mutations(mutations)
            pruned=pruned_by_links=links.prune_mutations(mutations,segments)
            pruned_segments=pruned.get('segments',[])
            for s in pruned_segments:
                logger.debug(f'Unregistering chainID {s.segname} because this entire segment was pruned')
                chainIDmanager.unregister_chain(s.segname)

            # Now any explicitly deleted ssbonds
            topomods=objmanager.get('topol',{})
            if 'ssbondsdelete' in topomods:
                for s in ssbonds:
                    if topomods['ssbondsdelete'].is_deleted(s):
                        ignored_ssbonds.append(ssbonds.remove(s))

            # finalize the objmanager
            ssbonds=objmanager.ingest(ssbonds,overwrite=True)
            links=objmanager.ingest(links,overwrite=True)
            grafts=objmanager.ingest(grafts,overwrite=True)
            patches=objmanager.ingest(patches,overwrite=True)
            
            segments.inherit_objs(objmanager)

            input_dict={
                'atoms':atoms,
                'residues':residues,
                'objmanager':objmanager,
                'segments':segments,
                'ignored':Namespace(residues=ignored_residues,links=ignored_links,ssbonds=ignored_ssbonds,seqadv=ignored_seqadvs),
                'pruned':Namespace(ssbonds=pruned_ssbonds,residues=pruned_by_links['residues'],links=pruned_by_links['links'],segments=pruned_by_links['segments'])
            }
        super().__init__(input_dict)

    def add_segment(self,seg):
        """
        Adds a segment to the asymmetric unit.

        Parameters
        ----------
        seg : Segment
            The segment to add to the asymmetric unit.
        """
        self.segments.append(seg)

    def ingest_grafts(self,grafts,residues,links):
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
            graft_chainID=g.chainID
            chain_residues=residues.filter(chainID=graft_chainID)
            # last_chain_residue_idx=residues.index(chain_residues[-1])
            next_available_resid=max([x.resseqnum for x in chain_residues])+1
            g_topomods=g.source_molecule.objmanager.get('topol',{})
            g_links=g_topomods.get('links',LinkList([]))
            g.source_seg=g.source_molecule.asymmetric_unit.segments.get(segname=g.source_chainID)
            # build the list of new residues this graft contributes; the index residue is not included!
            g.my_residues=ResidueList([])
            for residue in g.source_seg.residues:
                if residue>f'{g.source_resseqnum1}{g.source_insertion1}' and residue<=f'{g.source_resseqnum2}{g.source_insertion2}':
                    residue.set_resseqnum(next_available_resid)
                    residue.set_chainID(graft_chainID)
                    next_available_resid+=1
                    g.my_residues.append(residue)
            # only ingest links that are internal to this set of residues
            ingested_links=0
            for l in g_links:
                if l.residue1 in g.my_residues and l.residue2 in g.my_residues:
                    links.append(l)
                    ingested_links+=1
            logger.debug(f'Chain {graft_chainID} of raw asymmetric unit ingests {len(g.my_residues)} residues and {ingested_links} links from graft {g.id}')
                # residues.insert(last_chain_residue_idx+1,r)
                # last_chain_residue_idx+=1

    def set_coords(self,altstruct):
        """
        Sets the coordinates of the asymmetric unit from an alternative structure.
        
        Parameters
        ----------
        altstruct : dict
            The alternative structure containing atomic coordinates.
        """
        atoms=AtomList([Atom(p) for p in altstruct[Atom.PDB_keyword]])
        atoms.extend([Hetatm(p) for p in altstruct.get(Hetatm.PDB_keyword,[])])
        ters=TerList([Ter(p) for  p in altstruct.get(Ter.PDB_keyword,[])])
        atoms.reserialize()
        # atoms.adjustSerials(ters)
        altRes=ResidueList(atoms)
        overwrites=[]
        for ar in altRes:
            tr=self.residues.get(resname=ar.resname,resseqnum=ar.resseqnum,insertion=ar.insertion,chainID=ar.chainID)
            if tr:
                tr.atoms.overwritePositions(ar.atoms)
                overwrites.append(ar)
        logger.debug(f'set_coords: {len(overwrites)} residues overwritten')

