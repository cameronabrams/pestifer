# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A class for building the asymmetric unit from a PDB file

"""
import logging
from .mods import *
from argparse import Namespace
from .atom import AtomList,Atom,Hetatm
from .residue import ResidueList,Ter,TerList,EmptyResidue,EmptyResidueList
from .segment import SegmentList
from .basemod import AncestorAwareMod
from .modmanager import ModManager
from .stringthings import my_logger
from mmcif.api.PdbxContainers import DataContainer
from .config import *
from .util import write_residue_map
from .psf import PSFContents

logger=logging.getLogger(__name__)

class AsymmetricUnit(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['atoms','residues','modmanager']
    opt_attr=AncestorAwareMod.opt_attr+['segments','ignored','pruned']
    excludables={'resnames':'resname','chains':'chainID','segtypes':'segtype'}
    def __init__(self,**objs):
        if len(objs)==0:
            logger.debug('Generating an empty A.U.')
            input_dict={
                'atoms':AtomList([]),
                'residues':ResidueList([]),
                'modmanager':ModManager(),
                'segments':SegmentList({},ResidueList([]))
            }
        else: #
            pr=objs.get('parsed',None)
            sourcespecs=objs.get('sourcespecs',{})
            modmanager=objs.get('modmanager',ModManager())
            # logger.debug(f'User mods {type(mods)} at asymmetric unit: {mods.__dict__}')
            chainIDmanager=objs.get('chainIDmanager',None)
            psf=objs.get('psf','')

            missings=EmptyResidueList([])
            ssbonds=SSBondList([])
            links=LinkList([])
            ters=TerList([])
            seqadvs=SeqadvList([])
            if type(pr)==dict: # PDB format
                # minimal pr has ATOMS
                atoms=AtomList([Atom(p) for p in pr[Atom.PDB_keyword]])
                atoms.extend([Hetatm(p) for p in pr.get(Hetatm.PDB_keyword,[])])
                ters=TerList([Ter(p) for  p in pr.get(Ter.PDB_keyword,[])])
                atoms.adjustSerials(ters)
                modmanager.injest(ters)
                # mods.seqmods.terminals=ters
                if 'REMARK.465' in pr:
                    missings=EmptyResidueList([EmptyResidue(p) for p in pr['REMARK.465'].tables['MISSING']])
                ssbonds=SSBondList([SSBond(p) for p in pr.get(SSBond.PDB_keyword,[])])  
                seqadvs=SeqadvList([Seqadv(p) for p in pr.get(Seqadv.PDB_keyword,[])])
                links=LinkList([Link(p) for p in pr.get(Link.PDB_keyword,[])])
                links.remove_duplicates(fields=['chainID1','resseqnum1','insertion1','chainID2','resseqnum2','insertion2']) # some pdb files list links multiple times (2ins, i'm looking at you)
            elif type(pr)==DataContainer: # mmCIF format
                obj=pr.getObj(Atom.mmCIF_name)
                atoms=AtomList([Atom(CIFdict(obj,i)) for i in range(len(obj))])
                obj=pr.getObj(Seqadv.mmCIF_name)
                seqadvs=SeqadvList([Seqadv(CIFdict(obj,i)) for i in range(len(obj))])
                logger.debug(f'{len(seqadvs)} seqadv detected; typekeys: {[s.typekey for s in seqadvs]}')
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
                atoms.apply_psf_resnames(self.psf.atoms)
            seqmods=modmanager.get('seqmods',{})
            grafts=seqmods.get('grafts',GraftList([]))
            
            # Build the list of residues
            fromAtoms=ResidueList(atoms)
            fromEmptyResidues=ResidueList(missings)
            residues=fromAtoms+fromEmptyResidues
            # self.injest_grafts(grafts,residues,links)
            if sourcespecs.get('cif_residue_map_file',''):
                write_residue_map(residues.cif_residue_map(),sourcespecs['cif_residue_map_file'])
            residues.apply_segtypes()
            # apply seqmods
            if 'deletions' in seqmods:
                residues.deletion(seqmods['deletions'])
            if 'substitutions' in seqmods:
                new_seqadv,wearegone=residues.substitutions(seqmods['substitutions'])
                seqadvs.extend(new_seqadv)

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

            logger.debug(f'Segtypes present: {uniques["segtype"]}')
            # Delete any residues dictated by user-specified exclusions
            excludes=sourcespecs.get('exclude',{})
            thru_dict={resattr:excludes.get(yaml,[]) for yaml,resattr in self.excludables.items()}
            logger.debug(f'Exclusions: {thru_dict}')
            # delete residues that are in user-specified exclusions
            ignored_residues=residues.prune_exclusions(**thru_dict)
            # populates a specific mapping of PDB chainID to CIF label_asym_id
            # This is only meaningful if mmCIF input is used
            residues.map_chainIDs_label_to_auth()
            # initialize the chainID manager
            chainIDmanager.register_asymm_chains(residues.uniqattrs(['chainID'])['chainID'])
            logger.debug(f'Used chainIDs {chainIDmanager.Used}')

            # Give each Seqadv a residue identifier
            ignored_seqadvs=seqadvs.assign_residues(residues)
            ignored_ssbonds=ssbonds.assign_residues(residues)
            for b in ssbonds:
                logger.debug(f'after assign_residues {b}')
            more_ignored_residues,ignored_links=links.assign_residues(residues)
            for b in links:
                logger.debug(f'after assign_residues {b}')
            # a deleted link may create a "free" glycan; in this case
            # we should also delete its residues; problem is that 
            ignored_residues.extend(more_ignored_residues)
            ignored_grafts,more_ignored_residues,new_ignored_links=grafts.assign_residues(residues,links)
            ignored_residues.extend(more_ignored_residues)
            ignored_links.extend(new_ignored_links)
            
            if excludes or len(ignored_residues+ignored_seqadvs+ignored_ssbonds+ignored_links+ignored_grafts)>0:
                logger.debug(f'Exclusions result in deletion of:')
                logger.debug(f'    {len(ignored_residues)} residues; {len(residues)} remain;')
                logger.debug(f'    {len(ignored_seqadvs)} seqadvs; {len(seqadvs)} remain;')
                logger.debug(f'    {len(ignored_ssbonds)} ssbonds; {len(ssbonds)} remain;')
                logger.debug(f'    {len(ignored_links)} links; {len(links)} remain;')
                logger.debug(f'    {len(ignored_grafts)} grafts; {len(grafts)} remain.')
            
            # provide specifications of how to handle sequence issues
            # implied by PDB input
            seq_specs=sourcespecs.get('sequence',{})
            segments=SegmentList(seq_specs,residues,chainIDmanager)
            # this may have altered chainIDs for some residues.  So it is best
            # to be sure all mods that are residue-specific are updated
            seqadvs.update_attr_from_obj_attr('chainID','residue','chainID')
            ssbonds.update_attr_from_obj_attr('chainID1','residue1','chainID')
            ssbonds.update_attr_from_obj_attr('chainID2','residue2','chainID')
            links.update_attr_from_obj_attr('chainID1','residue1','chainID')
            links.update_attr_from_obj_attr('chainID2','residue2','chainID')
            grafts.update_attr_from_objlist_elem_attr('chainID','residues',0,'chainID')
            logger.debug(f'Segnames in A.U.: {",".join(segments.segnames)}')
            if segments.daughters:
                logger.debug(f'Daughter chains generated: {segments.daughters}')
            logger.debug(f'Used chainIDs {chainIDmanager.Used}')
            next_resseqnum=max([x.resseqnum for x in residues])+1
            for g in grafts:
                next_resseqnum=g.set_links(next_resseqnum)
                my_logger(str(g),logger.debug,fill=' ',just='<')
                links.extend(g.my_links)

            # at this point, we have built the asymmetric unit according to the intention of the 
            # author of the structure AND the intention of the user in excluding certain parts
            # of that structure (ligans, ions, chains, etc).  At this point, we should apply
            # any user-defined modifications

            # First, scan all seqadv's for relevant mutations to apply
            mutations=MutationList([])
            if seq_specs.get('fix_engineered_mutations',False):
                mutations.extend(MutationList([Mutation(s) for s in seqadvs if s.typekey=='engineered']))
            if seq_specs.get('fix_conflicts',False):
                mutations.extend(MutationList([Mutation(s) for s in seqadvs if s.typekey=='conflict']))
            mutations.extend(MutationList([Mutation(s) for s in seqadvs if s.typekey=='user']))
            # Now append these to the modmanager's mutations
            mutations=modmanager.injest(mutations)
            logger.debug(f'All mutations')
            for m in mutations:
                logger.debug(str(m))
            pruned_ssbonds=ssbonds.prune_mutations(mutations)
            pruned=pruned_by_links=links.prune_mutations(mutations,segments)
            pruned_segments=pruned.get('segments',[])
            for s in pruned_segments:
                logger.debug(f'Deactivating chainID {s.segname} because this entire segment was pruned')
                chainIDmanager.receive_chain(s.segname)

            # Now any added or deleted ssbonds
            ssbonds=modmanager.injest(ssbonds)
            ssbonds=modmanager['topomods']['ssbonds']
            if 'ssbondsdelete' in modmanager['topomods']:
                for s in ssbonds:
                    if modmanager['topomods']['ssbondsdelete'].is_deleted(s):
                        ignored_ssbonds.append(ssbonds.remove(s))

            ssbonds=modmanager.injest(ssbonds,overwrite=True)
            links=modmanager.injest(links)
            grafts=modmanager.injest(grafts,overwrite=True)
            
            segments.inherit_mods(modmanager)

            input_dict={
                'atoms':atoms,
                'residues':residues,
                'modmanager':modmanager,
                'segments':segments,
                'ignored':Namespace(residues=ignored_residues,links=ignored_links,ssbonds=ignored_ssbonds,seqadv=ignored_seqadvs),
                'pruned':Namespace(ssbonds=pruned_ssbonds,residues=pruned_by_links['residues'],links=pruned_by_links['links'],segments=pruned_by_links['segments'])
            }
        super().__init__(input_dict)

    def add_segment(self,seg):
        self.segments.append(seg)

    def injest_grafts(self,grafts,residues,links):
        # For each graft, add its residues beyond the graftpoint and the links involving those residues from the graft molecule as copies!
        # also need to reassign resids to these grafts and update any chainIDs
        for g in grafts:
            graft_chainID=g.chainID
            chain_residues=residues.filter(chainID=graft_chainID)
            # last_chain_residue_idx=residues.index(chain_residues[-1])
            next_available_resid=max([x.resseqnum for x in chain_residues])+1
            g_topomods=g.source_molecule.modmanager.get('topomods',{})
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
            # only injest links that are internal to this set of residues
            injested_links=0
            for l in g_links:
                if l.residue1 in g.my_residues and l.residue2 in g.my_residues:
                    links.append(l)
                    injested_links+=1
            logger.debug(f'Chain {graft_chainID} of raw asymmetric unit injests {len(g.my_residues)} residues and {injested_links} links from graft {g.id}')
                # residues.insert(last_chain_residue_idx+1,r)
                # last_chain_residue_idx+=1


    def set_coords(self,altstruct):
        atoms=AtomList([Atom(p) for p in altstruct[Atom.PDB_keyword]])
        atoms.extend([Hetatm(p) for p in altstruct.get(Hetatm.PDB_keyword,[])])
        ters=TerList([Ter(p) for  p in altstruct.get(Ter.PDB_keyword,[])])
        atoms.adjustSerials(ters)
        altRes=ResidueList(atoms)
        overwrites=[]
        for ar in altRes:
            tr=self.residues.get(resname=ar.resname,resseqnum=ar.resseqnum,insertion=ar.insertion,chainID=ar.chainID)
            if tr:
                tr.atoms.overwritePositions(ar.atoms)
                overwrites.append(ar)
        logger.debug(f'set_coords: {len(overwrites)} residues overwritten')

