# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A class for building the asymmetric unit from a PDB file

"""
import logging
from .mods import *
from argparse import Namespace
from .residue import ResidueList,AtomList,Atom,Hetatm,Ter,TerList,EmptyResidue,EmptyResidueList
from .segment import SegmentList
from .basemod import AncestorAwareMod
from .modcontainer import ModContainer
from mmcif.api.PdbxContainers import DataContainer
from .config import *
from .util import write_residue_map

logger=logging.getLogger(__name__)

class AsymmetricUnit(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['atoms','residues','mods']
    opt_attr=AncestorAwareMod.opt_attr+['segments','ignored','pruned']
    excludables={'resnames':'resname','chains':'chainID','segtypes':'segtype'}
    def __init__(self,**objs):
        if len(objs)==0:
            logger.debug('Generating an empty A.U.')
            input_dict={
                'atoms':AtomList([]),
                'residues':ResidueList([]),
                'mods':ModContainer(),
                'segments':SegmentList({},ResidueList([]))
            }
        else: #
            pr=objs.get('parsed',None)
            sourcespecs=objs.get('sourcespecs',{})
            mods=objs.get('usermodspecs',ModContainer())
            # logger.debug(f'User mods {type(mods)} at asymmetric unit: {mods.__dict__}')
            chainIDmanager=objs.get('chainIDmanager',None)
            
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
                mods.seqmods.terminals=ters
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

            # Build the list of residues
            fromAtoms=ResidueList(atoms)
            fromEmptyResidues=ResidueList(missings)
            residues=fromAtoms+fromEmptyResidues
            if sourcespecs.get('cif_residue_map_file',''):
                write_residue_map(residues.cif_residue_map(),sourcespecs['cif_residue_map_file'])
            residues.apply_segtypes()
            # apply seqmods
            if hasattr(mods.seqmods,'deletions'):
                residues.deletion(mods.seqmods.deletions)
            if hasattr(mods.seqmods,'substitutions'):
                new_seqadv,wearegone=residues.substitutions(mods.seqmods.substitutions)
                seqadvs.extend(new_seqadv)

            uniques=residues.uniqattrs(['segtype'],with_counts=True)
            nResolved=sum([1 for x in residues if x.resolved])
            nUnresolved=sum([1 for x in residues if not x.resolved])
            logger.debug(f'{len(residues)} total residues: {nResolved} resolved and {nUnresolved} unresolved')
            if hasattr(mods.seqmods,'deletions'):
                nResolvedDelete=len(fromAtoms)-nResolved
                nUnresolvedDelete=len(fromEmptyResidues)-nUnresolved
                logger.debug(f'Deletions removed {nResolvedDelete} resolved and {nUnresolvedDelete} unresolved residues')
            logger.debug(f'Segtypes present: {uniques["segtype"]}')
            # Delete any residues dictated by user-specified exclusions
            excludes=sourcespecs['exclude']
            thru_dict={resattr:excludes.get(yaml,[]) for yaml,resattr in self.excludables.items()}
            logger.debug(f'Exclusions: {thru_dict}')
            # delete residues that are in user-specified exclusions
            ignored_residues=residues.prune_exclusions(**thru_dict)
            # populates a specific mapping of PDB chainID to CIF label_asym_id
            # This is only meaningful if mmCIF input is used
            residues.map_chainIDs_label_to_auth()
            # initialize the chainID manager
            chainIDmanager.register_asymm_chains(residues.uniqattrs(['chainID'])['chainID'])
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
            
            if excludes:
                logger.debug(f'Exclusions result in deletion of:')
                logger.debug(f'    {len(ignored_residues)} residues; {len(residues)} remain;')
                logger.debug(f'    {len(ignored_seqadvs)} seqadvs; {len(seqadvs)} remain;')
                logger.debug(f'    {len(ignored_ssbonds)} ssbonds; {len(ssbonds)} remain; and')
                logger.debug(f'    {len(ignored_links)} links; {len(links)} remain.')

            seq_specs=sourcespecs['sequence']
            segments=SegmentList(seq_specs,residues,chainIDmanager)
            # this may have altered chainIDs for some residues.  So it is best
            # to be sure all mods that are residue-specific are updated
            seqadvs.update_attr_from_obj_attr('chainID','residue','chainID')
            ssbonds.update_attr_from_obj_attr('chainID1','residue1','chainID')
            ssbonds.update_attr_from_obj_attr('chainID2','residue2','chainID')
            links.update_attr_from_obj_attr('chainID1','residue1','chainID')
            links.update_attr_from_obj_attr('chainID2','residue2','chainID')
            logger.debug(f'Segnames in A.U.: {",".join(segments.segnames)}')
            if segments.daughters:
                logger.debug(f'Daughter chains generated: {segments.daughters}')

            # at this point, we have built the asymmetric unit according to the intention of the 
            # author of the structure AND the intention of the user in excluding certain parts
            # of that structure (ligans, ions, chains, etc).  At this point, we should apply
            # any user-defined modifications

            # First, any request to revert to database sequencing:
            mutations=MutationList([])
            if sourcespecs['sequence']['fix_engineered_mutations']:
                mutations.extend(MutationList([Mutation(s) for s in seqadvs if s.typekey=='engineered']))
            if sourcespecs['sequence']['fix_conflicts']:
                mutations.extend(MutationList([Mutation(s) for s in seqadvs if s.typekey=='conflict']))
            mutations.extend(MutationList([Mutation(s) for s in seqadvs if s.typekey=='user']))
            if hasattr(mods.seqmods,'mutations'):
                mutations.extend(mods.seqmods.mutations)
            logger.debug(f'mutations')
            for m in mutations:
                logger.debug(str(m))
            mods.seqmods.mutations=mutations
            pruned_ssbonds=ssbonds.prune_mutations(mutations)
            pruned=pruned_by_links=links.prune_mutations(mutations,segments)
            if hasattr(pruned,'segments'):
                for s in pruned['segments']:
                    logger.debug(f'Recovering chainID {s.segname}')
                    chainIDmanager.receive_chain(s.segname)
            # Now any added or deleted ssbonds
            if hasattr(mods.topomods,'ssbonds'):
                ssbonds.extend(mods.topomods.ssbonds)
            if hasattr(mods.topomods,'ssbondsdelete'):
                for s in ssbonds:
                    if mods.topomods.ssbondsdelete.is_deleted(s):
                        ignored_ssbonds.append(ssbonds.remove(s))
            mods.topomods.ssbonds=ssbonds
            # Now any added links
            if hasattr(mods.topomods,'links'):
                links.extend(mods.topomods.links)
            mods.topomods.links=links
            
            input_dict={
                'atoms':atoms,
                'residues':residues,
                'mods':mods,
                'segments':segments,
                'ignored':Namespace(residues=ignored_residues,links=ignored_links,ssbonds=ignored_ssbonds,seqadv=ignored_seqadvs),
                'pruned':Namespace(ssbonds=pruned_ssbonds,residues=pruned_by_links['residues'],links=pruned_by_links['links'],segments=pruned_by_links['segments'])
            }
        super().__init__(input_dict)