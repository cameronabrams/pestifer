"""

.. module:: asymmetricunit
   :synopsis: defines the Asymmetricunit class
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
from .mods import *
from argparse import Namespace
from .residue import ResidueList,AtomList,Atom,Hetatm,Ter,TerList
from .segment import SegmentList
from .basemod import AncestorAwareMod
from mmcif.api.PdbxContainers import DataContainer
from .config import *

logger=logging.getLogger(__name__)

class AsymmetricUnit(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['atoms','residues','seqmods','topomods']
    opt_attr=AncestorAwareMod.opt_attr+['segments','ignored']
    excludables={'resnames':'name','chains':'chainID','segtypes':'segtype'}
    def __init__(self,**objs):
        if len(objs)==0:
            logger.debug('Generating an empty A.U.')
            input_dict={
                'atoms':AtomList([]),
                'residues':ResidueList([]),
                'seqmods':Namespace(
                    engr_mutations=MutationList([]),
                    conflicts=MutationList([]),
                    ters=TerList([])
                    ),
                'topomods':Namespace(
                    SSBonds=SSBondList([]),
                    Links=LinkList([])
                    ),
                'segments':SegmentList({},ResidueList([]))
            }
        else: #
            pr=objs.get('parsed',None)
            sourcespecs=objs.get('sourcespecs',{})
            chainIDmanager=objs.get('chainIDmanager',None)
            missings=MissingList([])
            topomods=Namespace(ssbonds=SSBondList([]),links=LinkList([]))
            seqmods=Namespace(mutations=MutationList([]),conflicts=MutationList([]),ters=TerList([]))
            seqadvs=SeqadvList([])
            if type(pr)==dict: # PDB format
                # minimal pr has ATOMS
                atoms=AtomList([Atom(p) for p in pr['ATOM']])
                if 'HETATM' in pr:
                    atoms.extend([Hetatm(p) for p in pr['HETATM']])
                if 'TER' in pr:
                    seqmods.ters=TerList([Ter(p) for  p in pr['TER']])
                    atoms.adjustSerials(seqmods.ters) # VMD (1) ignors TERs and (2) *computes and assigns* serials
                if 'REMARK.465' in pr:
                    missings=MissingList([Missing(p) for p in pr['REMARK.465'].tables['MISSING']])
                if 'SSBOND' in pr:
                    topomods.ssbonds=SSBondList([SSBond(p) for p in pr['SSBOND']])  
                if 'SEQADV' in pr:
                    seqadvs=SeqadvList([Seqadv(p) for p in pr['SEQADV']])
                if 'LINK' in pr:
                    topomods.links=LinkList([Link(p) for p in pr['LINK']])
            elif type(pr)==DataContainer: # mmCIF format
                obj=pr.getObj('atom_site')
                atoms=AtomList([Atom(CIFdict(obj,i)) for i in range(len(obj))])
                obj=pr.getObj('struct_ref_seq_dif')
                seqadvs=SeqadvList([Seqadv(CIFdict(obj,i)) for i in range(len(obj))])
                obj=pr.getObj('pdbx_unobs_or_zero_occ_residues')
                missings=MissingList([Missing(CIFdict(obj,i)) for i in range(len(obj))])
                obj=pr.getObj('struct_conn')
                ssdicts=[]
                lndicts=[]
                for i in range(len(obj)):
                    conn_type_id=obj.getValue('conn_type_id',i)
                    if conn_type_id=='disulf':
                        ssdicts.append(CIFdict(obj,i))
                    elif conn_type_id=='covale':
                        lndicts.append(CIFdict(obj,i))
                topomods.ssbonds=SSBondList([SSBond(x) for x in ssdicts])
                topomods.links=LinkList([Link(x) for x in lndicts])

            # Build the list of residues
            fromAtoms=ResidueList(atoms)
            fromMissings=ResidueList(missings)
            residues=fromAtoms+fromMissings
            residues.apply_segtypes()
            uniques=residues.uniqattrs(['segtype'],with_counts=True)
            logger.debug(f'{len(residues)} total residues: {len(fromAtoms)} resolved and {len(fromMissings)} unresolved')
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
            ignored_ssbonds=topomods.ssbonds.assign_residues(residues)
            more_ignored_residues,ignored_links=topomods.links.assign_residues(residues)
            # a deleted link may create a "free" glycan; in this case
            # we should also delete its residues; problem is that 
            ignored_residues.extend(more_ignored_residues)
            
            if excludes:
                logger.debug(f'Exclusions result in deletion of:')
                logger.debug(f'    {len(ignored_residues)} residues; {len(residues)} remain;')
                logger.debug(f'    {len(ignored_seqadvs)} seqadvs; {len(seqadvs)} remain;')
                logger.debug(f'    {len(ignored_ssbonds)} ssbonds; {len(topomods.ssbonds)} remain; and')
                logger.debug(f'    {len(ignored_links)} links; {len(topomods.links)} remain.')

            seq_specs=sourcespecs['sequence']
            segments=SegmentList(seq_specs,residues,chainIDmanager)
            # this may have altered chainIDs for some residues.  So it is best
            # to be sure all mods that are residue-specific are updated
            seqadvs.update_attr_from_obj_attr('chainID','residue','chainID')
            seqmods.engr_mutations=MutationList([Mutation(s) for s in seqadvs if s.typekey=='engineered'])
            seqmods.conflicts=MutationList([Mutation(s) for s in seqadvs if s.typekey=='conflict'])
            topomods.ssbonds.update_attr_from_obj_attr('chainID1','residue1','chainID')
            topomods.ssbonds.update_attr_from_obj_attr('chainID2','residue2','chainID')
            topomods.links.update_attr_from_obj_attr('chainID1','residue1','chainID')
            topomods.links.update_attr_from_obj_attr('chainID2','residue2','chainID')
            logger.debug(f'Segnames in A.U.: {",".join(segments.segnames)}')
            if segments.daughters:
                logger.debug(f'Daughter chains generated: {segments.daughters}')
            input_dict={
                'atoms':atoms,
                'residues':residues,
                'seqmods':seqmods,
                'topomods':topomods,
                'segments':segments,
                'ignored':Namespace(residues=ignored_residues,links=ignored_links,ssbonds=ignored_ssbonds,seqadv=ignored_seqadvs)
            }
        super().__init__(input_dict)