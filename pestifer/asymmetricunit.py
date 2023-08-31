"""

.. module:: asymmetricunit
   :synopsis: defines the Asymmetricunit class
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
from .mods import *
from .residue import ResidueList,AtomList,Atom,Hetatm,Ter,TerList
from .segment import SegmentList
from .basemod import AncestorAwareMod
from mmcif.api.PdbxContainers import DataContainer

logger=logging.getLogger(__name__)

class AsymmetricUnit(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['Atoms','Residues','SSBonds','Mutations','Conflicts','Missings','Links']
    opt_attr=AncestorAwareMod.opt_attr+['Ters','Segments','chainIDs']
    # yaml-label: residue-attr-name
    excludables={'resnames':'name','chains':'chainID','segtypes':'segtype'}
    def __init__(self,*objs):
        if len(objs)==0:
            logger.debug('Generating an empty A.U.')
            input_dict={
                'Atoms':AtomList([]),
                'Residues':ResidueList([]),
                'SSBonds':SSBondList([]),
                'Mutations':MutationList([]),
                'Conflicts':MutationList([]),
                'Missings':MissingList([]),
                'Links':LinkList([]),
                'Ters':TerList([]),
                'Segments':SegmentList({},ResidueList([])),
                'chainIDs':[]
            }
        else: # p_struct,config,chainIDmanager,excludes
            pr=objs[0]
            config=objs[1]
            chainIDmanager=objs[2]
            excludes=objs[3]
            Missings=MissingList([])
            SSBonds=SSBondList([])
            Mutations=MutationList([])
            Conflicts=MutationList([])
            Seqadvs=SeqadvList([])
            Links=LinkList([])
            Ters=TerList([]) # no TER records in an mmCIF file!
            # Build the lists of
            # - Atoms
            # - Seqadvs (engineered mutations, conflicts, etc.)
            # - Missing residues
            # - Disulfides
            # - Links
            if type(pr)==dict: # PDB format
                assert config['rcsb_file_format']=='PDB'
                # minimal pr has ATOMS
                Atoms=AtomList([Atom(p) for p in pr['ATOM']])
                if 'HETATM' in pr:
                    Atoms.extend([Hetatm(p) for p in pr['HETATM']])
                if 'TER' in pr:
                    Ters=TerList([Ter(p) for  p in pr['TER']])
                    Atoms.adjustSerials(Ters) # VMD (1) ignors TERs and (2) *computes and assigns* serials
                if 'REMARK.465' in pr:
                    Missings=MissingList([Missing(p) for p in pr['REMARK.465'].tables['MISSING']])
                if 'SSBOND' in pr:
                    SSBonds=SSBondList([SSBond(p) for p in pr['SSBOND']])  
                if 'SEQADV' in pr:
                    Seqadvs=SeqadvList([Seqadv(p) for p in pr['SEQADV']])
                if 'LINK' in pr:
                    Links=LinkList([Link(p) for p in pr['LINK']])
            elif type(pr)==DataContainer: # mmCIF format
                assert config['rcsb_file_format']=='mmCIF'
                obj=pr.getObj('atom_site')
                Atoms=AtomList([Atom(CIFdict(obj,i)) for i in range(len(obj))])
                obj=pr.getObj('struct_ref_seq_dif')
                Seqadvs=SeqadvList([Seqadv(CIFdict(obj,i)) for i in range(len(obj))])
                obj=pr.getObj('pdbx_unobs_or_zero_occ_residues')
                Missings=MissingList([Missing(CIFdict(obj,i)) for i in range(len(obj))])
                obj=pr.getObj('struct_conn')
                ssdicts=[]
                lndicts=[]
                for i in range(len(obj)):
                    conn_type_id=obj.getValue('conn_type_id',i)
                    if conn_type_id=='disulf':
                        ssdicts.append(CIFdict(obj,i))
                    elif conn_type_id=='covale':
                        lndicts.append(CIFdict(obj,i))
                SSBonds=SSBondList([SSBond(x) for x in ssdicts])
                Links=LinkList([Link(x) for x in lndicts])

            # Build the list of residues
            fromAtoms=ResidueList(Atoms)
            fromMissings=ResidueList(Missings)
            Residues=fromAtoms+fromMissings
            Residues.apply_segtypes(config.get('Segtypes_by_Resnames',{}))
            uniques=Residues.uniqattrs(['segtype'],with_counts=True)
            logger.debug(f'{len(Residues)} total residues: {len(fromAtoms)} resolved and {len(fromMissings)} unresolved')
            logger.debug(f'Segtypes present: {uniques["segtype"]}')
            # Delete any residues dictated by user-specified exclusions
            thru_dict={resattr:excludes.get(yaml,[]) for yaml,resattr in self.excludables.items()}
            logger.debug(f'Exclusions: {thru_dict}')
            # delete residues that are in user-specified exclusions
            ignored_residues=Residues.prune_exclusions(**thru_dict)
            # populates a specific mapping of PDB chainID to CIF label_asym_id
            # This is only meaningful if mmCIF input is used
            Residues.map_chainIDs_label_to_auth()
            # initialize the chainID manager
            chainIDmanager.register_asymm_chains(Residues.uniqattrs(['chainID'])['chainID'])
            # Give each Seqadv a residue identifier
            ignored_seqadvs=Seqadvs.assign_residues(Residues)
            ignored_ssbonds=SSBonds.assign_residues(Residues)
            more_ignored_residues,ignored_links=Links.assign_residues(Residues)
            # a deleted link may create a "free" glycan; in this case
            # we should also delete its residues; problem is that 
            ignored_residues.extend(more_ignored_residues)
            
            if excludes:
                logger.debug(f'Exclusions result in deletion of:')
                logger.debug(f'    {len(ignored_residues)} residues; {len(Residues)} remain;')
                logger.debug(f'    {len(ignored_seqadvs)} seqadvs; {len(Seqadvs)} remain;')
                logger.debug(f'    {len(ignored_ssbonds)} ssbonds; {len(SSBonds)} remain; and')
                logger.debug(f'    {len(ignored_links)} links; {len(Links)} remain.')

            Segments=SegmentList(config,Residues,chainIDmanager)
            # this may have altered chainIDs for some residues.  So it is best
            # to be sure all mods that are residue-specific are updated
            Seqadvs.update_attr_from_obj_attr('chainID','residue','chainID')
            # we will only conside the impact of mutations on links/ssbonds 
            # in the context of generating the segment output, where 
            # - user mods are included, and
            # - user directives for acting on mutations/conflicts in the input header appear
            Mutations=MutationList([Mutation(s) for s in Seqadvs if 'ENGINEERED' in s.conflict])
            Conflicts=MutationList([Mutation(s) for s in Seqadvs if 'CONFLICT' in s.conflict])

            SSBonds.update_attr_from_obj_attr('chainID1','residue1','chainID')
            SSBonds.update_attr_from_obj_attr('chainID2','residue2','chainID')
            Links.update_attr_from_obj_attr('chainID1','residue1','chainID')
            Links.update_attr_from_obj_attr('chainID2','residue2','chainID')
            logger.debug(f'Segnames in A.U.: {",".join(Segments.segnames)}')
            if Segments.daughters:
                logger.debug(f'Daughter chains generated: {Segments.daughters}')
            input_dict={
                'Atoms':Atoms,
                'Residues':Residues,
                'Ters':Ters,
                'SSBonds':SSBonds,
                'Mutations':Mutations,
                'Conflicts':Conflicts,
                'Missings':Missings,
                'Links':Links,
                'Segments':Segments,
            }
        super().__init__(input_dict)