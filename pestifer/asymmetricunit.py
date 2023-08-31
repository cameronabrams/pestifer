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
            logger.debug(f'excludes: {excludes}')
            Missings=MissingList([])
            SSBonds=SSBondList([])
            Mutations=MutationList([])
            Conflicts=MutationList([])
            Seqadvs=SeqadvList([])
            Links=LinkList([])
            Ters=TerList([]) # no TER records in an mmCIF file!
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
            # if len(unresolved_sa)>0:
            #     logger.info('The following SEQADV record types found in input are not handled:')
            # for r in unresolved_sa:
            #     logger.info(r)
            fromAtoms=ResidueList(Atoms)
            fromMissings=ResidueList(Missings)
            # for r in fromMissings:
            #     logger.debug(f'missing {r.name}{r.resseqnum}{r.insertion} in {r.chainID} (auth {r.auth_comp_id}{r.auth_seq_id} in {r.auth_asym_id})')
            Residues=fromAtoms+fromMissings
            logger.debug(f'{len(Residues)} total residues: {len(fromAtoms)} resolved and {len(fromMissings)} unresolved')
            # populates a specific mapping of PDB chainID to CIF label_asym_id
            # This is only meaningful if mmCIF input is used
            Residues.map_chainIDs_label_to_auth()
            # initialize the chainID manager
            chainIDmanager.register_asymm_chains(Residues.unique_chainIDs())
            # Update all chainID attributes of Seqadvs; this is because
            # CIF files do not include the asym_id, only the author's chainID
            Seqadvs.setChainIDs(Residues)
            
            Mutations=MutationList([Mutation(s) for s in Seqadvs if 'ENGINEERED' in s.conflict])
            Conflicts=MutationList([Mutation(s) for s in Seqadvs if 'CONFLICT' in s.conflict])
            thru_dict={'name':excludes.get('resnames',[]),'chainID':excludes.get('chains',[])}
            logger.debug(f'Exclusions: {thru_dict}')
            # delete residues that are in user-specified exclusions
            ignored_residues=Residues.prune_exclusions(**thru_dict)
            if config!=None:
                Residues.apply_segtypes(config['Segtypes_by_Resnames'])
            # delete any disulfides or other bonds in which these deleted
            # residues appear
            attr_maps=[{'chainID1':'chainID','resseqnum1':'resseqnum','insertion1':'insertion'},{'chainID2':'chainID','resseqnum2':'resseqnum','insertion2':'insertion'}]
            SSBonds.prune(objlist=ignored_residues,attr_maps=attr_maps)
            Links.prune(objlist=ignored_residues,attr_maps=attr_maps)
            if config!=None:
                Links.apply_segtypes(config['Segtypes_by_Resnames'])
            Residues.update_links(Links,Atoms)

            Segments=SegmentList(config,Residues,chainIDmanager)
            # this may have altered segnames for some residues.  So it is best
            # to be sure all mods that are residue-specific are updated

            # chainIDs=Segments.segnames
            logger.debug(f'Segnames in A.U.: {",".join(Segments.segnames)}')
            # chainIDs.sort()
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
                # 'chainIDs':chainIDs
            }
        super().__init__(input_dict)