import logging
# from pidibble.pdbparse import PDBParser
# from .basemod import BaseMod
from .mods import *
from .residue import Residue,ResidueList,LinkList
from .segment import Segment,SegmentList
from .config import ConfigGetParam
from .basemod import CloneableMod,AncestorAwareMod

logger=logging.getLogger(__name__)

class AsymmetricUnit(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['Atoms','Residues','SSBonds','Mutations','Conflicts','Missings','Links']
    opt_attr=AncestorAwareMod.opt_attr+['Segments','chainIDs']
    @classmethod
    def from_pdb(cls,pr:dict):
        Atoms=AtomList([Atom.from_pdbrecord(p) for p in pr['ATOM']])
        if 'HETATM' in pr:
            Atoms.extend([Hetatm.from_pdbrecord(p) for p in pr['HETATM']])
        if 'REMARK.465' in pr:
            Missings=MissingList([Missing.from_pdbrecord(p) for p in pr['REMARK.465'].tables['MISSING']])
        else:
            Missings=MissingList([])
        Residues=ResidueList.from_atoms(Atoms)+ResidueList.from_missing(Missings)
        # Residues.sort()
        # print(f'Residues {len(Residues)}')
        if 'SSBOND' in pr:
            SSBonds=SSBondList([SSBond.from_pdbrecord(p) for p in pr['SSBOND']])
        else:
            SSBonds=SSBondList([])

        if 'SEQADV' in pr:
            ''' These are possible conflict keywords
            - Cloning artifact
            - Expression tag
            - Conflict
            - Engineered
            - Variant 
            - Insertion
            - Deletion
            - Microheterogeneity
            - Chromophore
            '''
            unresolved_sa=[x for x in pr['SEQADV'] if not ('ENGINEERED' in x.conflict or 'CONFLICT' in x.conflict)]
            logger.warning('The following SEQADV records are not handled:')
            for r in unresolved_sa:
                logger.warning(f'{r.conflict}')
            sa=[x for x in pr['SEQADV'] if 'ENGINEERED' in x.conflict]
            Mutations=MutationList([Mutation.from_seqadv(Seqadv.from_pdbrecord(p)) for p in sa])
            co=[x for x in pr['SEQADV'] if 'CONFLICT' in x.conflict]
            Conflicts=MutationList([Mutation.from_seqadv(Seqadv.from_pdbrecord(p)) for p in co])
        else:
            Mutations=MutationList([])
            Conflicts=MutationList([])
        if 'LINK' in pr:
            Links=LinkList([Link.from_pdbrecord(p) for p in pr['LINK']]).update_residues_atoms(Residues,Atoms)
        else:
            Links=LinkList([])
        Segments=SegmentList.from_residues(Residues)
        # TODO: Maybe have segments "claim" associated Mutations and Conflicts
        chainIDs=list(set([x.chainID for x in Residues]))
        chainIDs.sort()
        input_dict={
            'Atoms':Atoms,
            'Residues':Residues,
            'SSBonds':SSBonds,
            'Mutations':Mutations,
            'Conflicts':Conflicts,
            'Missings':Missings,
            'Links':Links,
            'Segments':Segments,
            'chainIDs':chainIDs
        }
        return cls(input_dict)