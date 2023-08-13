import logging
from .mods import *
from .residue import ResidueList,AtomList,Atom,Hetatm,Ter,TerList
from .segment import SegmentList
from .basemod import AncestorAwareMod

logger=logging.getLogger(__name__)

class AsymmetricUnit(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['Atoms','Residues','SSBonds','Mutations','Conflicts','Missings','Links']
    opt_attr=AncestorAwareMod.opt_attr+['Ters','Segments','chainIDs']
    @classmethod
    def from_pdb(cls,pr:dict):
        Atoms=AtomList([Atom.from_pdbrecord(p) for p in pr['ATOM']])
        if 'HETATM' in pr:
            Atoms.extend([Hetatm.from_pdbrecord(p) for p in pr['HETATM']])
        if 'TER' in pr:
            Ters=TerList([Ter.from_pdbrecord(p) for  p in pr['TER']])
            Atoms.adjustSerials(Ters) # VMD (1) ignors TERs and (2) *computes and assigns* serials
        else:
            Ters=TerList([])
        if 'REMARK.465' in pr:
            Missings=MissingList([Missing.from_pdbrecord(p) for p in pr['REMARK.465'].tables['MISSING']])
        else:
            Missings=MissingList([])
        Residues=ResidueList.from_atoms(Atoms)+ResidueList.from_missing(Missings)
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
            logger.info('The following SEQADV record types found in input are not handled:')
            for r in unresolved_sa:
                logger.info(f'{r.conflict}')
            sa=[x for x in pr['SEQADV'] if 'ENGINEERED' in x.conflict]
            Mutations=MutationList([Mutation.from_seqadv(Seqadv.from_pdbrecord(p)) for p in sa])
            co=[x for x in pr['SEQADV'] if 'CONFLICT' in x.conflict]
            Conflicts=MutationList([Mutation.from_seqadv(Seqadv.from_pdbrecord(p)) for p in co])
        else:
            Mutations=MutationList([])
            Conflicts=MutationList([])
        if 'LINK' in pr:
            Links=LinkList([Link.from_pdbrecord(p) for p in pr['LINK']])
            Residues.update_links(Links,Atoms)
        else:
            Links=LinkList([])
        Segments=SegmentList.from_residues(Residues)
        chainIDs=list(set([x.chainID for x in Residues]))
        chainIDs.sort()
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
            'chainIDs':chainIDs
        }
        return cls(input_dict)