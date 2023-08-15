import logging
from .mods import *
from .residue import ResidueList,AtomList,Atom,Hetatm,Ter,TerList
from .segment import SegmentList
from .basemod import AncestorAwareMod

logger=logging.getLogger(__name__)

class AsymmetricUnit(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['Atoms','Residues','SSBonds','Mutations','Conflicts','Missings','Links']
    opt_attr=AncestorAwareMod.opt_attr+['Ters','Segments','chainIDs']
    def __init__(self,pr:dict):
        # minimal pr has ATOMS
        Atoms=AtomList([Atom(p) for p in pr['ATOM']])
        if 'HETATM' in pr:
            Atoms.extend([Hetatm(p) for p in pr['HETATM']])
        if 'TER' in pr:
            Ters=TerList([Ter(p) for  p in pr['TER']])
            Atoms.adjustSerials(Ters) # VMD (1) ignors TERs and (2) *computes and assigns* serials
        else:
            Ters=TerList([])
        if 'REMARK.465' in pr:
            Missings=MissingList([Missing(p) for p in pr['REMARK.465'].tables['MISSING']])
        else:
            Missings=MissingList([])
        Residues=ResidueList(Atoms)+ResidueList(Missings)
        if 'SSBOND' in pr:
            SSBonds=SSBondList([SSBond(p) for p in pr['SSBOND']])
        else:
            SSBonds=SSBondList([])

        if 'SEQADV' in pr:
            unresolved_sa=[x for x in pr['SEQADV'] if not ('ENGINEERED' in x.conflict or 'CONFLICT' in x.conflict)]
            logger.info('The following SEQADV record types found in input are not handled:')
            for r in unresolved_sa:
                logger.info(f'{r.conflict}')
            sa=[x for x in pr['SEQADV'] if 'ENGINEERED' in x.conflict]
            Mutations=MutationList([Mutation(Seqadv(p)) for p in sa])
            co=[x for x in pr['SEQADV'] if 'CONFLICT' in x.conflict]
            Conflicts=MutationList([Mutation(Seqadv(p)) for p in co])
        else:
            Mutations=MutationList([])
            Conflicts=MutationList([])
        if 'LINK' in pr:
            Links=LinkList([Link(p) for p in pr['LINK']])
            Residues.update_links(Links,Atoms)
        else:
            Links=LinkList([])
        Segments=SegmentList(Residues)
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
        super().__init__(input_dict)