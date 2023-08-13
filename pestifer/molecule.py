"""

.. module:: molecule
   :synopsis: Manages all molecules
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
from pidibble.pdbparse import PDBParser
from .basemod import AncestorAwareMod
from .asymmetricunit import AsymmetricUnit
from .bioassemb import BioAssembList
from .stringthings import ByteCollector
from .mods import MutationList
from .config import ConfigGetParam
logger=logging.getLogger(__name__)

class Molecule(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['molid','pdb_code','asymmetric_unit','biological_assemblies','pdb_entry']
    opt_attr=AncestorAwareMod.opt_attr+['active_biological_assembly']
    _molcounter=0
    def __init__(self,input_dict):
        Molecule._molcounter+=1
        super().__init__(input_dict)

    @classmethod
    def from_rcsb(cls,**options):
        pdb_code=options.get('pdb_code',None)
        if not pdb_code:
            return None
        p_struct=PDBParser(PDBcode=pdb_code).parse().parsed
        input_dict={
            'molid': Molecule._molcounter,
            'pdb_code': pdb_code,
            'pdb_entry': p_struct,
            'asymmetric_unit': AsymmetricUnit.from_pdb(p_struct),
            'biological_assemblies': BioAssembList.from_pdb(p_struct)
        }
        inst=cls(input_dict)
        inst.asymmetric_unit.claim_descendants(inst,0)
        inst.biological_assemblies.claim_descendants(inst,0)
        return inst
    
    def activate_biological_assembly(self,index,chainIDmanager):
        biological_assembly=[x for x in self.biological_assemblies if x.index==index]
        assert biological_assembly!=[],f'No biological assembly "{index}" found.'
        assert len(biological_assembly)==1,f'No unique biological assembly "{index}" found.'
        self.active_biological_assembly=biological_assembly[0]
        logger.info(f'Activating biological assembly {self.active_biological_assembly.name} (idx {index})')
        ba=self.active_biological_assembly
        auChainIDs=self.asymmetric_unit.chainIDs
        ba.biomt[0].chainIDmap={c:c for c in auChainIDs}
        for biomt in ba.biomt[1:]:
            biomt.chainIDmap=chainIDmanager.generate_map(auChainIDs)            
        return self
    
    def write_TcL(self,user_mods):
        # TODO: merge the A.U.'s mods into the mods passed in from config
        B=ByteCollector()
        au=self.asymmetric_unit
        ba=self.active_biological_assembly
        allmods=user_mods.copy()
        if not 'Mutations' in allmods:
            allmods['Mutations']=MutationList([])
        if not 'Conflicts' in allmods:
            allmods['Conflicts']=MutationList([])
        for k in ba.biomt[0].chainIDmap.keys():
            if au.Mutations:
                allmods['Mutations'].extend(au.Mutations.get(chainID=k))
            if au.Conflicts:
                allmods['Conflicts'].extend(au.Conflicts.get(chainID=k))
        # TODO: Delete any SSBonds or Links that contain residues that are mutated
        au.SSBonds.prune_mutations(user_mods.get('Mutations',MutationList([])))
        if ConfigGetParam('Fix_engineered_mutations'):
            au.SSBonds.prune_mutations(au.Mutations)
        if ConfigGetParam('Fix_conflicts'):
            au.SSBonds.prune_mutations(au.Conflicts)
        for biomt in ba.biomt:
            B.banner(f'TRANSFORM {biomt.index} BEGINS')
            B.banner('The following mappings of A.U. chains is used:')
            for k,v in biomt.chainIDmap.items():
                B.comment(f'{k}: {v}')
            B.banner('SEGMENTS FOLLOW')
            B.write(au.Segments.write_TcL(biomt,allmods))
            B.banner('DISU PATCHES FOLLOW')
            B.write(au.SSBonds.write_TcL(biomt,allmods))
            B.banner('LINK PATCHES FOLLOW')
            B.write(au.Links.write_TcL(biomt,allmods))
            B.banner(f'TRANSFORM {biomt.index} ENDS')
        return B.byte_collector

