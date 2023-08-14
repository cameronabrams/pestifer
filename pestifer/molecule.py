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
    req_attr=AncestorAwareMod.req_attr+['molid','source','asymmetric_unit','biological_assemblies','pdb_entry']
    opt_attr=AncestorAwareMod.opt_attr+['active_biological_assembly']
    _molcounter=0
    def __init__(self,**options):
        reset=options.get('reset_counter',False)
        if reset:
            Molecule._molcounter=0
        else:
            Molecule._molcounter+=1
        source=options.get('source',None)
        if not source:
            input_dict={
                'molid': Molecule._molcounter,
                'source': None,
                'pdb_entry': None,
                'asymmetric_unit': None,
                'biological_assemblies': None
            }
            super().__init__(input_dict)
        p_struct=PDBParser(PDBcode=source).parse().parsed
        input_dict={
            'molid': Molecule._molcounter,
            'source': source,
            'pdb_entry': p_struct,
            'asymmetric_unit': AsymmetricUnit(p_struct),
            'biological_assemblies': BioAssembList(p_struct,reset=True)
        }
        super().__init__(input_dict)
        self.asymmetric_unit.claim_descendants(self,0)
        self.biological_assemblies.claim_descendants(self,0)
    
    def activate_biological_assembly(self,index,chainIDmanager):
        biological_assembly=[x for x in self.biological_assemblies if x.index==index]
        assert biological_assembly!=[],f'No biological assembly "{index:d}" found.'
        assert len(biological_assembly)==1,f'No unique biological assembly "{index}" found.'
        self.active_biological_assembly=biological_assembly[0]
        logger.info(f'Activating biological assembly {self.active_biological_assembly.name} (idx {index})')
        ba=self.active_biological_assembly
        auChainIDs=self.asymmetric_unit.chainIDs
        ba.biomt[0].chainIDmap={c:c for c in auChainIDs}
        for biomt in ba.biomt[1:]:
            biomt.chainIDmap=chainIDmanager.generate_map(auChainIDs)            
        return self
    
    def write_TcL(self,B:ByteCollector,user_mods,file_collector=None):
        # TODO: merge the A.U.'s mods into the mods passed in from config
        au=self.asymmetric_unit
        ba=self.active_biological_assembly
        allmods=user_mods.copy()
        if not 'Mutations' in allmods:
            allmods['Mutations']=MutationList([])
        if not 'Conflicts' in allmods:
            allmods['Conflicts']=MutationList([])
        for k in ba.biomt[0].chainIDmap.keys():
            if au.Mutations:
                allmods['Mutations'].extend(au.Mutations.filter(chainID=k))
            if au.Conflicts:
                allmods['Conflicts'].extend(au.Conflicts.filter(chainID=k))
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
            au.Segments.write_TcL(B,biomt,allmods,file_collector=file_collector)
            B.banner('DISU PATCHES FOLLOW')
            B.write(au.SSBonds.write_TcL(biomt,allmods))
            B.banner('LINK PATCHES FOLLOW')
            B.write(au.Links.write_TcL(biomt,allmods))
            B.banner(f'TRANSFORM {biomt.index} ENDS')
        return 0

