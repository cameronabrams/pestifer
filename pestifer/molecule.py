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
from .stringthings import my_logger,ByteCollector

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
        ba.biomt[0].chainIDmap={}
        for biomt in ba.biomt[1:]:
            biomt.chainIDmap=chainIDmanager.generate_map(auChainIDs)            
        return self
    
    def write_TcL(self,mods):
        B=ByteCollector()
        au=self.asymmetric_unit
        ba=self.active_biological_assembly
        for biomt in ba.biomt:
            B.comment(f'TRANSFORM {biomt.index} BEGINS')
            if biomt.chainIDmap:
                B.comment('The following mappings of A.U. chains is used:')
                for k,v in biomt.chainIDmap.items():
                    B.addline(f'#   {k}: {v}')
            B.write(au.Segments.write_TcL(biomt,mods))
            B.write(au.SSBonds.write_TcL(biomt,mods))
            B.write(au.Links.write_TcL(biomt,mods))
            B.comment(f'TRANSFORM {biomt.index} ENDS')
        return B.byte_collector

