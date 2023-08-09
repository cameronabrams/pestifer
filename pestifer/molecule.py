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

logger=logging.getLogger(__name__)

class Molecule(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['molid','pdb_code','asymmetric_unit','biological_assemblies','pdb_entry']
    # opt_attr=[]
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
