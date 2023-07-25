"""

.. module:: molecule
   :synopsis: Manages all molecules
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from .rcsb import PDBParser,MolData,MolMeta
class Molecule:
    _molcounter=0
    def __init__(self):
        self.molid=Molecule._molcounter
        Molecule._molcounter+=1
    @classmethod
    def from_rcsb(cls,**options):
        pdb_code=options.get('pdb_code',None)
        if not pdb_code:
            return None
        inst=cls()
        inst.pdb_code=pdb_code
        inst.pdb_parser=PDBParser(PDBcode=pdb_code)
        inst.pdb_parser.fetch()

        # inst.Meta,inst.Data=inst.pdb_parser.parse()

        return inst
