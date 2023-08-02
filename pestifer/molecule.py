"""

.. module:: molecule
   :synopsis: Manages all molecules
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from pidibble.pdbparse import PDBParser
from .mods import *
from .chain import Chain
from .residue import Residue,ResidueList
from .config import ConfigGetParam

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
        inst.pdb_entry=PDBParser(PDBcode=pdb_code).parse()
        inst.build(**options)
        return inst
    
    def build(self,**options):
        """build builds complete molecule description from parsed pdb_entry
        """
        pr=self.pdb_entry.parsed
        atprs=pr['ATOM']
        self.Atoms=AtomList([Atom.from_pdbrecord(p) for p in atprs])
        het=pr['HETATM']
        self.Atoms.extend([Hetatm.from_pdbrecord(p) for p in het])
        ssb=pr['SSBOND']
        self.SSBonds=SSBondList([SSBond.from_pdbrecord(p) for p in ssb])
        sa=[x for x in pr['SEQADV'] if 'MUTATION' in x.conflict]
        self.Mutations=MutationList([Mutation.from_seqadv(Seqadv.from_pdbrecord(p)) for p in sa])
        co=[x for x in pr['SEQADV'] if 'CONFLICT' in x.conflict]
        self.Conflicts=MutationList([Mutation.from_seqadv(Seqadv.from_pdbrecord(p)) for p in co])
        mi=[x for x in pr['REMARK.465'].tables['MISSING']]
        self.Missing=MissingList([Missing.from_pdbrecord(p) for p in mi])
        li=[x for x in pr['LINK']]
        self.Links=LinkList([Link.from_pdbrecord(p) for p in li])
        self._MakeResidues()
        self._MakeChains()

    def _MakeResidues(self):
        segtype_by_resname=ConfigGetParam('Segtypes_by_Resnames')
        ''' make residues from atoms '''
        R=[]
        a=self.Atoms[0]
        r=Residue.from_atom(a)
        R.append(r)
        for a in self.Atoms[1:]:
            if not any([x.add_atom(a) for x in R]):
                r=Residue.from_atom(a)
                R.append(r)
        ''' insert missing residues '''
        for m in self.Missing:
            R.append(Residue.from_missing(m))
        for r in R:
            r.segtype=segtype_by_resname.get(r.name,'')
        self.Residues=ResidueList(R)

    def _MakeChains(self):
        r=self.Residues[0]
        c=Chain.from_residue(r,parent_molecule=self)
        self.Chains=[]
        self.Chains.append(c)
        for r in self.Residues[1:]:
            if not any([x.add_residue(r) for x in self.Chains]):
                c=Chain.from_residue(r,parent_molecule=self)
                self.Chains.append(c)
