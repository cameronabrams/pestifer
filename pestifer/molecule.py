"""

.. module:: molecule
   :synopsis: Manages all molecules
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
from pidibble.pdbparse import PDBParser
from .mods import *
from .residue import Residue,ResidueList
from .segment import Segment,SegmentList
from .config import ConfigGetParam
from .bioassemb import BioAssemb

logger=logging.getLogger(__name__)

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
        li=pr['LINK']
        self.Links=LinkList([Link.from_pdbrecord(p) for p in li])
        self._MakeResidues()
        self._MakeLinks()
        self._MakeSegments()
        # self._MakeChains()
        self._MakeBiologicalAssemblies()

    @classmethod
    def clone(cls,inst):
        # nchains=len(inst.Chains)
        # if nchains>len(inst.available_chainIDs):
        #     logger.warning(f'There are not enough chain IDs available to clone this molecule')
        #     return
        copy=cls()
        copy.is_clone_of=inst
        copy.chainID_map=inst.generate_chainIDmap()

    def generate_chainIDmap(self):
        chainIDmap={}
        for ch in [x.chainID for x in self.Chains]:
            nch=self.available_chainIDs[0]
            chainIDmap[ch]=nch
            self.claim_chainID(nch)
        return chainIDmap


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

    def _MakeLinks(self):
        """for each Link, apply res1.linkTo(res2) 
        """
        for link in self.Links:
            link.residue1=self.Residues.get(chainID=link.chainID1,resseqnum=link.resseqnum1)
            link.residue2=self.Residues.get(chainID=link.chainID2,resseqnum=link.resseqnum2)
            a1=self.Atoms.get(chainID=link.chainID1,resseqnum=link.resseqnum1,name=link.name1)
            if hasattr(a1,'len'):
                a1=a1.get(altloc=link.altloc1)
            link.atom1=a1
            a2=self.Atoms.get(chainID=link.chainID2,resseqnum=link.resseqnum2,name=link.name2)
            if hasattr(a1,'len'):
                a2=a2.get(altloc=link.altloc2)
            link.atom2=a2
            link.residue1.linkTo(link.residue2,link)

    def _MakeSegments(self):
        self.Segments=SegmentList([])
        # make a single segment for each protein chain
        protein_res=self.Residues.get(segtype='PROTEIN')
        protein_chains=[]
        for r in protein_res:
            if not r.chainID in protein_chains:
                protein_chains.append(r.chainID)
        for pc in protein_chains:
            pc_res=self.Residues.get(segtype='PROTEIN',chain=pc)
            self.Segments.append(Segment.from_residue_list(pc_res,parent_molecule=self))
        


    # def _MakeChains(self):
    #     r=self.Residues[0]
    #     c=Chain.from_residue(r,parent_molecule=self)
    #     self.Chains=[]
    #     self.Chains.append(c)
    #     for r in self.Residues[1:]:
    #         if not any([x.add_residue(r) for x in self.Chains]):
    #             c=Chain.from_residue(r,parent_molecule=self)
    #             self.Chains.append(c)

    def _MakeBiologicalAssemblies(self):
        self.BioAssemb=[BioAssemb.from_asymmetricUnit()] # 0th BioAssemb is the asymmetric unit
        # each instance of REMARK.350.BIOMOLECULEn.TRANSFORMm (n, m integers 1, 2, ...)
        # has one or more transformations.
        P=self.pdb_entry.parsed
        barecs=[P[x] for x in P if ('REMARK.350.BIOMOLECULE' in x and 'TRANSFORM' in x)]
        for rec in barecs:
            self.BioAssemb.append(BioAssemb.from_PDBRecord(rec))

    # def claim_chainID(self,ch):
    #     self.claimed_chainIDs.append(ch)
    #     self.available_chainIDs.remove(ch)