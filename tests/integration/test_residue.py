import unittest
from pestifer.residue import Residue, EmptyResidue
from pestifer.config import Config, res_123
from pestifer.molecule import Molecule
from pestifer.util.cifutil import CIFdict, CIFload
from pestifer.objs.insertion import InsertionList, Insertion
from io import StringIO
from pidibble.pdbparse import PDBParser
import yaml

class TestResidue(unittest.TestCase):
    def get_source_dict(self,pdbid):
        source=f"""
source:
    biological_assembly: 0
    exclude: {{}}
    file_format: PDB
    id: {pdbid}
    sequence:
        fix_conflicts: true
        fix_engineered_mutations: true
        include_terminal_loops: false
    loops:
        declash:
        maxcycles: 20
        min_loop_length: 4
        sac_res_name: GLY
"""
        f=StringIO(source)
        directive=yaml.safe_load(f)
        return directive
    def test_from_atom(self):
        c=Config()
        directive=self.get_source_dict('1gc1')
        m=Molecule(source=directive["source"])
        au=m.asymmetric_unit
        residues=[]
        a=au.atoms[0]
        r=Residue(a)
        self.assertEqual(a.chainID,r.chainID)
        residues.append(r)
        for a in au.atoms[1:]:
            if not any([x.add_atom(a) for x in residues]):
                r=Residue(a)
                residues.append(r)
        self.assertEqual(len(residues),1538)

    def test_from_shortcode(self):
        c=Config()
        shortcode='1:C:R525A'
        R=Residue(shortcode)
        self.assertEqual(R.chainID,'C')
        self.assertEqual(R.resseqnum,525)
        self.assertEqual(R.model,1)
        self.assertEqual(R.resname,'ARG')
        self.assertEqual(R.insertion,'A')

    def test_segtypes(self):
        c=Config()
        directive=self.get_source_dict('4tvp')
        m=Molecule(source=directive["source"])
        au=m.asymmetric_unit
        p=au.residues.filter(segtype='protein')
        pc=[]
        for r in p:
            if not r.chainID in pc:
                pc.append(r.chainID)
        self.assertEqual(pc,['G', 'B', 'L', 'H', 'D', 'E'])
        i=au.residues.filter(segtype='ion')
        ic=[]
        for r in i:
            if not r.chainID in ic:
                ic.append(r.chainID)
        self.assertEqual(ic,['X', 'Y', 'Z'])
        w=au.residues.filter(segtype='water')
        wc=[]
        for r in w:
            if not r.chainID in wc:
                wc.append(r.chainID)
        self.assertEqual(wc,['a', 'b', 'c', 'd', 'e', 'f'])
        g=au.residues.filter(segtype='glycan')
        gc=[]
        for r in g:
            if not r.chainID in gc:
                gc.append(r.chainID)
        self.assertEqual(gc,['A', 'C', 'F', 'I', 'J', 'K', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W'])
        l=au.residues.filter(segtype='ligand')
        lc=[]
        for r in l:
            if not r.chainID in lc:
                lc.append(r.chainID)
        self.assertEqual(lc,[])
        o=au.residues.filter(segtype='other')
        oc=[]
        for r in o:
            if not r.chainID in oc:
                oc.append(r.chainID)
        self.assertEqual(oc,[])

    def test_residuelist(self):
        c=Config()
        directive=self.get_source_dict('6m0j')
        m=Molecule(source=directive["source"])
        au=m.asymmetric_unit
        r=au.residues.get_residue(resseqnum=427,chainID='A')
        self.assertEqual(r.resname,'ASP')
        a=au.residues.get_atom('OD2',resseqnum=427,chainID='A')
        self.assertEqual(a.serial,3302)

    def test_insertions(self):
        c=Config()
        directive=self.get_source_dict('6m0j')
        m=Molecule(source=directive["source"])
        au=m.asymmetric_unit
        R=au.residues
        orig_numres=len(R)
        I=InsertionList([Insertion('A,427,GRETA')])
        R.apply_insertions(I)
        self.assertEqual(len(R),orig_numres+5)
        for i,n in zip('ABCDE','GRETA'):
            r=R.get(chainID='A',resseqnum=427,insertion=i)
            self.assertEqual(r.resname,res_123[n])
            self.assertEqual(r.resolved,False)
            self.assertEqual(r.segtype,'protein')

    def test_relations(self):
        c=Config()
        R_low=Residue('1:C:R525')
        R_hi=Residue('1:C:R525A')
        self.assertTrue(R_low<R_hi)
        R_same=Residue('1:C:R525')
        self.assertTrue(R_low<=R_same)
        str_hi='525A'
        self.assertTrue(R_low<str_hi)
        self.assertTrue(str_hi>R_low)
        str_same='525'
        self.assertFalse(R_low<str_same)
        self.assertTrue(R_low.same_resid(str_same))
    
class TestEmptyResidue(unittest.TestCase):
    def setUp(self):
        return super().setUp()
    def test_pdbrecord(self):
        p=PDBParser(PDBcode='4i53').parse()
        self.assertTrue('REMARK.465' in p.parsed)
        r=p.parsed['REMARK.465'].tables['MISSING'][0]
        m=EmptyResidue(r)
        self.assertEqual(m.model,'')
        self.assertEqual(m.resname,'VAL')
        self.assertEqual(m.chainID,'A')
        self.assertEqual(m.resseqnum,44)
        self.assertEqual(m.insertion,'')
        r=p.parsed['REMARK.465'].tables['MISSING'][14]
        m=EmptyResidue(r)
        self.assertEqual(m.model,'')
        self.assertEqual(m.resname,'GLY')
        self.assertEqual(m.chainID,'A')
        self.assertEqual(m.resseqnum,457)
        self.assertEqual(m.insertion,'B')
    def test_cifdict(self):
        p=CIFload('4i53')
        obj=p.getObj('pdbx_unobs_or_zero_occ_residues')
        m=EmptyResidue(CIFdict(obj,0))
        self.assertEqual(m.model,1)
        self.assertEqual(m.resname,'VAL')
        self.assertEqual(m.chainID,'A')
        self.assertEqual(m.resseqnum,1)
        self.assertEqual(m.insertion,'')
        
    def test_shortcode(self):
        c=Config()
        shortcode='2:X:W427'
        E=EmptyResidue(shortcode)
        self.assertEqual(E.model,2)
        self.assertEqual(E.chainID,'X')
        self.assertEqual(E.resseqnum,427)
        self.assertEqual(E.resname,'TRP')
        self.assertEqual(E.insertion,'')
        shortcode='Y:S111'
        E=EmptyResidue(shortcode)
        self.assertEqual(E.model,1)
        self.assertEqual(E.chainID,'Y')
        self.assertEqual(E.resseqnum,111)
        self.assertEqual(E.resname,'SER')
        self.assertEqual(E.insertion,'')
        shortcode='Y666'
        E=EmptyResidue(shortcode)
        self.assertEqual(E.model,1) # default
        self.assertEqual(E.chainID,'A') # default
        self.assertEqual(E.resseqnum,666)
        self.assertEqual(E.resname,'TYR')
        self.assertEqual(E.insertion,'')
