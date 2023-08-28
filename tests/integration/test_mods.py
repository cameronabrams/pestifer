from pestifer.mods import Mutation, Missing, Seqadv, SSBond, Crot, Link
from pestifer.cifutil import CIFdict
from pidibble.pdbparse import PDBParser
from pestifer.bioassemb import Transform
from pestifer.config import Config
from pestifer.stringthings import ByteCollector
import unittest
# import pytest

class TestMutation(unittest.TestCase):
    def setUp(self):
        self.input_dict={
            'chainID':'A',
            'origresname':'PHE',
            'newresname':'TYR',
            'resseqnum':'123',
            'insertion':'A',
            'source':'USER'
        }
        self.shortcode='A:PHE,123A,TYR'
        return super().setUp()
    def test_init(self):
        m=Mutation(self.input_dict)
        self.assertTrue(type(m)==Mutation)
    def test_shortcode(self):
        m=Mutation(self.shortcode)
        self.assertTrue(m.chainID=='A' and m.resseqnumi=='123A')
    def test_clone(self):
        m=Mutation(self.input_dict)
        n=Mutation.clone(m,chainID=m.chainID)
        self.assertTrue(m==n and not m is n)
    def test_clone2(self):
        m=Mutation(self.input_dict)
        newchainid=chr((ord(m.chainID)+1)%26)
        n=Mutation.clone(m,chainID=newchainid)
        self.assertTrue(n.chainID==newchainid and n.resseqnum==m.resseqnum)

class TestMissing(unittest.TestCase):
    def setUp(self):
        self.cifdict=CIFdict({
            'pdb_ins_code':'?',
            'pdb_model_num':'1',
            'auth_comp_id':'GLU',
            'auth_asym_id':'G',
            'auth_seq_id':'185'
        })
        return super().setUp()
    def test_pdbrecord(self):
        p=PDBParser(PDBcode='4i53').parse()
        self.assertTrue('REMARK.465' in p.parsed)
        r=p.parsed['REMARK.465'].tables['MISSING'][0]
        m=Missing(r)
        self.assertEqual(m.model,'')
        self.assertEqual(m.resname,'VAL')
        self.assertEqual(m.chainID,'A')
        self.assertEqual(m.resseqnum,44)
        self.assertEqual(m.insertion,'')
        r=p.parsed['REMARK.465'].tables['MISSING'][14]
        m=Missing(r)
        self.assertEqual(m.model,'')
        self.assertEqual(m.resname,'GLY')
        self.assertEqual(m.chainID,'A')
        self.assertEqual(m.resseqnum,457)
        self.assertEqual(m.insertion,'B')
    def test_cifdict(self):
        m=Missing(self.cifdict)
        self.assertEqual(m.model,1)
        self.assertEqual(m.resname,'GLU')
        self.assertEqual(m.chainID,'G')
        self.assertEqual(m.resseqnum,185)
        self.assertEqual(m.insertion,' ')


class TestSSBond(unittest.TestCase):
    def test_shortcode(self):
        ss=SSBond('D_378-B_379')
        self.assertEqual(ss.chainID1,'D')
        self.assertEqual(ss.chainID2,'B')
        self.assertEqual(ss.resseqnum1,378)
        self.assertEqual(ss.resseqnum2,379)
    def test_psfgen(self):
        B=ByteCollector()
        ss=SSBond('D_378-B_379')
        transform=Transform()
        transform.chainIDmap={'D':'D','B':'E'}
        ss.write_TcL(B,transform)
        self.assertEqual(str(B),'patch DISU D:378 E:379\n')

class TestCrot(unittest.TestCase):
    def test_phi(self):
        for t in ['PHI','PSI','OMEGA']:
            cr=Crot.from_shortcode(f'{t},A,123,124,180.0')
            self.assertEqual(cr.angle,t)
            self.assertEqual(cr.chainID,'A')
            self.assertEqual(cr.resseqnum1,123)
            self.assertEqual(cr.resseqnum2,124)
            self.assertEqual(cr.degrees,180.0)
    def test_chi(self):
        for t in ['CHI1','CHI2']:
            cr=Crot.from_shortcode(f'{t},A,123,180.0')
            self.assertEqual(cr.angle,t)
            self.assertEqual(cr.chainID,'A')
            self.assertEqual(cr.resseqnum1,123)
            self.assertEqual(cr.resseqnum2,-1)
            self.assertEqual(cr.degrees,180.0)
    def test_glycan(self):
        cr=Crot.from_shortcode('GLYCAN,S,123,N,1123,O1,180.0')
        self.assertEqual(cr.angle,'GLYCAN')
        self.assertEqual(cr.segname,'S')
        self.assertEqual(cr.resseqnum1,123)
        self.assertEqual(cr.atom1,'N')
        self.assertEqual(cr.resseqnum2,1123)
        self.assertEqual(cr.atom2,'O1')
        self.assertEqual(cr.degrees,180.0)
    def test_link(self):
        cr=Crot.from_shortcode('LINK,S,123,N,T,1123,O1,180.0')
        self.assertEqual(cr.angle,'LINK')
        self.assertEqual(cr.segname1,'S')
        self.assertEqual(cr.segname2,'T')
        self.assertEqual(cr.resseqnum1,123)
        self.assertEqual(cr.atom1,'N')
        self.assertEqual(cr.resseqnum2,1123)
        self.assertEqual(cr.atom2,'O1')
        self.assertEqual(cr.degrees,180.0)
    def test_angleijk(self):
        cr=Crot.from_shortcode('ANGLEIJK,S,123,N,T,1123,O1,1123,C1,180.0')
        self.assertEqual(cr.angle,'ANGLEIJK')
        self.assertEqual(cr.segnamei,'S')
        self.assertEqual(cr.segnamejk,'T')
        self.assertEqual(cr.resseqnumi,123)
        self.assertEqual(cr.atomi,'N')
        self.assertEqual(cr.resseqnumj,1123)
        self.assertEqual(cr.atomj,'O1')
        self.assertEqual(cr.resseqnumk,1123)
        self.assertEqual(cr.atomk,'C1')
        self.assertEqual(cr.degrees,180.0)

class TestLink(unittest.TestCase):
    def test_link_create(self):
        self.config=Config('user_config.yaml')
        input_dict={
            'name1': 'NH',
            'resname1':'ASN',
            'chainID1':'A',
            'resseqnum1':1,
            'insertion1':'',
            'name2': 'C6',
            'resname2':'AMAN',
            'chainID2':'A',
            'resseqnum2':2,
            'insertion2':'',}
        l=Link(input_dict)
        l.map_attr('segtype1','resname1',self.config['Segtypes_by_Resnames'])
        l.map_attr('segtype2','resname2',self.config['Segtypes_by_Resnames'])
        self.assertEqual(l.segtype1,'PROTEIN')
        self.assertEqual(l.segtype2,'GLYCAN')