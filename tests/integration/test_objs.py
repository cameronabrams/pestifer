from pestifer.objs.mutation import Mutation
from pestifer.objs.ter import Ter
from pestifer.objs.seqadv import Seqadv
from pestifer.objs.ssbond import SSBond
from pestifer.objs.crot import Crot
from pestifer.objs.link import Link
from pestifer.objs.insertion import Insertion
from pestifer.util.cifutil import CIFdict, CIFload
from pidibble.pdbparse import PDBParser
from pestifer.bioassemb import Transform
from pestifer.config import Config, segtype_of_resname
from pestifer.stringthings import ByteCollector
from pidibble.pdbparse import PDBParser
import unittest
# import pytest

class TestTer(unittest.TestCase):
    def test_ter_fromdict(self):
        d={'serial':-1,'chainID':'A','resname':'ABC','resseqnum':1,'insertion':''}
        t=Ter(d)
        self.assertEqual(type(t),Ter)
        self.assertEqual(d['serial'],t.serial)
        self.assertEqual(d['chainID'],t.chainID)
    def test_ter_frompdbrecord(self):
        p=PDBParser(PDBcode='4zmj').parse()
        pr=p.parsed['TER'][0]
        t=Ter(pr)
        self.assertEqual(type(t),Ter)

class TestSeqadv(unittest.TestCase):
    def test_seqadv_fromdict(self):
        d={a:0 for a in Seqadv.req_attr}
        d['typekey']='user'
        s=Seqadv(d)
        self.assertEqual(type(s),Seqadv)
        self.assertEqual(d['idCode'],s.idCode)
    def test_seqadv_frompdbrecord(self):
        p=PDBParser(PDBcode='4zmj').parse()
        pr=p.parsed['SEQADV'][0]
        s=Seqadv(pr)
        self.assertEqual(type(s),Seqadv)
    def test_seqadv_fromcifdict(self):
        p=CIFload('4zmj')
        obj=p.getObj(Seqadv.mmCIF_name)
        d=CIFdict(obj,0)
        s=Seqadv(d)
        self.assertEqual(type(s),Seqadv)


class TestMutation(unittest.TestCase):
    def test_mutation_fromdict(self):
        input_dict={
            'chainID':'A',
            'origresname':'PHE',
            'newresname':'TYR',
            'resseqnum':'123',
            'insertion':'A',
            'typekey':'user'
        }
        m=Mutation(input_dict)
        self.assertTrue(type(m)==Mutation)

    def test_mutation_fromshortcode(self):
        c=Config()
        shortcode='A:PHE,123A,TYR'
        m=Mutation(shortcode)
        self.assertTrue(m.chainID=='A' and m.resseqnum==123 and m.insertion=='A')
        another_format='A:F11T' #must call Config() to set up the dicts
        m=Mutation(another_format)
        self.assertEqual(m.chainID,'A')
        self.assertEqual(m.origresname,'PHE')
        self.assertEqual(m.newresname,'THR')
        self.assertEqual(m.resseqnum,11)
        self.assertEqual(m.insertion,'')

    def test_mutation_fromseqadv(self):
        p=PDBParser(PDBcode='4zmj').parse()
        pr=p.parsed['SEQADV'][0]
        s=Seqadv(pr)
        m=Mutation(s)
        self.assertEqual(type(m),Mutation)

    def test_mutation_clone_exact(self):
        input_dict={
            'chainID':'A',
            'origresname':'PHE',
            'newresname':'TYR',
            'resseqnum':'123',
            'insertion':'A',
            'typekey':'user'
        }
        m=Mutation(input_dict)
        n=m.clone(chainID=m.chainID)
        self.assertTrue(m==n and not m is n)
        
    def test_mutation_clone_nonexact(self):
        input_dict={
            'chainID':'A',
            'origresname':'PHE',
            'newresname':'TYR',
            'resseqnum':'123',
            'insertion':'A',
            'typekey':'user'
        }
        m=Mutation(input_dict)
        newchainid=chr((ord(m.chainID)+1)%26)
        n=Mutation.clone(m,chainID=newchainid)
        self.assertTrue(n.chainID==newchainid and n.resseqnum==m.resseqnum)


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
            cr=Crot(f'{t},A,123,124,180.0')
            self.assertEqual(cr.angle,t)
            self.assertEqual(cr.chainID,'A')
            self.assertEqual(cr.resseqnum1,123)
            self.assertEqual(cr.resseqnum2,124)
            self.assertEqual(cr.degrees,180.0)
    def test_chi(self):
        for t in ['CHI1','CHI2']:
            cr=Crot(f'{t},A,123,180.0')
            self.assertEqual(cr.angle,t)
            self.assertEqual(cr.chainID,'A')
            self.assertEqual(cr.resseqnum1,123)
            self.assertEqual(cr.resseqnum2,-1)
            self.assertEqual(cr.degrees,180.0)
    def test_angleijk(self):
        cr=Crot('ANGLEIJK,S,123,N,T,1123,O1,1123,C1,180.0')
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
        l.map_attr('segtype1','resname1',segtype_of_resname)
        l.map_attr('segtype2','resname2',segtype_of_resname)
        self.assertEqual(l.segtype1,'protein')
        self.assertEqual(l.segtype2,'glycan')

class TestInsertion(unittest.TestCase):
    def test_insertion_create(self):
        I=Insertion('C,111,GRETA')
        self.assertEqual(I.chainID,'C')