from pestifer.mods import ModsContainer, BaseMod,Atom,Mutation, Missing, Seqadv, SSBond, Crot
import unittest
import pytest

class TestBaseMod(unittest.TestCase):
    def test1(self):
        # required attributes are not in mutual exclusives
        class tbm(BaseMod):
            req_attr=['R1','R2','R3']
            opt_attr=['O1','O2']
            alt_attr=[('R1','O2')] # wrong!
        with pytest.raises(AssertionError) as E:
            a=tbm({'R1':1,'R2':2,'R3':3})
        self.assertEqual(E.type,AssertionError)
    def test2(self):
        # all required attributes have values
        class tbm(BaseMod):
            req_attr=['R1','R2','R3']
        with pytest.raises(AssertionError) as E:
            a=tbm({'R1':1,'R2':2})
        self.assertEqual(E.type,AssertionError)
    def test3(self):
        # all mutual exclusivity requirements are met
        class tbm(BaseMod):
            req_attr=['R1','R2','R3']
            opt_attr=['O1','O2']
            alt_attr=[('O1','O2')]
        with pytest.raises(AssertionError) as E:
            a=tbm({'R1':1,'R2':2,'R3':3,'O1':4,'O2':5})
        self.assertEqual(E.type,AssertionError)
    def test4(self):
        # all attributes with limited set of possible values have valid values
        class tbm(BaseMod):
            req_attr=['R1','R2','R3']
            opt_attr=['O1','O2']
            alt_attr=[('O1','O2')]
            attr_choices={'R1':[7,8,9]}
        with pytest.raises(AssertionError) as E:
            a=tbm({'R1':1,'R2':2,'R3':3})
        self.assertEqual(E.type,AssertionError)
    def test5(self):
        # for each optional attribute present, all required dependent attributes are also present
        class tbm(BaseMod):
            req_attr=['R1','R2','R3']
            opt_attr=['O1','O2']
            opt_attr_deps={'O1':['O2']}
        with pytest.raises(AssertionError) as E:
            a=tbm({'R1':1,'R2':2,'R3':3,'O1':4})
        self.assertEqual(E.type,AssertionError)
class TestMutation(unittest.TestCase):
    def setUp(self):
        self.input_dict={
            'chainID':'A',
            'origresname':'PHE',
            'newresname':'TYR',
            'resseqnum':'123',
            'insertion':'A',
            'resseqnumi':'123A'
        }
        self.shortcode='A:PHE,123A,TYR'
        return super().setUp()
    def test_init(self):
        m=Mutation(self.input_dict)
        self.assertTrue(type(m)==Mutation)
    def test_shortcode(self):
        m=Mutation.from_shortcode(self.shortcode)
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
        self.pdbrecord='REMARK 465     GLU G   185A  '
        self.cifdict={
            'pdb_ins_code':'?',
            'pdb_model_num':'1',
            'auth_comp_id':'GLU',
            'auth_asym_id':'G',
            'auth_seq_id':'185'
        }
        return super().setUp()
    def test_pdbrecord(self):
        m=Missing.from_pdbrecord(self.pdbrecord)
        self.assertEqual(m.model,'')
        self.assertEqual(m.resname,'GLU')
        self.assertEqual(m.chainID,'G')
        self.assertEqual(m.resseqnum,185)
        self.assertEqual(m.insertion,'A')
    def test_cifdict(self):
        m=Missing.from_cifdict(self.cifdict)
        self.assertEqual(m.model,1)
        self.assertEqual(m.resname,'GLU')
        self.assertEqual(m.chainID,'G')
        self.assertEqual(m.resseqnum,185)
        self.assertEqual(m.insertion,' ')
    def test_psfgenlines(self):
        m=Missing.from_pdbrecord(self.pdbrecord)
        testln=m.psfgen_lines()
        expectedln=[f'     residue {m.resname} {m.resseqnum}{m.insertion} {m.chainID}']
        self.assertEqual(testln,expectedln)
    def test_clone(self):
        m=Missing.from_pdbrecord(self.pdbrecord)
        n=Missing.clone(m,chainID=m.chainID)
        self.assertTrue(m==n and not m is n)
    def test_clone2(self):
        m=Missing.from_pdbrecord(self.pdbrecord)
        newchainid=chr((ord(m.chainID)+1)%26)
        n=Missing.clone(m,chainID=newchainid)
        self.assertTrue(n.chainID==newchainid and n.resseqnum==m.resseqnum)

class TestAtom(unittest.TestCase):
    def setUp(self):
        return super().setUp()
    def test_pdbrecord(self):
        pr1='ATOM    981  N   GLU G 164A     10.462 130.792 -22.224  1.00 70.92           N  '
        a1=Atom.from_pdbrecord(pdbrecord=pr1)
        self.assertEqual(a1.serial,981)
        self.assertEqual(a1.name,'N')
        self.assertEqual(a1.resname,'GLU')
        self.assertEqual(a1.chainID,'G')
        self.assertEqual(a1.resseqnum,164)
        self.assertEqual(a1.insertion,'A')
        self.assertEqual(a1.x,10.462)
        self.assertEqual(a1.y,130.792)
        self.assertEqual(a1.occ,1.00)
        self.assertEqual(a1.beta,70.92)
        self.assertEqual(a1.elem,'N')
        self.assertEqual(a1.charge,'')
    def test_clone2(self):
        pr1='ATOM    981  N   GLU G 164A     10.462 130.792 -22.224  1.00 70.92           N  '
        a=Atom.from_pdbrecord(pr1)
        newchainid=chr((ord(a.chainID)+1)%26)
        b=Atom.clone(a,chainID=newchainid)
        self.assertTrue(b.chainID==newchainid and a.serial==b.serial)
    def test_pdbline(self):
        pr1='ATOM    981  N   GLU G 164A     10.462 130.792 -22.224  1.00 70.92           N  '
        a=Atom.from_pdbrecord(pr1)
        self.assertEqual(a.pdb_line(),pr1)

class TestSeqadv(unittest.TestCase):
    def setUp(self):
        self.test_records=[
            'SEQADV 4ZMJ ASN G  332  UNP  Q2N0S6    THR   330 ENGINEERED MUTATION            ',
            'SEQADV 5ZMJ CYS G  501  UNP  Q2N0S6    ALA   498 ENGINEERED MUTATION            ',
            'SEQADV 6ZMJ ARG G  508A UNP  Q2N0S6              EXPRESSION TAG                 ',
            'SEQADV 7ZMJ ARG G  509  UNP  Q2N0S6              EXPRESSION TAG                 '
        ]
        self.sq=[Seqadv.from_pdbrecord(x) for x in self.test_records]

    def test_insertions(self):
        expected_icodes=[' ',' ','A',' ']
        self.assertTrue(self.sq[i].insertion==expected_icodes[i] for i in range(len(expected_icodes)))

    def test_pdbcode(self):
        expected_pdbcodes=['4ZMJ','5ZMJ','6ZMJ','7ZMJ']
        self.assertTrue(self.sq[i].idCode==expected_pdbcodes[i] for i in range(len(expected_pdbcodes)))
    
    def test_db(self):
        expected_dbres=['THR','ALA','   ','   ']
        self.assertTrue(self.sq[i].dbRes==expected_dbres[i] for i in range(len(expected_dbres)))

    def test_pdbline(self):
        for i in range(len(self.sq)):
            print('input:     [',self.test_records[i],']')
            print('pdb_line():[',self.sq[i].pdb_line(),']')
        self.assertTrue(any([self.test_records[i]==self.sq[i].pdb_line() for i in range(len(self.sq))]))
    
    def test_clone(self):
        newchainid=chr((ord(self.sq[0].chainID)+1)%26)
        b=Seqadv.clone(self.sq[0],chainID=newchainid)
        self.assertTrue(b.chainID==newchainid and self.sq[0].idCode==b.idCode)

    def test_mutations(self):
        mutlist=[Mutation.from_seqdav(x) for x in self.sq]
        for i in [0,1]:
            self.assertEqual(mutlist[i].chainID,self.sq[i].chainID)
            self.assertEqual(mutlist[i].origresname,self.sq[i]. resname)
            self.assertEqual(mutlist[i].newresname,self.sq[i].dbRes)
            self.assertEqual(mutlist[i].insertion,self.sq[i].insertion)
        self.assertFalse(mutlist[2])
        self.assertFalse(mutlist[3])

class TestSSBond(unittest.TestCase):
    def setUp(self):
        self.pdbline1="SSBOND   4 CYS D  378    CYS B  379                          1555   2655  2.04"
        self.pdbline2="SSBOND   8 CYS D  378    CYS B  178                          2555   2655  2.04"
        return super().setUp()
    def test_pdbrecord(self):
        ss=SSBond.from_pdbrecord(self.pdbline1)
        self.assertEqual(ss.serial_number,4)
        self.assertEqual(ss.resname1,'CYS')
        self.assertEqual(ss.resname2,'CYS')
        self.assertEqual(ss.chainID1,'D')
        self.assertEqual(ss.chainID2,'B')
        self.assertEqual(ss.resseqnum1,378)
        self.assertEqual(ss.resseqnum2,379)
        self.assertEqual(ss.sym1,'1555')
        self.assertEqual(ss.sym2,'2655')
        self.assertEqual(ss.length,2.04)
    def test_shortcode(self):
        ss=SSBond.from_shortcode('D_378-B_379')
        self.assertEqual(ss.chainID1,'D')
        self.assertEqual(ss.chainID2,'B')
        self.assertEqual(ss.resseqnum1,378)
        self.assertEqual(ss.resseqnum2,379)
    def test_psfgen(self):
        ss=SSBond.from_shortcode('D_378-B_379')
        self.assertEqual(ss.psfgen_lines(),['patch DISU D:378 B:379'])
    def test_pdbline(self):
        ss=SSBond.from_pdbrecord(self.pdbline2)
        self.assertEqual(ss.pdb_line(),self.pdbline2)
    def test_clone(self):
        ss=SSBond.from_pdbrecord(self.pdbline2)
        ss2=SSBond.clone(ss,chainID2='X')
        self.assertEqual(ss.chainID1,ss2.chainID1)
        self.assertEqual(ss2.chainID2,'X')

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