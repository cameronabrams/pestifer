from pestifer.mods import Atom,Mutation, Missing, Seqadv, SSBond, Crot
from pestifer.basemod import BaseMod, ModList, StateBound
from pidibble.pdbparse import PDBParser
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

    def test6(self):
        class tbm(BaseMod):
            req_attr=['a','b']
        id1={'a':1,'b':2}
        bm1=tbm(id1)
        bm2=BaseMod({})
        self.assertFalse(bm1==bm2)

    def test_binary_logic(self):
        class tbm(BaseMod):
            req_attr=['a','b']
        bm1=tbm({'a':1,'b':2})
        bm2=tbm({'a':0,'b':2})
        self.assertTrue(bm2<bm1)
        bm1=tbm({'a':1,'b':2})
        bm2=tbm({'a':0,'b':3})
        self.assertFalse(bm2<bm1)
        self.assertTrue(bm2.weak_lt(bm1,['a']))
        self.assertFalse(bm2.weak_lt(bm1,['b']))
        self.assertFalse(bm2.weak_lt(bm1,[]))
        bm2=tbm({'a':1,'b':2})
        self.assertTrue(bm1==bm2)
        self.assertFalse(bm2<bm1)
        self.assertFalse(bm1<bm2)
        class xbm(BaseMod):
            req_attr=['a','b']
        bm3=xbm({'a':1,'b':2})
        self.assertFalse(bm3==bm1)
        self.assertFalse(bm3<bm1)
        self.assertFalse(bm1<bm3)

    def test_uniquify(self):
        L=ModList([])
        class tbm(BaseMod):
            req_attr=['a','b']
        L.append(tbm({'a':1,'b':1})) # 0
        L.append(tbm({'a':1,'b':1})) # 1
        L.append(tbm({'a':2,'b':1})) # 2
        L.append(tbm({'a':3,'b':1})) # 3
        L.append(tbm({'a':4,'b':1})) # 4
        L.append(tbm({'a':4,'b':1})) # 5
        L.append(tbm({'a':4,'b':1})) # 6
        L.append(tbm({'a':1,'b':1})) # 7
        L.puniquify(['a'])
        self.assertEqual(L[0].a,1)
        self.assertTrue(not hasattr(L[0],'_ORIGINAL_'))
        self.assertEqual(L[1].a,5)
        self.assertTrue(hasattr(L[1],'_ORIGINAL_'))
        self.assertEqual(L[1]._ORIGINAL_['a'],1)        
        self.assertEqual(L[5].a,6)
        self.assertTrue(hasattr(L[5],'_ORIGINAL_'))
        self.assertEqual(L[5]._ORIGINAL_['a'],4)        
        self.assertEqual(L[6].a,7)
        self.assertTrue(hasattr(L[6],'_ORIGINAL_'))
        self.assertEqual(L[6]._ORIGINAL_['a'],4)        
        self.assertEqual(L[7].a,8)
        self.assertTrue(hasattr(L[7],'_ORIGINAL_'))
        self.assertEqual(L[7]._ORIGINAL_['a'],1)        

    def test_statebounds(self):
        L=ModList([])
        class tbm(BaseMod):
            req_attr=['a','b']
        L.append(tbm({'a':[],'b':1})) # 0
        L.append(tbm({'a':[],'b':1})) # 1
        L.append(tbm({'a':[2,3,4],'b':1})) # 2
        L.append(tbm({'a':[5,6],'b':1})) # 3
        L.append(tbm({'a':[7,8,9,10],'b':1})) # 4
        L.append(tbm({'a':[],'b':1})) # 5
        L.append(tbm({'a':[11,12],'b':1})) # 6
        L.append(tbm({'a':[13],'b':1})) # 7
        b=L.state_bounds(lambda x: 'RESOLVED' if len(x.a)>0 else 'MISSING')
        self.assertEqual(b[0],StateBound({'state':'MISSING','bounds':[0,1]}))    
        self.assertEqual(b[1],StateBound({'state':'RESOLVED','bounds':[2,4]}))   
        self.assertEqual(b[2],StateBound({'state':'MISSING','bounds':[5,5]}))
        self.assertEqual(b[3],StateBound({'state':'RESOLVED','bounds':[6,7]}))
        r=[c.pstr() for c in b]
        self.assertEqual(r,['MISSING(2)','RESOLVED(3)','MISSING(1)','RESOLVED(2)'])

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
"""  
REMARK 465                                                                      
REMARK 465   M RES C SSSEQI                                                     
REMARK 465     VAL A    44                                                      
REMARK 465     TRP A    45                                                      
REMARK 465     LYS A    46                                                      
REMARK 465     GLU A    47                                                      
REMARK 465     ALA A    48                                                      
REMARK 465     ASN A   317                                                      
REMARK 465     GLY A   318                                                      
REMARK 465     GLY A   319                                                      
REMARK 465     SER A   320                                                      
REMARK 465     GLY A   321                                                      
REMARK 465     SER A   322                                                      
REMARK 465     GLY A   323                                                      
REMARK 465     GLY A   324                                                      
REMARK 465     GLY A   457A                                                     
REMARK 465     GLY A   457B
"""
class TestMissing(unittest.TestCase):
    def setUp(self):
        self.cifdict={
            'pdb_ins_code':'?',
            'pdb_model_num':'1',
            'auth_comp_id':'GLU',
            'auth_asym_id':'G',
            'auth_seq_id':'185'
        }
        return super().setUp()
    def test_pdbrecord(self):
        p=PDBParser(PDBcode='4i53').parse()
        self.assertTrue('REMARK.465' in p.parsed)
        r=p.parsed['REMARK.465'].tables['MISSING'][0]
        m=Missing.from_pdbrecord(r)
        self.assertEqual(m.model,'')
        self.assertEqual(m.resname,'VAL')
        self.assertEqual(m.chainID,'A')
        self.assertEqual(m.resseqnum,44)
        self.assertEqual(m.insertion,'')
        r=p.parsed['REMARK.465'].tables['MISSING'][14]
        m=Missing.from_pdbrecord(r)
        self.assertEqual(m.model,'')
        self.assertEqual(m.resname,'GLY')
        self.assertEqual(m.chainID,'A')
        self.assertEqual(m.resseqnum,457)
        self.assertEqual(m.insertion,'B')
    def test_cifdict(self):
        m=Missing.from_cifdict(self.cifdict)
        self.assertEqual(m.model,1)
        self.assertEqual(m.resname,'GLU')
        self.assertEqual(m.chainID,'G')
        self.assertEqual(m.resseqnum,185)
        self.assertEqual(m.insertion,' ')
#     def test_psfgenlines(self):
#         m=Missing.from_pdbrecord(self.pdbrecord)
#         testln=m.psfgen_lines()
#         expectedln=[f'     residue {m.resname} {m.resseqnum}{m.insertion} {m.chainID}']
#         self.assertEqual(testln,expectedln)
#     def test_clone(self):
#         m=Missing.from_pdbrecord(self.pdbrecord)
#         n=Missing.clone(m,chainID=m.chainID)
#         self.assertTrue(m==n and not m is n)
#     def test_clone2(self):
#         m=Missing.from_pdbrecord(self.pdbrecord)
#         newchainid=chr((ord(m.chainID)+1)%26)
#         n=Missing.clone(m,chainID=newchainid)
#         self.assertTrue(n.chainID==newchainid and n.resseqnum==m.resseqnum)

# class TestAtom(unittest.TestCase):
#     def setUp(self):
#         return super().setUp()
#     def test_pdbrecord(self):
#         pr1='ATOM    981  N   GLU G 164A     10.462 130.792 -22.224  1.00 70.92           N  '
#         a1=Atom.from_pdbrecord(pdbrecord=pr1)
#         self.assertEqual(a1.serial,981)
#         self.assertEqual(a1.name,'N')
#         self.assertEqual(a1.resname,'GLU')
#         self.assertEqual(a1.chainID,'G')
#         self.assertEqual(a1.resseqnum,164)
#         self.assertEqual(a1.insertion,'A')
#         self.assertEqual(a1.x,10.462)
#         self.assertEqual(a1.y,130.792)
#         self.assertEqual(a1.occ,1.00)
#         self.assertEqual(a1.beta,70.92)
#         self.assertEqual(a1.elem,'N')
#         self.assertEqual(a1.charge,'')
#     def test_clone2(self):
#         pr1='ATOM    981  N   GLU G 164A     10.462 130.792 -22.224  1.00 70.92           N  '
#         a=Atom.from_pdbrecord(pr1)
#         newchainid=chr((ord(a.chainID)+1)%26)
#         b=Atom.clone(a,chainID=newchainid)
#         self.assertTrue(b.chainID==newchainid and a.serial==b.serial)
#     def test_pdbline(self):
#         pr1='ATOM    981  N   GLU G 164A     10.462 130.792 -22.224  1.00 70.92           N  '
#         a=Atom.from_pdbrecord(pr1)
#         self.assertEqual(a.pdb_line(),pr1)

# class TestSeqadv(unittest.TestCase):
#     def setUp(self):
#         self.test_records=[
#             'SEQADV 4ZMJ ASN G  332  UNP  Q2N0S6    THR   330 ENGINEERED MUTATION            ',
#             'SEQADV 5ZMJ CYS G  501  UNP  Q2N0S6    ALA   498 ENGINEERED MUTATION            ',
#             'SEQADV 6ZMJ ARG G  508A UNP  Q2N0S6              EXPRESSION TAG                 ',
#             'SEQADV 7ZMJ ARG G  509  UNP  Q2N0S6              EXPRESSION TAG                 '
#         ]
#         self.sq=[Seqadv.from_pdbrecord(x) for x in self.test_records]

#     def test_insertions(self):
#         expected_icodes=[' ',' ','A',' ']
#         self.assertTrue(self.sq[i].insertion==expected_icodes[i] for i in range(len(expected_icodes)))

#     def test_pdbcode(self):
#         expected_pdbcodes=['4ZMJ','5ZMJ','6ZMJ','7ZMJ']
#         self.assertTrue(self.sq[i].idCode==expected_pdbcodes[i] for i in range(len(expected_pdbcodes)))
    
#     def test_db(self):
#         expected_dbres=['THR','ALA','   ','   ']
#         self.assertTrue(self.sq[i].dbRes==expected_dbres[i] for i in range(len(expected_dbres)))

#     def test_pdbline(self):
#         for i in range(len(self.sq)):
#             print('input:     [',self.test_records[i],']')
#             print('pdb_line():[',self.sq[i].pdb_line(),']')
#         self.assertTrue(any([self.test_records[i]==self.sq[i].pdb_line() for i in range(len(self.sq))]))
    
#     def test_clone(self):
#         newchainid=chr((ord(self.sq[0].chainID)+1)%26)
#         b=Seqadv.clone(self.sq[0],chainID=newchainid)
#         self.assertTrue(b.chainID==newchainid and self.sq[0].idCode==b.idCode)

#     def test_mutations(self):
#         mutlist=[Mutation.from_seqdav(x) for x in self.sq]
#         for i in [0,1]:
#             self.assertEqual(mutlist[i].chainID,self.sq[i].chainID)
#             self.assertEqual(mutlist[i].origresname,self.sq[i]. resname)
#             self.assertEqual(mutlist[i].newresname,self.sq[i].dbRes)
#             self.assertEqual(mutlist[i].insertion,self.sq[i].insertion)
#         self.assertFalse(mutlist[2])
#         self.assertFalse(mutlist[3])

class TestSSBond(unittest.TestCase):
    # def test_pdbrecord(self):
    #     ss=SSBond.from_pdbrecord(self.pdbline1)
    #     self.assertEqual(ss.serial_number,4)
    #     self.assertEqual(ss.resname1,'CYS')
    #     self.assertEqual(ss.resname2,'CYS')
    #     self.assertEqual(ss.chainID1,'D')
    #     self.assertEqual(ss.chainID2,'B')
    #     self.assertEqual(ss.resseqnum1,378)
    #     self.assertEqual(ss.resseqnum2,379)
    #     self.assertEqual(ss.sym1,'1555')
    #     self.assertEqual(ss.sym2,'2655')
    #     self.assertEqual(ss.length,2.04)
    def test_shortcode(self):
        ss=SSBond.from_shortcode('D_378-B_379')
        self.assertEqual(ss.chainID1,'D')
        self.assertEqual(ss.chainID2,'B')
        self.assertEqual(ss.resseqnum1,378)
        self.assertEqual(ss.resseqnum2,379)
    def test_psfgen(self):
        ss=SSBond.from_shortcode('D_378-B_379')
        self.assertEqual(ss.psfgen_lines(),['patch DISU D:378 B:379'])
    # def test_pdbline(self):
    #     ss=SSBond.from_pdbrecord(self.pdbline2)
    #     self.assertEqual(ss.pdb_line(),self.pdbline2)
    # def test_clone(self):
    #     ss=SSBond.from_pdbrecord(self.pdbline2)
    #     ss2=SSBond.clone(ss,chainID2='X')
    #     self.assertEqual(ss.chainID1,ss2.chainID1)
    #     self.assertEqual(ss2.chainID2,'X')

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
