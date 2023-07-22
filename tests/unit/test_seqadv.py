from pestifer.seqadv import Seqadv
import unittest

class TestSeqadv(unittest.TestCase):
    def setUp(self):
        self.test_records=[
            'SEQADV 4ZMJ ASN G  332  UNP  Q2N0S6    THR   330 ENGINEERED MUTATION            ',
            'SEQADV 5ZMJ CYS G  501  UNP  Q2N0S6    ALA   498 ENGINEERED MUTATION            ',
            'SEQADV 6ZMJ ARG G  508A UNP  Q2N0S6              EXPRESSION TAG                 ',
            'SEQADV 7ZMJ ARG G  509  UNP  Q2N0S6              EXPRESSION TAG                 '
        ]
        self.sq=[Seqadv(x) for x in self.test_records]

    def test_icode(self):
        expected_icodes=[' ',' ','A',' ']
        self.assertTrue(self.sq[i].iCode==expected_icodes[i] for i in range(len(expected_icodes)))

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
        newone=self.sq[0].clone(self.sq[0].chainID)
        self.assertTrue((not newone is self.sq[0]) and (newone==self.sq[0]))