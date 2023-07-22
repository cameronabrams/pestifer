from pestifer.mods import ModsContainer, Mutation
import unittest

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
    def test_init(self):
        m=Mutation(self.input_dict)
        self.assertTrue(type(m)==Mutation)
    def test_shortcode(self):
        m=Mutation({'shortcode':self.shortcode})
        self.assertTrue(m.chainID=='A' and m.resseqnumi=='123A')
    def test_clone(self):
        m=Mutation(self.input_dict)
        n=m.clone()
        self.assertTrue(m==n and not m is n)
