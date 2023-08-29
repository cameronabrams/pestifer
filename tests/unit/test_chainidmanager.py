import unittest
from pestifer.chainids import ChainIDManager

class TestChainIDManager(unittest.TestCase):
    def test_pdb(self):
        cm=ChainIDManager()
        UC=[chr(x) for x in range(ord('A'),ord('A')+26)]
        LC=[chr(x) for x in range(ord('a'),ord('a')+26)]
        DC=[str(x) for x in range(10)]
        res=UC+LC+DC
        self.assertEqual(cm.OrderedSupply,res)
    def test_cif(self):
        cm=ChainIDManager(format='mmCIF')
        self.assertEqual(len(cm.OrderedSupply),27*26)
        first_few=['A','B','C','D']
        self.assertEqual(cm.OrderedSupply[:len(first_few)],first_few)
        items_26_up=['AA','BA','CA','DA']
        self.assertEqual(cm.OrderedSupply[26:26+len(items_26_up)],items_26_up)
        items_52_up=['AB','BB','CB','DB']
        self.assertEqual(cm.OrderedSupply[52:52+len(items_52_up)],items_52_up)
