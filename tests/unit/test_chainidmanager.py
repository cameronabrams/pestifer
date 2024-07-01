import unittest
from pestifer.chainidmanager import ChainIDManager

class TestChainIDManager(unittest.TestCase):
    def test_pdb(self):
        cm=ChainIDManager()
        UC=[chr(x) for x in range(ord('A'),ord('A')+26)]
        LC=[chr(x) for x in range(ord('a'),ord('a')+26)]
        DC=[str(x) for x in range(10)]
        res=UC+LC+DC
        self.assertEqual(cm.Unused,res)
    def test_cif(self):
        cm=ChainIDManager(format='mmCIF')
        self.assertEqual(len(cm.Unused),27*26)
        first_few=['A','B','C','D']
        self.assertEqual(cm.Unused[:len(first_few)],first_few)
        items_26_up=['AA','BA','CA','DA']
        self.assertEqual(cm.Unused[26:26+len(items_26_up)],items_26_up)
        items_52_up=['AB','BB','CB','DB']
        self.assertEqual(cm.Unused[52:52+len(items_52_up)],items_52_up)

    def test_register(self):
        cm=ChainIDManager()
        for c in ['A','B','x']:
            cm.check(c)
        self.assertTrue(not any([x in cm.Unused for x in ['A','B','x']]))
        for c in ['A','B','x']:
            self.assertTrue(c in cm.Used)

    def test_unregister(self):
        cm=ChainIDManager()
        for c in ['A','B','x']:
            cm.check(c)
        cm.unregister_chain('A')
        self.assertFalse('A' in cm.Used)
        for c in ['B','x']:
            self.assertTrue(c in cm.Used)

    def test_reserved(self):
        cm=ChainIDManager(transform_reserves={'A':['G','Q'],'B':['H','R']})
        self.assertEqual(cm.ReservedUnused,['G','Q','H','R'])
        cm.check('A')
        cm.check('B')
        c=cm.next_unused_chainID()
        self.assertEqual(c,'C')
        d=cm.next_reserved_chainID('A')
        self.assertEqual(d,'G')
        d=cm.next_reserved_chainID('A')
        self.assertEqual(d,'Q')

    def test_remap_chainIDs(self):
        cm=ChainIDManager()
        self.assertTrue(cm.remap=={})
        cm=ChainIDManager(remap={'A':'X','B':'Y','C':'Z'})
        self.assertEqual(cm.remap['A'],'X')
        
    def test_generate_next_map(self):
        cm=ChainIDManager(transform_reserves={'A':['G','Q'],'B':['H','R']})
        cm.check('A')
        cm.check('B')
        nm=cm.generate_next_map(['A','B'])
        self.assertEqual(nm['A'],'G')
        self.assertEqual(nm['B'],'H')
        nm=cm.generate_next_map(['A','B'])
        self.assertEqual(nm['A'],'Q')
        self.assertEqual(nm['B'],'R')

