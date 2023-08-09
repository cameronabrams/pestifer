from pestifer.util import *
import unittest

class TestUtil(unittest.TestCase):
    def test_special_update(self):
        d1={'A_LIST':[1,2,3],'B_DICT':{'k1':'v1','k2':'v2'},'C_scal':11}
        d2={'A_LIST':[4,5,6],'B_DICT':{'k3':'v3'},'C_scal':12,'D_new':'Hello'}
        d1=special_update(d1,d2)
        self.assertEqual(d1['A_LIST'],[1,2,3,4,5,6])
        self.assertEqual(d1['B_DICT'],{'k1':'v1','k2':'v2','k3':'v3'})
        self.assertEqual(d1['C_scal'],12)
        self.assertEqual(d1['D_new'],'Hello')
    def test_reduce_intlist(self):
        l=[1,2,3,4,5,7,8,9,10,12]
        L=reduce_intlist(l)
        self.assertEqual(L,'1 to 5 7 to 10 12')
        l=list(range(100))
        l.remove(23)
        l.remove(73)
        l.remove(74)
        l.remove(75)
        L=reduce_intlist(l)
        self.assertEqual(L,'0 to 22 24 to 72 76 to 99')
