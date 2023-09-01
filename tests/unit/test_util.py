from pestifer.util import *
import unittest
from pestifer.tasks import *
import os

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
        l=[1,2]
        L=reduce_intlist(l)
        self.assertEqual(L,'1 2')
    def test_inspect_classes(self):
        cls=inspect_classes('pestifer.tasks',use_yaml_headers_as_keys=True)
        self.assertTrue(cls['psfgen'],PsfgenTask)

    def test_replace(self):
        starting_dict={
            'SomeIntegers':[1,2,'$(VARIABLE3)'],
            'SomeWords':{'subdata':['$(VARIABLE2)','hello']},
            'YetMoreData':[0.1,'$(VARIABLE4)'],
            'NestyMcNestface':{'level1':{'level2':{'ALIST':'$(VAR6)'}},'level1-str':'$(VAR5)'},
            'SearchReplace': '$(VAR7)!!!'
        }

        expected_dict={
            'SomeIntegers':[1,2,3],
            'SomeWords':{'subdata':['word','hello']},
            'YetMoreData':[0.1,0.99],
            'NestyMcNestface':{'level1':{'level2':{'ALIST':['this','is','a','list']}},'level1-str':'A string'},
            'SearchReplace':'replaced!!!'
        }

        replacements={
            'VARIABLE1':'SomeIntegers',
            'VARIABLE3':3,
            'VARIABLE2':'word',
            'VARIABLE4':0.99,
            'VAR5':'A string',
            'VAR6': ['this','is','a','list'],
            'VAR7': 'replaced'
        }
        for s,r in replacements.items():
            replace(starting_dict,s,r)
        self.assertEqual(starting_dict,expected_dict)
    
    def test_is_periodic(self):
        self.assertFalse(is_periodic(None,'no.xsc'))
        self.assertTrue(is_periodic(None,'yes.xsc'))
        self.assertTrue(is_periodic('yes_cell.tcl',None))

    def test_get_version(self):
        self.assertEqual(get_version(),'1.0.1')