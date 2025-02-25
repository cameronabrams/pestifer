from pestifer.util.util import *
import unittest
from pestifer.tasks.psfgen import *
import os
from pestifer.objs.mutation import Mutation
from pestifer import objs

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
        cls,dum=inspect_classes('pestifer.tasks.psfgen',use_yaml_headers_as_keys=True)
        self.assertEqual(cls['psfgen'],PsfgenTask)

    def test_inspect_package_dir(self):
        x,y=inspect_package_dir(os.path.dirname(objs.__file__),key='List')
        self.assertEqual(len(x),15)
        self.assertEqual(x['Mutation'],Mutation)

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
        self.assertFalse(is_periodic('no.xsc'))
        self.assertTrue(is_periodic('yes.xsc'))

    def test_split_list(self):
        mL=[{'a':x,'b':x+1} for x in range(2,12,2)]
        orig_len=len(mL)
        idx=4
        the_id=id(mL[idx])
        dL=split_list(mL,idx)
        # first element in dL is the same object that was idx'th element of mL
        self.assertTrue(id(dL[0])==the_id)
        self.assertTrue(len(mL)==idx)
        self.assertTrue(len(dL)==(orig_len-idx))

        mL=[{'a':x,'b':x+1} for x in range(2,12,2)]
        idx=4
        c1mL=mL[:idx]
        c1dL=mL[idx:]
        # does not copy
        self.assertTrue(mL[idx] is c1dL[0])

    def test_cell_from_xsc(self):
        box,orig=cell_from_xsc('test.xsc')
        self.assertTrue(box.shape==(3,3))
        self.assertTrue(orig.shape==(3,))
        self.assertEqual(box[0][0],191.901209569)
        self.assertEqual(box[0][1],0)
        self.assertEqual(box[0][2],0)
        self.assertEqual(box[1][1],191.901209569)
        self.assertEqual(box[1][0],0)
        self.assertEqual(box[1][2],0)
        self.assertEqual(box[2][2],206.131438627)
        self.assertEqual(box[2][0],0)
        self.assertEqual(box[2][1],0)
        self.assertEqual(orig[0],0)
        self.assertEqual(orig[1],0)
        self.assertEqual(orig[2],0)
