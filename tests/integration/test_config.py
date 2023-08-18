"""

.. module:: test_config
   :synopsis: tests pestifer.Config()
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
import platform
import os
import pytest
from pestifer.config import Config
from pestifer.resourcemanager import ResourceManager
from pestifer.util import replace

class A:
    d={}
    def __init__(self,idict,**cdict):
        self.d=idict
        self.d.update(**cdict)
        A.d.update(**cdict)

class ClassTest(unittest.TestCase):
    def test_1(self):
        a=A({'a':1,'b':2},c=3,d=4)
        self.assertEqual(A.c['c'],3)
        b=A({},c=6)
        self.assertEqual(A.c['c'],6)

class ConfigTest(unittest.TestCase):
    def setUp(self):
        plat=platform.system()
        if plat=='Linux':
            self.user_home=os.environ['HOME']
        elif plat=='Windows':
            self.user_home=os.environ['HOMEPATH']
        elif plat=='Darwin':
            self.user_home=os.environ['HOME']
        else:
            raise(Exception,f'platform {plat} not recognized')
        self.starting_dict={
            'SomeIntegers':[1,2,'$(VARIABLE3)'],
            'SomeWords':{'subdata':['$(VARIABLE2)','hello']},
            'YetMoreData':[0.1,'$(VARIABLE4)'],
            'NestyMcNestface':{'level1':{'level2':{'ALIST':'$(VAR6)'}},'level1-str':'$(VAR5)'},
            'SearchReplace': '$(VAR7)!!!'
        }

        self.expected_dict={
            'SomeIntegers':[1,2,3],
            'SomeWords':{'subdata':['word','hello']},
            'YetMoreData':[0.1,0.99],
            'NestyMcNestface':{'level1':{'level2':{'ALIST':['this','is','a','list']}},'level1-str':'A string'},
            'SearchReplace':'replaced!!!'
        }
    def test_replace(self):
        # print('test_replace here',os.getcwd())
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
            replace(self.starting_dict,s,r)
        self.assertEqual(self.starting_dict,self.expected_dict)
    def test_config(self):
        c=Config(self.userinputs)
        expected_charmmdir=os.path.join(self.user_home,'test_charmm')
        self.assertEqual(c['CHARMMDIR'],expected_charmmdir)

    def test_resids(self):
        c=Config(self.userinputs)
        self.assertEqual(c['PDB_1char_to_3char_Resnames']['S'],'SER')
        self.assertEqual(c['PDB_to_CHARMM_Resnames']['NAG'],'BGNA')
    def test_seqtypes(self):
        c=Config(self.userinputs)
        self.assertEqual(c['Segtypes_by_Resnames']['HOH'],'WATER')
        self.assertEqual(c['Segtypes_by_Resnames']['NAG'],'GLYCAN')
        self.assertEqual(c['Segtypes_by_Resnames']['BGNA'],'GLYCAN')
        self.assertEqual(c['Segtypes_by_Resnames']['PRO'],'PROTEIN')
    def test_resources(self):
        c=Config(self.userinputs)
        r=ResourceManager(c)
        namd2=r.namd2
        self.assertEqual(namd2,c.defs['NAMD2'])
        plat=r.Platform
        self.assertEqual(plat,platform.system())

    def test_dictlikeness(self):
        c=Config(self.userinputs)
        p=c.get('AssHat','Cannot find an asshat')
        self.AssertEqual(p,'Cannot find an asshat')
        p=c.get('steps')
        self.assertTrue(type(p)==list)
        self.assertEqual(len(p),3)
