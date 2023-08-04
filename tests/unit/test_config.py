"""

.. module:: test_config
   :synopsis: tests pestifer.Config()
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
import platform
import os
import pytest
from pestifer.config import ConfigSetup, replace, ConfigGetParam
from pestifer.resourcemanager import ResourcesGet

class A:
    d={}
    def __init__(self,idict,**cdict):
        self.d=idict
        self.d.update(**cdict)
        A.d.update(**cdict)

class ClassTest(unittest.TestCase):
    def test_1(self):
        a=A({'a':1,'b':2},c=3,d=4)
        self.assertEqual(A.d['c'],3)
        b=A({},c=6)
        self.assertEqual(A.d['c'],6)

class ConfigTest(unittest.TestCase):
    def setUp(self):
        # print('setUp is here',os.getcwd())
        self.userinputs='example.yaml'
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
        c=ConfigSetup(self.userinputs)
        d=c.defs
        expected_charmmdir=os.path.join(self.user_home,'my_charmm')
        self.assertEqual(d['CHARMMDIR'],expected_charmmdir)

    def test_modreads(self):
        c=ConfigSetup(self.userinputs)
        d=c.defs
        expected_results={
            'Query':'All good!',
            'AnotherQuery':'Yep still good',
            'Mutations':
                [
                {'chainID':'A','resseqnum':123,'newresname':'ARG'},
                'G,555,PHE'
                ]
                }
        test_results={'Query':'NotFound','AnotherQuery':'NotFound'}
        for step in d['BuildSteps']:
            if 'mods' in step:
                test_results.update(step['mods'])
        self.assertDictEqual(test_results,expected_results)
    def test_resids(self):
        c=ConfigSetup(self.userinputs)
        d=c.defs
        self.assertEqual(d['PDB_1char_to_3char_Resnames']['S'],'SER')
        self.assertEqual(d['PDB_to_CHARMM_Resnames']['NAG'],'BGNA')
    def test_seqtypes(self):
        c=ConfigSetup(self.userinputs)
        d=c.defs
        self.assertEqual(d['Segtypes_by_Resnames']['HOH'],'WATER')
        self.assertEqual(d['Segtypes_by_Resnames']['NAG'],'GLYCAN')
        self.assertEqual(d['Segtypes_by_Resnames']['BGNA'],'GLYCAN')
        self.assertEqual(d['Segtypes_by_Resnames']['PRO'],'PROTEIN')
    def test_resources(self):
        c=ConfigSetup(self.userinputs)
        namd2=ResourcesGet('namd2')
        self.assertEqual(namd2,c.defs['NAMD2'])
        plat=ResourcesGet('Platform')
        self.assertEqual(plat,platform.system())

    def test_globality(self):
        c=ConfigSetup(self.userinputs)
        p=ConfigGetParam('BuildSteps')
        self.assertTrue(type(p)==list)
        self.assertEqual(len(p),3)
        p=ConfigGetParam('CHARMM_to_PDB_Resnames')
        self.assertTrue(type(p)==dict)
        self.assertEqual(p['AFUC'],'FUC')
        p=ConfigGetParam('PDB_to_CHARMM_Resnames')
        self.assertTrue(type(p)==dict)
        self.assertEqual(p['FUC'],'AFUC')
