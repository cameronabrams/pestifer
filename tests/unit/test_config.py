"""

.. module:: test_config
   :synopsis: tests pestifer.Config()
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
import platform
import os
from pestifer.config import Config, replace
from pestifer.resources import ResourceManager

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
        self.replacements={
            'VARIABLE1':'SomeIntegers',
            'VARIABLE3':3,
            'VARIABLE2':'word',
            'VARIABLE4':0.99,
            'VAR5':'A string',
            'VAR6': ['this','is','a','list'],
            'VAR7': 'replaced'
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
        for s,r in self.replacements.items():
            replace(self.starting_dict,s,r)
        self.assertEqual(self.starting_dict,self.expected_dict)
    def test_config(self):
        r=ResourceManager()
        c=Config(self.userinputs,r)
        d=c.data
        expected_charmmdir=os.path.join(self.user_home,'my_charmm')
        self.assertEqual(d['CHARMMDIR'],expected_charmmdir)
    def test_modreads(self):
        r=ResourceManager()
        c=Config(self.userinputs,r)
        d=c.data
        expected_results={'Query':'All good!','AnotherQuery':'Yep still good'}
        test_results={'Query':'NotFound','AnotherQuery':'NotFound'}
        for step in d['BuildSteps']:
            if 'mods' in step:
                for mod in step['mods']:
                    test_results.update(mod)
        self.assertDictEqual(test_results,expected_results)
    def test_resids(self):
        r=ResourceManager()
        c=Config(self.userinputs,r)
        d=c.data
        self.assertEqual(d['PDB_1char_to_3char_Resnames']['S'],'SER')
        self.assertEqual(d['PDB_to_CHARMM_Resnames']['NAG'],'BGNA')
    def test_seqtypes(self):
        r=ResourceManager()
        c=Config(self.userinputs,r)
        d=c.data
        self.assertEqual(d['Segtypes_by_Resnames']['HOH'],'WATER')
        self.assertEqual(d['Segtypes_by_Resnames']['NAG'],'GLYCAN')
        self.assertEqual(d['Segtypes_by_Resnames']['BGNA'],'GLYCAN')
        self.assertEqual(d['Segtypes_by_Resnames']['PRO'],'PROTEIN')
