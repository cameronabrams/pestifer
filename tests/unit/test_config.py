"""

.. module:: test_config
   :synopsis: tests pestifer.Config()
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
import importlib.resources
import platform
import os
from pestifer.config import Config, replace
from pestifer.resources import ResourceManager

class ConfigTest(unittest.TestCase):
    def setUp(self):
        self.userinputs=str(importlib.resources.files('tests.unit').joinpath(f'fixtures/example.yaml'))
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
        for s,r in self.replacements.items():
            replace(self.starting_dict,s,r)
        self.assertEqual(self.starting_dict,self.expected_dict)
    def test_config(self):
        r=ResourceManager()
        c=Config(self.userinputs,r)
        d=c.data
        expected_charmmdir=os.path.join(self.user_home,'my_charmm')
        self.assertEqual(d['CHARMMDIR'],expected_charmmdir)
        