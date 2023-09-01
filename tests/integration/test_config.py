"""

.. module:: test_config
   :synopsis: tests pestifer.Config()
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
import platform
import os
import yaml
from pestifer.config import Config

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
        self.userinputs='user_config.yaml'

    def test_config(self):
        c=Config(self.userinputs)
        expected_charmmdir=os.path.join(self.user_home,'charmm')
        self.assertEqual(c.user_charmm_path,expected_charmmdir)
        expected_charmmtoppardir=os.path.join(self.user_home,'charmm','toppar')
        self.assertEqual(c.user_charmm_toppar_path,expected_charmmtoppardir)

    def test_resids(self):
        c=Config(self.userinputs)
        self.assertEqual(c['PDB_1char_to_3char_Resnames']['S'],'SER')
        self.assertEqual(c['PDB_to_CHARMM_Resnames']['NAG'],'BGLCNA')
    def test_seqtypes(self):
        c=Config(self.userinputs)
        self.assertEqual(c['Segtypes_by_Resnames']['HOH'],'WATER')
        self.assertEqual(c['Segtypes_by_Resnames']['NAG'],'GLYCAN')
        self.assertEqual(c['Segtypes_by_Resnames']['BGLCNA'],'GLYCAN')
        self.assertEqual(c['Segtypes_by_Resnames']['ALA'],'PROTEIN')
        self.assertEqual(c.segtype('PHE'),'PROTEIN')

    def test_dictlikeness(self):
        c=Config(self.userinputs)
        p=c.get('AssHat','Cannot find an asshat')
        self.assertEqual(p,'Cannot find an asshat')
        p=c.get('steps',[])
        with open(self.userinputs,'r') as f:
            d=yaml.safe_load(f)
        self.assertTrue(type(p)==list)
        self.assertEqual(p,d['steps'])
        self.assertEqual(len(p),2)

    def test_boolean(self):
        c=Config()
        self.assertTrue(c)