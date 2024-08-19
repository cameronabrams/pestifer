# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
# from pestifer.config import ResourceManager
# from pestifer.charmmtop import getResis, makeBondGraph, getMasses, CharmmMasses
from pestifer.PestiferResources.lmsd.make_database import *
# import glob
# import os
import logging
# import gzip
# import networkx as nx
# import yaml
logger=logging.getLogger(__name__)

from rdkit.Chem import ForwardSDMolSupplier
from rdkit import Chem

class TestCHARMMtop(unittest.TestCase):

    def test_resnames(self):
        DB=LMSDDatabase()
        DB.yaml_names()
        
    def test_lipid_charmmtop(self):
        DB=LMSDDatabase()
        result=make_charmm_pdb('DPPC',DB)
        self.assertTrue(result!=-1)

    def test_do_psfgen(self):
        DB=LMSDDatabase()
        result=do_psfgen('DPPC',DB)