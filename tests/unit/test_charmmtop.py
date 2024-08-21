# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
# from pestifer.config import ResourceManager
# from pestifer.charmmtop import getResis, makeBondGraph, getMasses, CharmmMasses
from pestifer.PestiferResources.charmmff.pdb.make_database import *
# import glob
# import os
import logging
import shutil
# import gzip
# import networkx as nx
# import yaml
logger=logging.getLogger(__name__)

class TestCHARMMtop(unittest.TestCase):

    def test_do_psfgen(self):
        DB=CharmmResiDatabase(streams=['lipid'])
        lname='CHL1'
        if not os.path.exists(lname):
            os.mkdir(lname)
        else:
            shutil.rmtree(lname)
            os.mkdir(lname)
        os.chdir(lname)
        result=do_psfgen(lname,DB)