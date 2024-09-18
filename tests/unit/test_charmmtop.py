# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
# from pestifer.config import ResourceManager
# from pestifer.charmmtop import getResis, makeBondGraph, getMasses, CharmmMasses
from pestifer.charmmtop import *
# import glob
# import os
import logging
import shutil
# import gzip
# import networkx as nx
# import yaml
logger=logging.getLogger(__name__)

class TestCHARMMtop(unittest.TestCase):

    def test_detect_structure_CHL1(self):
        DB=CharmmResiDatabase()
        DB.add_stream('lipid')
        resid='CHL1'
        charmm_topfile=DB.get_charmm_topfile(resid)
        topo=DB.get_topo(resid)
        heads,tails,shortest_paths=topo.head_tail_atoms()
        logger.debug(f'{heads} {tails} {shortest_paths}')
        self.assertTrue(heads==['O3'])
        self.assertTrue(tails==['C27'])

    def test_detect_structure_DPPC(self):
        DB=CharmmResiDatabase()
        DB.add_stream('lipid')
        resid='DPPC'
        charmm_topfile=DB.get_charmm_topfile(resid)
        topo=DB.get_topo(resid)
        heads,tails,shortest_paths=topo.head_tail_atoms()
        logger.debug(f'{heads} {tails} {shortest_paths}')
        dist1=shortest_paths[heads[0]][tails[0]]
        dist2=shortest_paths[heads[0]][tails[1]]
        # logger.debug(f'dist1 {dist1} dist2 {dist2}')
        self.assertEqual(dist1,25)
        self.assertEqual(dist2,26)