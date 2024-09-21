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
        topo=DB.get_topo(resid)
        self.assertTrue(topo.metadata['substream']=='cholesterol')
        topo.lipid_annotate()
        heads=topo.annotation['heads']
        tails=topo.annotation['tails']
        self.assertTrue(heads==['C3'])
        self.assertTrue(tails==['C27'])
    
    def test_detect_structure_CHM1(self):
        DB=CharmmResiDatabase()
        DB.add_stream('lipid')
        resid='CHM1'
        topo=DB.get_topo(resid)
        self.assertTrue(topo.metadata['substream']=='cholesterol')
        topo.lipid_annotate()
        heads=topo.annotation['heads']
        tails=topo.annotation['tails']
        self.assertTrue(heads==['C1'])
        self.assertTrue(tails==['C6'])

    def test_detect_structure_DPPC(self):
        DB=CharmmResiDatabase()
        DB.add_stream('lipid')
        resid='DPPC'
        topo=DB.get_topo(resid)
        self.assertTrue(topo.metadata['substream']=='')
        topo.lipid_annotate()
        heads=topo.annotation['heads']
        tails=topo.annotation['tails']
        shortest_paths=topo.annotation['shortest_paths']
        logger.debug(f'{heads} {tails} {shortest_paths}')
        dist1=shortest_paths[heads[0]][tails[0]]
        dist2=shortest_paths[heads[0]][tails[1]]
        # logger.debug(f'dist1 {dist1} dist2 {dist2}')
        self.assertEqual(dist1,25)
        self.assertEqual(dist2,26)

    def test_detect_structure_SDS(self):
        DB=CharmmResiDatabase()
        DB.add_stream('lipid')
        resid='SDS'
        topo=DB.get_topo(resid)
        self.assertTrue(topo.metadata['substream']=='detergent')
        topo.lipid_annotate()
        heads=topo.annotation['heads']
        tails=topo.annotation['tails']
        self.assertTrue(heads==['S'])
        self.assertTrue(tails==['C12'])

    def test_detect_structure_TOCL1(self):
        DB=CharmmResiDatabase()
        DB.add_stream('lipid')
        resid='TOCL1'
        topo=DB.get_topo(resid)
        self.assertTrue(topo.metadata['substream']=='cardiolipin')
        topo.lipid_annotate()
        heads=topo.annotation['heads']
        tails=topo.annotation['tails']
        self.assertEqual(len(tails),4)
        self.assertTrue(heads==['C2'])
        self.assertTrue(tails==['CA18', 'CB18', 'CC18', 'CD18'])