import unittest
# from pestifer.config import ResourceManager
# from pestifer.charmmtop import getResis, makeBondGraph, getMasses, CharmmMasses
from pestifer.colvars import *
from pestifer.scriptwriters import *
# import glob
# import os
import logging
# import gzip
# import networkx as nx
# import yaml
logger=logging.getLogger(__name__)

class TestColvars(unittest.TestCase):
    def test_write_colvars(self):
        pdb='DPPC.pdb'
        colvar_specs={
            'groups': {
                'head1': {'atomnames': ['C2']},
                'head2': {'atomnames': ['N']},
                'tail1': {'atomnames': ['C216']},
                'tail2': {'atomnames': ['C316']}
                },
            'distances': {
                'tail1_tail2': {'groups': ['tail1','tail2']},
                'head1_tail1': {'groups': ['head1','tail1']},
                'head1_tail2':{'groups': ['head1','tail2']}
                },
            'harmonics': {
                'tail1_tail2_attract': {
                    'colvars': ['tail1_tail2'],
                    'forceConstant': [1.0],
                    'distance':[4.0]},
                'head1_tail12_repulse': {
                    'colvars': ['head1_tail1','head1_tail2'],
                    'forceConstant': [-1.0, -1.0],
                    'distance':[20.0,20.0]}
            }
        }
        writer=Filewriter()
        writer.newfile(f'cv.inp')
        colvar_writer(colvar_specs,writer,pdb=pdb)
        writer.writefile()