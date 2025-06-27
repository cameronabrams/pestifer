# Author: Cameron F. Abrams

import unittest
import os

from pestifer.core.pestifer import mdplot  # this is the subcommand
from argparse import Namespace
class Test_MDPlot(unittest.TestCase):
    def test_eplot(self):
        inputs=os.listdir('inputs')
        logs=[os.path.join('inputs',f) for f in inputs if f.endswith('.lxg')]
        logs.sort()
        xsts=[os.path.join('inputs',f) for f in inputs if f.endswith('.xst')]
        xsts.sort()
        args=Namespace(logs=logs,xsts=xsts,
                       basename='testmdplot',profiles=[],
                       figsize=[9,6],traces=['density',['a_x','b_y','c_z']])
        mdplot(args)

        self.assertTrue(os.path.exists('testmdplot-energy.csv'))
        self.assertTrue(os.path.exists('testmdplot-cell.csv'))
        self.assertTrue(os.path.exists('testmdplot-density.png'))
        self.assertTrue(os.path.exists('testmdplot-a_x-b_y-c_z.png'))
        
        # os.remove('states.yaml')
        # os.remove('testmdplot-density.png')
        # os.remove('testsave.csv')
        # os.remove('xst-testsave.csv')
        # os.remove('testmdplot-a_x-b_y-c_z.png')



