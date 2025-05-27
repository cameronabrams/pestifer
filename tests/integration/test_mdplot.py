# Author: Cameron F. Abrams

import unittest
import glob
import os

from pestifer.pestifer import mdplot  # this is the subcommand
from argparse import Namespace
class Test_MDPlot(unittest.TestCase):
    def test_eplot(self):
        logs=glob.glob('*.log')
        logs.sort()
        xsts=glob.glob('*.xst')
        xsts.sort()
        args=Namespace(logs=logs,xsts=xsts,
                       basename='testmdplot',
                       savedata='testsave.csv',figsize=[9,6],traces=['density',['a_x','b_y','c_z']])
        mdplot(args)

        self.assertTrue(os.path.exists('testsave.csv'))
        self.assertTrue(os.path.exists('xst-testsave.csv'))
        self.assertTrue(os.path.exists('testmdplot-density.png'))
        self.assertTrue(os.path.exists('testmdplot-a_x-b_y-c_z.png'))
        
        # os.remove('states.yaml')
        # os.remove('testmdplot-density.png')
        # os.remove('testsave.csv')
        # os.remove('xst-testsave.csv')
        # os.remove('testmdplot-a_x-b_y-c_z.png')



