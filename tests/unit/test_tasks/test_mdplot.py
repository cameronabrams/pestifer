# Author: Cameron F. Abrams

import unittest
import os

from pestifer.subcommands.mdplot_subcommand import MDPlotSubcommand
from argparse import Namespace

class Test_MDPlot(unittest.TestCase):
    def test_eplot(self):
        inputs = os.listdir('inputs')
        logs = [os.path.join('inputs', f) for f in inputs if f.endswith('.lxg')]
        logs.sort()
        xsts = [os.path.join('inputs', f) for f in inputs if f.endswith('.xst')]
        xsts.sort()
        args = Namespace(logs=logs, xsts=xsts,
                         basename='testmdplot', profiles=[], profiles_per_block=1,
                         figsize=[9, 6], timeseries=['density'], timecoseries=['a_x', 'b_y', 'c_z'],
                         colormap='viridis', colormap_direction=-1)
        file_artifacts = MDPlotSubcommand.func(args)
        expected_artifact_count = 26
        self.assertEqual(len(file_artifacts), expected_artifact_count)

        for artifact in file_artifacts:
            self.assertTrue(artifact.exists())
            os.remove(artifact.path)
        