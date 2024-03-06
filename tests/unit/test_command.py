# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import unittest
from pestifer.command import Command,CondaCheck
from pestifer.config import Config
import numpy
import os

class TestCommand(unittest.TestCase):
    def test_echo(self):
        expected_stdout='123\n'
        c=Command('echo 123')
        c.run()
        self.assertEqual(expected_stdout,c.stdout)

    def test_condafy(self):
        pass

class TestCondaCheck(unittest.TestCase):
    def test_condacheck(self):
        home=os.environ['HOME']
        has_conda=False
        if os.path.exists(os.path.join(home,'.condarc')):
            has_conda=True
        c=CondaCheck()
        if not has_conda:
            self.assertFalse(c.conda_root)
        else:
            self.assertTrue(c.conda_root!=None)
            self.assertTrue(len(c.conda_envs)>0)
    
    def test_getpackageversion(self):
        c=CondaCheck()
        real_numpy_version=str(numpy.__version__)
        test_numpy_version=c.get_package_version(c.active_env,'numpy')
        self.assertEqual(real_numpy_version,test_numpy_version)
        test_ambertools_version=c.get_package_version(c.active_env,'ambertools',from_list=True)
        self.assertEqual(test_ambertools_version,'23.3')

class TestAmbertools(unittest.TestCase):
    def test_ambertools_available_by_env_search(self):
        c=Config('noenv.yaml')
        # this is specific to MY installation
        self.assertEqual(c['user']['ambertools']['available'],True)
        self.assertEqual(c['user']['ambertools']['local'],False)
        self.assertEqual(c['user']['ambertools']['venv'],'pmmg')
    
    def test_ambertools_available_by_explicit_env(self):
        c=Config('env.yaml')
        # this is specific to MY installation
        self.assertEqual(c['user']['ambertools']['available'],True)
        self.assertEqual(c['user']['ambertools']['local'],False)
        self.assertEqual(c['user']['ambertools']['venv'],'pmmg')

    def test_no_ambertools(self):
        c=Config()
        self.assertFalse('available' in c['user']['ambertools'])
        self.assertTrue(c['user']['ambertools']['is_unnecessary'])