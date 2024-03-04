# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
import os
from pestifer.condacheck import CondaCheck
from pestifer.config import Config
import numpy
class TestCondaCheck(unittest.TestCase):
    def test_condacheck(self):
        home=os.environ['HOME']
        has_conda=False
        if os.path.exists(os.path.join(home,'.condarc')):
            has_conda=True
        c=CondaCheck()
        if not has_conda:
            self.assertFalse(c.has_conda)
        else:
            self.assertTrue(c.has_conda!=None)
            self.assertTrue(len(c.conda_envs)>0)
    
    def test_getpackageversion(self):
        c=CondaCheck()
        real_numpy_version=str(numpy.__version__)
        test_numpy_version=c.get_package_version(c.current_conda_env,'numpy')
        self.assertEqual(real_numpy_version,test_numpy_version)

class TestPackmolMemgenAvailable(unittest.TestCase):
    def test_packmolmemgenavailable(self):
        c=Config()
        # this is specific to MY installation
        self.assertEqual(c['user']['packmol_memgen']['available'],True)
        self.assertEqual(c['user']['packmol_memgen']['local'],False)
        self.assertEqual(c['user']['packmol_memgen']['venv'],'pmmg')