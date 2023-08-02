"""

.. module:: test_assets
   :synopsis: tests pestifer.Resources()
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest

from pestifer.resourcemanager import ResourceManager as RM
from pestifer.resourcemanager import excludes
from pestifer import Resources
import os
import glob

class ResourceManagerTest(unittest.TestCase):
    def test_RM2(self):
        self.assertEqual(RM.ResourcesRoot,os.path.dirname(Resources.__file__))
    def test_RM3(self):
        m=glob.glob(os.path.dirname(Resources.__file__)+'/*')
        self.assertEqual(RM.ResourceFullPaths,m)
    def test_RM4(self):
        config_path=[x for x in RM.ResourceFullPaths if 'config' in x][0]
        self.assertEqual(RM.ResourcePaths['config'],config_path)
    def test_RM5(self):
        self.assertEqual(len(RM.ResourceDirs),2)
    def test_RM5(self):
        self.assertTrue('platform' in RM.System)
        self.assertTrue('user_home' in RM.System)
    # def setUp(self):
    #     self.expected_subdirs=['bash','charmm','config','tcl','templates']
    #     self.expected_tcl_lib='tcl/modules/lib/bondstruct.so'

    # def test_libraries(self):
    #     """test_libraries test to see if installed resources has the expected top-level 
    #     subdirectories
    #     """
    #     s1=set(self.expected_subdirs)
    #     s2=set(RM.resource_paths)
    #     self.assertEqual(s1,s2)

    # def test_bondstruct_so(self):
    #     self.assertTrue(any([self.expected_tcl_lib in str(fn) for fn in RM.files]))

