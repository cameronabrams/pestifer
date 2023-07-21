"""

.. module:: test_assets
   :synopsis: tests pestifer.Resources()
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest

from pestifer.resources import ResourceManager

class ResourceManagerTest(unittest.TestCase):
    def setUp(self):
        self.expected_subdirs=['bash','charmm','config','tcl','templates']
        self.expected_tcl_lib='tcl/modules/lib/bondstruct.so'

    def test_libraries(self):
        """test_libraries test to see if installed resources has the expected top-level 
        subdirectories
        """
        l=ResourceManager()
        s1=set(self.expected_subdirs)
        s2=set(l.resource_paths)
        self.assertEqual(s1,s2)

    def test_bondstruct_so(self):
        l=ResourceManager()
        self.assertTrue(any([self.expected_tcl_lib in str(fn) for fn in l.files]))

