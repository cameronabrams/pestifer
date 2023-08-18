"""

.. module:: test_assets
   :synopsis: tests pestifer.Resources()
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
from pestifer.resourcemanager import ResourceManager
import platform

class ResourceManagerTest(unittest.TestCase):
    def test_RM2(self):
        user_options={
            'CHARMMDIR':'/home/cfa/charmm',
            'CHARMRUN':'/usr/local/bin/charmrun',
            'NAMD2':'/usr/local/bin/namd2',
            'VMD':'/usr/local/bin/vmd',
        }
        r=ResourceManager(user_options)
        self.assertEqual(r.namd2,user_options['NAMD2'])
        self.assertEqual(r.plat,platform.system())
        # print(ResourcesInfo())
        # self.assertTrue(False)
    # def test_libraries(self):
    #     """test_libraries test to see if installed resources has the expected top-level 
    #     subdirectories
    #     """
    #     s1=set(self.expected_subdirs)
    #     s2=set(RM.resource_paths)
    #     self.assertEqual(s1,s2)

    # def test_bondstruct_so(self):
    #     self.assertTrue(any([self.expected_tcl_lib in str(fn) for fn in RM.files]))

