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
        r=ResourceManager().ApplyUserOptions(user_options)
        self.assertEqual(r.namd2,user_options['NAMD2'])
        self.assertEqual(r.platform,platform.system())

