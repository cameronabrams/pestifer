import unittest
import os
from pestifer.core.config import Config
from pestifer.scripters import NAMDScripter

class TestNAMDScripter(unittest.TestCase):

    def test_charmm(self):
        c = Config().configure_new()
        p: NAMDScripter = c.get_scripter('namd')
        self.assertIsInstance(p, NAMDScripter)
        self.assertEqual(p.namd_version, 3)
        c.RM.charmmff_content.clean_local_charmmff_files()
