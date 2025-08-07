import os
import unittest

from pestifer.core.config import Config
from pestifer.scripters import PsfgenScripter, NAMDScripter

class TestPsfgen(unittest.TestCase):

    def test_header(self):
        c=Config()
        p=c.get_scripter('psfgen')
        self.assertIsInstance(p, PsfgenScripter)
        p.newscript()
        p.B.banner('TESTING')
        p.writescript('test-state')
        self.assertTrue(os.path.isfile('pestifer-script.tcl'))
        os.remove(p.basename+'.tcl')
        p.newscript('testing')
        p.B.banner('TESTING')
        p.writescript('test-state')
        self.assertTrue(os.path.isfile('testing.tcl'))
        os.remove(p.basename+'.tcl')
        self.assertFalse(os.path.exists('testing.tcl'))
        c.RM.charmmff_content.clean_local_charmmff_files()

    def test_charmm(self):
        c=Config()
        p=c.get_scripter('namd')
        self.assertIsInstance(p, NAMDScripter)
        self.assertEqual(p.namd_version, 3)
        c.RM.charmmff_content.clean_local_charmmff_files()


