import unittest
from pestifer.scripters import PsfgenScripter, NAMDScripter
from pestifer.core.config import Config

import os
class TestPsfgen(unittest.TestCase):

    def test_header(self):
        c=Config()
        p=PsfgenScripter(c)
        p.newscript()
        p.B.banner('TESTING')
        p.writescript('test-state')
        self.assertTrue(os.path.isfile('pestifer-script.tcl'))
        p.newscript('testing')
        p.B.banner('TESTING')
        p.writescript('test-state')
        self.assertTrue(os.path.isfile('testing.tcl'))
        os.remove(p.basename+'.tcl')
        self.assertFalse(os.path.exists('testing.tcl'))
        c.RM.charmmff_content.clean_local_charmmff_files()

    def test_charmm(self):
        c=Config()
        p=NAMDScripter(c)
        assert p.namd_version==3
        c.RM.charmmff_content.clean_local_charmmff_files()


