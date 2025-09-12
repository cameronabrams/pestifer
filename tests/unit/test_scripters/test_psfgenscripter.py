import os
import unittest

from pestifer.core.config import Config
from pestifer.scripters import PsfgenScripter

class TestPsfgenScripter(unittest.TestCase):

    def setUp(self):
        self.c = Config().configure_new()
        self.p: PsfgenScripter = self.c.get_scripter('psfgen')

    def tearDown(self):
        self.c.RM.charmmff_content.clean_local_charmmff_files()

    def test_header(self):
        self.assertIsInstance(self.p, PsfgenScripter)
        self.p.newscript()
        self.p.B.banner('TESTING')
        self.p.writescript('test-state')
        self.assertTrue(os.path.isfile('pestifer-script.tcl'))
        os.remove(self.p.basename+'.tcl')
        self.p.newscript('testing')
        self.p.B.banner('TESTING')
        self.p.writescript('test-state')
        self.assertTrue(os.path.isfile('testing.tcl'))
        os.remove(self.p.basename+'.tcl')
        self.assertFalse(os.path.exists('testing.tcl'))
    
    def test_fetch_charmff_files(self):
        self.p.fetch_standard_charmm_topologies()
        self.assertTrue(os.path.isfile('top_all36_prot.rtf'))
        self.assertTrue(os.path.isfile('toppar_water_ions.str'))

