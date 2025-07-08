import unittest
from pestifer.core.scripters import PsfgenScripter, NAMDScripter
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

    def test_atomselect_macros(self):
        c=Config()
        p=PsfgenScripter(c)
        if os.path.exists('atomselect_macros.tcl'):
            os.remove('atomselect_macros.tcl')
        p.newfile('atomselect_macros.tcl')
        p.addline('# this is a test -- should write four update_atomselect_macro commands')
        p.atomselect_macros()
        p.writefile()
        with open('atomselect_macros.tcl','r') as f:
            lines=f.read().split('\n')
        nmacros=0
        for l in lines:
            if len(l)==0:
                continue
            w=l.split()
            if w[0]=='update_atomselect_macro':
                nmacros+=1
        self.assertEqual(nmacros,5) # glycan, ligand, ion, lipid, and nucleicacid

