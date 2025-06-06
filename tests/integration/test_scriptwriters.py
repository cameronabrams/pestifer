import unittest
from pestifer.scriptwriters import Psfgen, NAMD
from pestifer.config import Config

import os
class TestPsfgen(unittest.TestCase):
    def test_header(self):
        c=Config()
        p=Psfgen(c)
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
        p=NAMD(c)
        assert p.namd_version==3
        c.RM.charmmff_content.clean_local_charmmff_files()
    def test_atomselect_macros(self):
        c=Config()
        p=Psfgen(c)
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
        self.assertEqual(nmacros,4)

    # def test_write_mol(self):
    #     c=ConfigSetup('example.yaml')
    #     chainIDmanager=ChainIDManager()
    #     m=Molecule(source='1gc1').activate_biological_assembly(1,chainIDmanager)
    #     p=Psfgen(c.resman)
    #     p.beginscript()
    #     p.describe_molecule(m,{})
    #     p.endscript()
    #     p.write()
    #     C=Command(f'grep -c "END PESTIFER" {p.script_name}')
    #     rc=C.run()
    #     self.assertEqual(rc,0)
    #     self.assertTrue('1' in C.stdout)
    # def test_clean(self):
    #     c=ConfigSetup('example.yaml')
    #     p=Psfgen(c.resman)
    #     p.write()
    #     p.cleanfiles()
    #     self.assertTrue(not os.path.exists(p.script_name))
