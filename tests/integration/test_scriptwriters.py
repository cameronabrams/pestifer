import unittest
from pestifer.command import Command
from pestifer.molecule import Molecule
from pestifer.scriptwriters import Psfgen, NAMD2
from pestifer.config import Config
from pestifer.chainidmanager import ChainIDManager

import os
class TestPsfgen(unittest.TestCase):
    def test_header(self):
        c=Config()
        p=Psfgen(c)
        p.newscript()
        p.B.banner('TESTING')
        p.writescript()
        self.assertTrue(os.path.isfile('pestifer-script.tcl'))
        p.newscript('testing')
        p.B.banner('TESTING')
        p.writescript()
        self.assertTrue(os.path.isfile('testing.tcl'))
        os.remove(p.basename+'.tcl')
        self.assertFalse(os.path.exists('testing.tcl'))
    def test_charmm(self):
        c=Config()
        p=NAMD2(c)
        print(p.standard_charmmff_parfiles)
        for x in p.standard_charmmff_parfiles:
            self.assertTrue('pestifer/PestiferResources' in x)
    def test_atomselect_macros(self):
        c=Config()
        p=Psfgen(c)
        if os.path.exists('atomselect_macros.tcl'):
            os.remove('atomselect_macros.tcl')
        p.newfile('atomselect_macros.tcl')
        p.addline('# this is a test -- should write atomselect macros')
        p.atomselect_macros()
        p.writefile()
        with open('atomselect_macros.tcl','r') as f:
            lines=f.read().split('\n')
        nmacros=0
        for l in lines:
            if len(l)==0:
                continue
            w=l.split()
            if w[0]=='atomselect' and w[1]=='macro':
                nmacros+=1
        self.assertEqual(nmacros,3)

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
