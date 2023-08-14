import unittest
from pestifer.command import Command
from pestifer.molecule import Molecule
from pestifer.psfgen import Psfgen
from pestifer.config import ConfigSetup, ConfigGetParam
from pestifer.chainids import ChainIDManager

import os
class TestPsfgen(unittest.TestCase):
    def test_header(self):
        c=ConfigSetup('example.yaml')
        n=ConfigGetParam('psfgen_script_name')
        p=Psfgen(c.resman)
        self.assertEqual(p.script_name,n)
        p.write()
        self.assertTrue(os.path.isfile(p.script_name))
    def test_write_mol(self):
        c=ConfigSetup('example.yaml')
        chainIDmanager=ChainIDManager()
        m=Molecule(source='1gc1').activate_biological_assembly(1,chainIDmanager)
        p=Psfgen(c.resman)
        p.beginscript()
        p.describe_molecule(m,{})
        p.endscript()
        p.write()
        C=Command(f'grep -c "END PESTIFER" {p.script_name}')
        rc=C.run()
        self.assertEqual(rc,0)
        self.assertTrue('1' in C.stdout)
    def test_clean(self):
        c=ConfigSetup('example.yaml')
        p=Psfgen(c.resman)
        p.write()
        p.cleanfiles()
        self.assertTrue(not os.path.exists(p.script_name))
