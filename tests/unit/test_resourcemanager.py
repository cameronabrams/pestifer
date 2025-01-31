import unittest
import glob
import os
import yaml
from pestifer.resourcemanager import ResourceManager
from pestifer import resources

class TestResourceManager(unittest.TestCase):
    def test_resource_manager(self):
        RM=ResourceManager()
        self.assertEqual(RM.resources_path,os.path.dirname(resources.__file__))
        for r in RM.base_resources:
            self.assertTrue(os.path.exists(RM.resource_path[r]))
    
    def test_resource_get_example_yaml(self):
        RM=ResourceManager()
        c=RM.get_example_yaml_by_index(1)
        self.assertEqual(os.path.basename(c),'01-bpti.yaml')

    def test_resource_ycleptic(self):
        RM=ResourceManager()
        self.assertEqual(os.path.basename(RM.ycleptic_config),"base.yaml")
        self.assertEqual(os.path.basename(RM.get_ycleptic_config()),"base.yaml")

    def test_resource_charmmff(self):
        RM=ResourceManager()
        self.assertTrue(os.path.exists(RM.get_charmmff_customdir()))
        self.assertTrue(os.path.exists(RM.get_charmmff_toppardir()))
        self.assertTrue(os.path.exists(RM.get_charmmff_pdbdir()))

    def test_resource_tcl(self):
        RM=ResourceManager()
        self.assertTrue(os.path.exists(RM.get_tcldir()))
        self.assertTrue(os.path.exists(RM.get_tcl_pkgdir()))
        self.assertTrue(os.path.exists(RM.get_tcl_scriptsdir()))

    def test_resource_pdb_path(self):
        RM=ResourceManager()
        t1=RM.get_charmmff_pdb_path('PSM')
        self.assertTrue(t1 != None)
        t2=RM.get_charmmff_pdb_path('PXM')
        self.assertTrue(t2 == None)
