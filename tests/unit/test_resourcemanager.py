import unittest
import glob
import os
import yaml
from pestifer.resourcemanager import ResourceManager
from pestifer.pdbcollection import PDBInput
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

    def test_resource_tcl(self):
        RM=ResourceManager()
        self.assertTrue(os.path.exists(RM.get_tcldir()))
        self.assertTrue(os.path.exists(RM.get_tcl_pkgdir()))
        self.assertTrue(os.path.exists(RM.get_tcl_scriptsdir()))

    def test_resource_pdb_path(self):
        RM=ResourceManager()
        t2=RM.get_pdb('PSM')
        self.assertTrue(t2!=None)
        self.assertEqual(type(t2),PDBInput)
        self.assertTrue(os.path.exists(t2.conformers[0]))
        self.assertTrue(len(t2.info)>0)
        self.assertTrue('parameters' in t2.info)
        self.assertTrue('conformers' in t2.info)
        self.assertTrue('charge' in t2.info)
        self.assertTrue('reference-atoms' in t2.info)
        t2=RM.get_pdb('PXM')
        self.assertTrue(t2 == None)
        t3=RM.get_pdb('FAKE')
        self.assertTrue(t3==None)
        RM.pdb_collection.registercollection('pdb_depot','user')
        self.assertTrue('user' in RM.pdb_collection.collections)
        t3=RM.get_pdb('FAKE')
        self.assertEqual(type(t3),PDBInput)
        self.assertTrue(os.path.exists(t3.conformers[0]))
