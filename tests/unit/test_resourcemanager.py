import unittest
import os
from pestifer.core.resourcemanager import ResourceManager
from pestifer import resources

class TestResourceManager(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Ensure the resource manager is initialized before any tests run
        cls.RM = ResourceManager()

    def test_resource_manager(self):
        self.assertEqual(self.RM.resources_path,os.path.dirname(resources.__file__))
        for r in self.RM.base_resources:
            self.assertTrue(os.path.exists(self.RM.resource_path[r]))
    
    def test_resource_get_example_yaml(self):
        c=self.RM.example_manager.checkout_example_yaml(1)
        self.assertEqual(os.path.basename(c),'bpti1.yaml')
        self.assertTrue(os.path.exists(c))
        os.remove(c)

    def test_resource_ycleptic(self):
        self.assertEqual(os.path.basename(self.RM.ycleptic_config),"base.yaml")
        self.assertEqual(os.path.basename(self.RM.get_ycleptic_config()),"base.yaml")

    def test_resource_tcl(self):
        self.assertTrue(os.path.exists(self.RM.get_tcldir()))
        self.assertTrue(os.path.exists(self.RM.get_tcl_pkgdir()))
        self.assertTrue(os.path.exists(self.RM.get_tcl_scriptsdir()))


