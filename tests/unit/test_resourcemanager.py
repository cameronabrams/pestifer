import unittest
import os
from pestifer.resourcemanager import ResourceManager
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
        c=self.RM.get_example_yaml_by_index(1)
        self.assertEqual(os.path.basename(c),'01-bpti.yaml')

    def test_resource_ycleptic(self):
        self.assertEqual(os.path.basename(self.RM.ycleptic_config),"base.yaml")
        self.assertEqual(os.path.basename(self.RM.get_ycleptic_config()),"base.yaml")

    def test_resource_tcl(self):
        self.assertTrue(os.path.exists(self.RM.get_tcldir()))
        self.assertTrue(os.path.exists(self.RM.get_tcl_pkgdir()))
        self.assertTrue(os.path.exists(self.RM.get_tcl_scriptsdir()))

    def test_pdb_repository(self):
        PD=self.RM.charmmff_content.pdb_repository
        self.assertTrue(PD!=None)
        c=PD.checkout('PSM')
        self.assertTrue(c!=None)
        ch=c.get_charge()
        self.assertTrue(ch==0.0)
        p=c.get_parameters()
        self.assertTrue('toppar_all36_lipid_sphingo.str' in p)
        self.assertTrue(len(c.info['conformers'])==10)
        self.assertEqual(c.info['conformers'][0]['head-tail-length'],26.779)
        self.assertEqual(c.info['conformers'][0]['max-internal-length'],30.398)
        c.get_pdb(0)
        self.assertTrue(os.path.exists('PSM-00.pdb'))
        os.remove('PSM-00.pdb')
        c.get_pdb(0,noh=True)
        self.assertTrue(os.path.exists('PSM-00-noh.pdb'))
        os.remove('PSM-00-noh.pdb')
        c=PD.checkout('TIP3')
        self.assertTrue(c!=None)
        self.assertTrue(c.info=={})
        c.get_pdb(0)
        self.assertTrue(os.path.exists('TIP3.pdb'))
        c=PD.checkout('FAKE')
        self.assertTrue(c==None)
        PD.registercollection('pdb_depot','user')
        c=PD.checkout('FAKE')
        self.assertTrue(c!=None)
        c.get_pdb(0)
        self.assertTrue(os.path.exists('FAKE-00.pdb'))
        os.remove('FAKE-00.pdb')

