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
        for r in self.RM._base_resources:
            self.assertTrue(os.path.exists(self.RM.resource_path[r]))
    
    def test_resource_get_example_yaml(self):
        e1 = self.RM.example_manager.checkout_example(1)
        self.assertEqual(os.path.basename(e1.scriptname),'bpti1.yaml')
        self.assertTrue(os.path.exists(e1.scriptname))
        os.remove(e1.scriptname)

    def test_resource_tcl(self):
        self.assertTrue(os.path.exists(self.RM.get_tcldir()))
        self.assertTrue(os.path.exists(self.RM.get_tcl_pkgdir()))
        self.assertTrue(os.path.exists(self.RM.get_tcl_scriptsdir()))

    def test_lookup_resname(self):
        # a lipid: defined in the topologies and present in the PDB repository
        popc = self.RM.lookup_resname('popc')  # case-insensitive
        self.assertEqual(popc['resname'], 'POPC')
        self.assertTrue(popc['in_topology'])
        self.assertEqual(popc['kind'], 'residue (RESI)')
        self.assertTrue(popc['topfile'])
        self.assertEqual(popc['segtype'], 'lipid')
        self.assertTrue(popc['in_pdbrepository'])
        self.assertGreaterEqual(popc['nconformers'], 1)
        # a patch: defined in the topologies (PRES) but has no coordinates
        aspp = self.RM.lookup_resname('ASPP')
        self.assertTrue(aspp['in_topology'])
        self.assertEqual(aspp['kind'], 'patch (PRES)')
        self.assertTrue(aspp['is_patch'])          # a pure patch...
        self.assertFalse(aspp['in_pdbrepository'])  # ...so the PDB repo is not searched
        # a name that exists nowhere
        bogus = self.RM.lookup_resname('ZZZZ')
        self.assertFalse(bogus['in_topology'])
        self.assertIsNone(bogus['topfile'])
        self.assertFalse(bogus['in_pdbrepository'])
        # richer detail for a PDB-repo lipid: charge present, head-to-tail length > 0
        self.assertIsNotNone(popc['charge'])
        self.assertGreater(popc['head_tail_length'], 0.0)

    def test_search_resnames(self):
        allnames = self.RM.all_resnames()
        self.assertIn('POPC', allnames)
        self.assertIn('CHL1', allnames)
        matches = self.RM.search_resnames('popc')  # case-insensitive substring
        self.assertIn('POPC', matches)
        self.assertTrue(all('POPC' in m for m in matches))
        self.assertEqual(matches, sorted(matches))
        self.assertEqual(self.RM.search_resnames('ZZZZZZ'), [])


