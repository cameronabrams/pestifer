import unittest
import os
import tempfile
from pathlib import Path
from unittest import mock
from pestifer.core.resourcemanager import ResourceManager
from pestifer.core.errors import PestiferError
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
        # a residue defined in a stream/substream file: the cached index must be authoritative
        # (this one is missed by a naive topology-file scan)
        oct_ = self.RM.lookup_resname('OCT')
        self.assertTrue(oct_['in_topology'])
        self.assertEqual(oct_['kind'], 'residue (RESI)')
        # a name that exists nowhere
        bogus = self.RM.lookup_resname('ZZZZ')
        self.assertFalse(bogus['in_topology'])
        self.assertIsNone(bogus['topfile'])
        self.assertFalse(bogus['in_pdbrepository'])
        # richer detail for a PDB-repo lipid: charge present, head-to-tail length > 0
        self.assertIsNotNone(popc['charge'])
        self.assertGreater(popc['head_tail_length'], 0.0)

    def test_lookup_resname_source(self):
        # 'standard' = native to the CHARMM release
        self.assertEqual(self.RM.lookup_resname('POPC')['source'], 'standard')
        # 'custom' = a pestifer built-in custom file (charmmff/<ver>/custom/)
        self.assertEqual(self.RM.lookup_resname('PO4')['source'], 'custom')
        # 'user' = a user-custom directory (e.g. ~/.pestifer/toppar or a configured searchpath)
        import tempfile
        from pestifer.core.resourcemanager import ResourceManager
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'mylig.str'), 'w') as f:
                f.write('* test custom stream\n*\nRESI TESTLIG   0.000\nGROUP\nATOM C1 CG331 -0.27\n')
            rm = ResourceManager(charmmff_config={'user_custom': {'searchpath': [d]}})
            info = rm.lookup_resname('TESTLIG')
            self.assertTrue(info['in_topology'])
            self.assertEqual(info['source'], 'user')

    def test_add_custom_residue(self):
        # exercise the content mutation with its real-resource side effects redirected:
        # the custom file goes to a temp dir, the segtype edit and cache clear are stubbed
        d = tempfile.mkdtemp()
        src = os.path.join(d, 'mylig.str')
        with open(src, 'w') as f:
            f.write('* my ligand\n*\nRESI ZZL1  0.000\nGROUP\nATOM C1 CG331 -0.27\n')
        customdir = os.path.join(d, 'custom')
        with mock.patch.object(self.RM, 'get_charmmff_customdir', return_value=customdir), \
             mock.patch('pestifer.core.labels.register_resnames_segtype',
                        return_value=(Path('/fake/labels.py'), ['ZZL1'])), \
             mock.patch('pestifer.util.cacheable_object.CacheableObject.clear_cache',
                        return_value=[]):
            result = self.RM.add_custom_residue(src, segtype='ligand')
        self.assertEqual(result['resnames'], ['ZZL1'])
        self.assertEqual(result['segtype'], 'ligand')
        self.assertTrue(os.path.exists(os.path.join(customdir, 'mylig.str')))
        # both the copied file and the (stubbed) labels.py are staged for a commit
        self.assertIn(os.path.join(customdir, 'mylig.str'), result['touched_paths'])
        self.assertIn('/fake/labels.py', result['touched_paths'])

    def test_add_custom_residue_rejects_collision(self):
        # a RESI name that already exists in the force field must be refused unless forced
        d = tempfile.mkdtemp()
        src = os.path.join(d, 'dup.str')
        with open(src, 'w') as f:
            f.write('RESI POPC  0.000\nGROUP\nATOM C1 CG331 -0.27\n')
        with mock.patch.object(self.RM, 'get_charmmff_customdir', return_value=os.path.join(d, 'custom')):
            with self.assertRaises(PestiferError):
                self.RM.add_custom_residue(src, segtype='ligand')

    def test_add_custom_residue_rejects_pres_only(self):
        d = tempfile.mkdtemp()
        src = os.path.join(d, 'patch.str')
        with open(src, 'w') as f:
            f.write('PRES ZZP1  0.000\nATOM C1 CG331 -0.27\n')
        with mock.patch.object(self.RM, 'get_charmmff_customdir', return_value=os.path.join(d, 'custom')):
            with self.assertRaises(PestiferError):
                self.RM.add_custom_residue(src, segtype='ligand')

    def test_search_resnames(self):
        allnames = self.RM.all_resnames()
        self.assertIn('POPC', allnames)
        self.assertIn('CHL1', allnames)
        matches = self.RM.search_resnames('popc')  # case-insensitive substring
        self.assertIn('POPC', matches)
        self.assertTrue(all('POPC' in m for m in matches))
        self.assertEqual(matches, sorted(matches))
        self.assertEqual(self.RM.search_resnames('ZZZZZZ'), [])


