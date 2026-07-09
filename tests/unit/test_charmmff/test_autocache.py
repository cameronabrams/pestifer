# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Tests for the on-demand PDB-repository generation cache (pestifer.charmmff.autocache).

Unit tests (no NAMD): cache path helpers, cache-hit short-circuit, and the solvate task's
generate-on-miss policy (toggle + not-in-force-field guards).  The actual box build is a slow
NAMD equilibration and is exercised in test_solvate.py's slow integration tier.
"""
import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

import pytest
import yaml

from pestifer.charmmff import autocache


def _write_fake_box_entry(collection_dir: Path, resname: str):
    """Create a minimal ``kind: box`` cache entry (info.yaml + placeholder psf/pdb).

    The directory-collection reader only reads file *contents*; it does not parse them, so
    placeholder psf/pdb are sufficient to exercise registration/checkout wiring."""
    entry = collection_dir / resname
    entry.mkdir(parents=True)
    (entry / f'{resname}-box.psf').write_text('PSF\n')
    (entry / f'{resname}-box.pdb').write_text('END\n')
    (entry / 'info.yaml').write_text(yaml.dump({
        'kind': 'box', 'resname': resname, 'nmol': 8, 'box_edge': 12.0, 'density': 0.9,
        'key_atom': 'C1', 'quality': 'auto',
        'psf': f'{resname}-box.psf', 'pdb': f'{resname}-box.pdb',
    }))


def _fake_cc(release_key='feb26', release_str='February2026'):
    """A stand-in CHARMMFFContent carrying just what autocache derives from it."""
    cc = mock.Mock()
    cc.charmmff_path = f'/x/{release_key}'
    cc.release_str = release_str
    return cc


class TestCachePaths(unittest.TestCase):

    def setUp(self):
        self._tmp = TemporaryDirectory()
        self._patch = mock.patch.object(autocache, 'PDBCACHE_ROOT', Path(self._tmp.name))
        self._patch.start()

    def tearDown(self):
        self._patch.stop()
        self._tmp.cleanup()

    def test_cache_release_root(self):
        self.assertEqual(autocache.cache_release_root('feb26'),
                         Path(self._tmp.name) / 'feb26')

    def test_cached_collection_dirs_empty(self):
        self.assertEqual(autocache.cached_collection_dirs('feb26'), [])

    def test_cached_collection_dirs_skips_dotdirs_and_files(self):
        root = autocache.cache_release_root('feb26')
        (root / 'solvent').mkdir(parents=True)
        (root / 'lipid').mkdir()
        (root / '.tmp').mkdir()          # hidden work/lock dirs must be skipped
        (root / 'note.txt').write_text('x')  # files must be skipped
        dirs = [p.name for p in autocache.cached_collection_dirs('feb26')]
        self.assertEqual(dirs, ['lipid', 'solvent'])

    def test_has_entry(self):
        coll = autocache.cache_release_root('feb26') / 'solvent'
        (coll / 'DMSO').mkdir(parents=True)
        self.assertFalse(autocache._has_entry(coll, 'DMSO'))
        (coll / 'DMSO' / 'info.yaml').write_text('kind: box\n')
        self.assertTrue(autocache._has_entry(coll, 'DMSO'))


class TestEnsureSolventBoxCacheHit(unittest.TestCase):
    """A pre-existing cache entry must be returned without invoking the (slow) builder."""

    def setUp(self):
        self._tmp = TemporaryDirectory()
        self._patch = mock.patch.object(autocache, 'PDBCACHE_ROOT', Path(self._tmp.name))
        self._patch.start()

    def tearDown(self):
        self._patch.stop()
        self._tmp.cleanup()

    def test_cache_hit_does_not_build(self):
        coll = autocache.cache_release_root('feb26') / 'solvent'
        (coll / 'GLYE').mkdir(parents=True)
        (coll / 'GLYE' / 'info.yaml').write_text(yaml.dump({'kind': 'box', 'resname': 'GLYE'}))
        with mock.patch('pestifer.charmmff.make_solvent_box.make_solvent_box') as m:
            got = autocache.ensure_solvent_box('GLYE', _fake_cc())
            m.assert_not_called()
        self.assertEqual(got, coll)

    def test_lipid_cache_hit_does_not_build(self):
        coll = autocache.cache_release_root('feb26') / 'lipid'
        (coll / 'POPX').mkdir(parents=True)
        (coll / 'POPX' / 'info.yaml').write_text(yaml.dump({'kind': 'molecule', 'conformers': []}))
        with mock.patch('pestifer.charmmff.make_pdb_collection.do_resi') as m:
            got = autocache.ensure_lipid_conformer('POPX', _fake_cc())
            m.assert_not_called()
        self.assertEqual(got, coll)


class TestSolvateGenerateOnMissPolicy(unittest.TestCase):
    """SolvateTask._generate_solvent_box: toggle + not-in-force-field guards (no build)."""

    def _task(self):
        from pestifer.tasks.solvate import SolvateTask
        return SolvateTask.__new__(SolvateTask)

    def _cc(self, *, toggle=True, in_ff=True):
        cc = mock.Mock()
        cc.generate_missing_coordinates = toggle
        cc.__contains__ = mock.Mock(return_value=in_ff)
        cc.charmmff_path = '/x/feb26'
        return cc

    def test_toggle_off_raises_without_building(self):
        from pestifer.core.errors import PestiferError
        task = self._task()
        cc = self._cc(toggle=False)
        with mock.patch('pestifer.charmmff.autocache.ensure_solvent_box') as m:
            with self.assertRaises(PestiferError) as ctx:
                task._generate_solvent_box('DMSO', cc, mock.Mock())
            m.assert_not_called()
        self.assertIn('generate_missing_coordinates', str(ctx.exception))

    def test_not_in_force_field_raises(self):
        # a solvent absent from the force field surfaces as ValueError from ensure_solvent_box
        # (built against a fresh, fully-provisioned CC), wrapped as PestiferError
        from pestifer.core.errors import PestiferError
        task = self._task()
        cc = self._cc(toggle=True)
        with mock.patch('pestifer.charmmff.autocache.ensure_solvent_box',
                        side_effect=ValueError('RESI NOPE is not defined in the CHARMM force field')):
            with self.assertRaises(PestiferError):
                task._generate_solvent_box('NOPE', cc, mock.Mock())

    def test_miss_builds_and_registers(self):
        task = self._task()
        cc = self._cc(toggle=True, in_ff=True)
        repo = mock.Mock()
        # after add_resource, the solvent is now present
        repo.__contains__ = mock.Mock(return_value=True)
        with mock.patch('pestifer.charmmff.autocache.ensure_solvent_box',
                        return_value=Path('/cache/feb26/solvent')) as m:
            task._generate_solvent_box('DMSO', cc, repo)
        m.assert_called_once_with('DMSO', cc)
        repo.add_resource.assert_called_once_with('/cache/feb26/solvent')


class TestBilayerGenerateOnMissPolicy(unittest.TestCase):
    """Bilayer._generate_missing_species: toggle + not-in-force-field guards (no build)."""

    def _bilayer(self, *, toggle=True, in_ff=True):
        from pestifer.molecule.bilayer import Bilayer
        b = Bilayer.__new__(Bilayer)
        cc = mock.Mock()
        cc.generate_missing_coordinates = toggle
        cc.__contains__ = mock.Mock(return_value=in_ff)
        b.charmmffcontent = cc
        return b

    def test_toggle_off_raises_without_building(self):
        from pestifer.core.errors import PestiferBuildError
        b = self._bilayer(toggle=False)
        with mock.patch('pestifer.charmmff.autocache.ensure_lipid_conformer') as m:
            with self.assertRaises(PestiferBuildError) as ctx:
                b._generate_missing_species('POPX', mock.Mock())
            m.assert_not_called()
        self.assertIn('generate_missing_coordinates', str(ctx.exception))

    def test_not_in_force_field_raises(self):
        from pestifer.core.errors import PestiferBuildError
        b = self._bilayer(toggle=True, in_ff=False)
        with mock.patch('pestifer.charmmff.autocache.ensure_lipid_conformer') as m:
            with self.assertRaises(PestiferBuildError):
                b._generate_missing_species('NOPE', mock.Mock())
            m.assert_not_called()

    def test_miss_builds_and_registers(self):
        b = self._bilayer(toggle=True, in_ff=True)
        repo = mock.Mock()
        repo.__contains__ = mock.Mock(return_value=True)  # present after registration
        with mock.patch('pestifer.charmmff.autocache.ensure_lipid_conformer',
                        return_value=Path('/cache/feb26/lipid')) as m:
            b._generate_missing_species('POPX', repo)
        m.assert_called_once_with('POPX', b.charmmffcontent)
        repo.add_resource.assert_called_once_with('/cache/feb26/lipid')


class TestProvisionAutoRegistersCache(unittest.TestCase):
    """A previously generated cache collection must be auto-registered at provision time, so
    cached entries are available on a subsequent run with no regeneration."""

    def setUp(self):
        self._tmp = TemporaryDirectory()
        self._patch = mock.patch.object(autocache, 'PDBCACHE_ROOT', Path(self._tmp.name))
        self._patch.start()

    def tearDown(self):
        self._patch.stop()
        self._tmp.cleanup()

    def test_cached_entry_registered_and_checkoutable(self):
        from pestifer.core.resourcemanager import ResourceManager
        RM = ResourceManager(charmmff_config={})
        CC = RM.charmmff_content
        release_key = os.path.basename(str(CC.charmmff_path))
        _write_fake_box_entry(autocache.cache_release_root(release_key) / 'solvent', 'FAKESOLV')
        CC.provision_pdbrepository()
        repo = CC.pdbrepository
        self.assertIn('FAKESOLV', repo)
        entry = repo.checkout('FAKESOLV')
        self.assertTrue(entry.is_box())
        # shipped collections must still be present
        self.assertIn('DMSO', repo)


class TestEnsureBuildIntegration(unittest.TestCase):
    """Slow: actually build cache entries for real force-field residues absent from the shipped
    repository, and confirm the cache-hit short-circuit on the second call."""

    def _live_cc(self):
        from pestifer.core.resourcemanager import ResourceManager
        CC = ResourceManager(charmmff_config={}).charmmff_content
        CC.provision()
        return CC

    @pytest.mark.slow
    def test_build_benz_box(self):
        with TemporaryDirectory() as tmp:
            with mock.patch.object(autocache, 'PDBCACHE_ROOT', Path(tmp)):
                cc = self._live_cc()
                coll = autocache.ensure_solvent_box('BENZ', cc, nmol=216, npt_steps=500,
                                                    minimize_steps=200)
                info = yaml.safe_load((coll / 'BENZ' / 'info.yaml').read_text())
                self.assertEqual(info['kind'], 'box')
                self.assertEqual(info['quality'], 'auto')
                self.assertGreater(info['box_edge'], 0)
                with mock.patch('pestifer.charmmff.make_solvent_box.make_solvent_box') as m:
                    coll2 = autocache.ensure_solvent_box('BENZ', cc)
                    m.assert_not_called()
                self.assertEqual(coll2, coll)

    @pytest.mark.slow
    def test_build_lipid_conformer(self):
        # POPC is shipped; generate conformers for it into an isolated cache to exercise the
        # kind:molecule path end-to-end without depending on an unshipped lipid name
        with TemporaryDirectory() as tmp:
            with mock.patch.object(autocache, 'PDBCACHE_ROOT', Path(tmp)):
                cc = self._live_cc()
                coll = autocache.ensure_lipid_conformer('POPC', cc, nsamples=2, sample_steps=200,
                                                        minimize_steps=100)
                info = yaml.safe_load((coll / 'POPC' / 'info.yaml').read_text())
                self.assertEqual(info.get('kind', 'molecule'), 'molecule')
                self.assertEqual(info['quality'], 'auto')
                self.assertGreaterEqual(len(info['conformers']), 1)
                self.assertIn('reference-atoms', info)
                with mock.patch('pestifer.charmmff.make_pdb_collection.do_resi') as m:
                    coll2 = autocache.ensure_lipid_conformer('POPC', cc)
                    m.assert_not_called()
                self.assertEqual(coll2, coll)


if __name__ == '__main__':
    unittest.main()
