# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""``source.obstacles``: atoms a build leaves out but its modelers must keep clear of.

Example 25 builds its Env trimer with the 17b Fabs excluded, because the Fabs are completed
separately and merged in afterwards.  The Env build models the V3 tip, which 5vn3 does not
resolve -- and with the Fabs absent, nothing told the loop closer where they would be.  In every
replica of the 3.29.8 and 3.30.0 sweeps the modeled tip overlapped a Fab, and in one the merge
locked a Tyr ring and a Phe ring through each other.

These tests drive the real call sites -- the CCD closer, the tail modeler and the loop declash --
and capture what reaches the modeler, because a correct ``obstacle_coords`` that a call site
does not use is exactly the failure this repository keeps a list of.
"""
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
from pidibble.pdbparse import PDBParser

from pestifer.core.errors import PestiferBuildError
from pestifer.core.objmanager import ObjManager
from pestifer.molecule.asymmetricunit import AsymmetricUnit
from pestifer.molecule.chainidmanager import ChainIDManager
from pestifer.molecule.molecule import Molecule
from pestifer.molecule.transform import Transform
from pestifer.tasks.ligate import LigateTask
from pestifer.tasks.psfgen import PsfgenTask

OBS = np.array([[100.0, 0.0, 0.0], [101.0, 0.0, 0.0]])


class TestObstaclesAreReadFromTheSource(unittest.TestCase):

    def _au(self, sourcespecs):
        pdb = Path(__file__).parents[2] / 'inputs' / '6pti.pdb'
        parsed = PDBParser(filepath=pdb, input_format='PDB').parse().parsed
        return AsymmetricUnit(parsed=parsed, chainIDmanager=ChainIDManager(format='PDB'),
                              objmanager=ObjManager(), sourcespecs=sourcespecs,
                              source_format='PDB')

    def test_excluded_atoms_can_be_obstacles(self):
        """The usual case: left out of the build, still in the way."""
        au = self._au({'exclude': ["resname == 'HOH'"], 'obstacles': ["resname == 'HOH'"]})
        self.assertGreater(len(au.obstacles), 0)
        self.assertTrue(all(a.resname == 'HOH' for a in au.obstacles.data))
        built = {r.resname for s in au.segments for r in s.residues.data}
        self.assertNotIn('HOH', built)

    def test_only_heavy_atoms(self):
        au = self._au({'obstacles': ["resname == 'CYS'"]})
        names = {a.name for a in au.obstacles.data}
        self.assertIn('SG', names)
        self.assertFalse(any(n.startswith('H') for n in names))

    def test_a_selection_matching_nothing_stops_the_build(self):
        """An obstacle list that matches nothing protects nothing, silently."""
        with self.assertRaises(PestiferBuildError):
            self._au({'obstacles': ["chainID == 'Z'"]})

    def test_none_declared_is_empty(self):
        self.assertEqual(len(self._au({}).obstacles), 0)


class TestObstaclesFollowTheAssembly(unittest.TestCase):

    def _mol(self, transforms):
        atoms = [SimpleNamespace(x=1.0, y=2.0, z=3.0, chainID='H'),
                 SimpleNamespace(x=4.0, y=5.0, z=6.0, chainID='L')]
        m = Molecule.__new__(Molecule)
        m.asymmetric_unit = SimpleNamespace(obstacles=_Atoms(atoms))
        m.active_biological_assembly = SimpleNamespace(transforms=transforms)
        return m

    def test_every_image_carries_its_obstacles(self):
        shift = np.identity(4)
        shift[:3, 3] = [10.0, 0, 0]
        T0 = Transform(index=0, tmat=np.identity(4), applies_chainIDs=[])
        T1 = Transform(index=1, tmat=shift, applies_chainIDs=['H'])
        xyz = self._mol([T0, T1]).obstacle_coords()
        self.assertEqual(xyz.shape, (3, 3), 'two from the identity image, one (chain H) from the other')
        self.assertTrue(any(np.allclose(r, [11.0, 2.0, 3.0]) for r in xyz))

    def test_none_declared_is_none(self):
        m = self._mol([])
        m.asymmetric_unit.obstacles = _Atoms([])
        self.assertIsNone(m.obstacle_coords())


class _Atoms:
    def __init__(self, data):
        self.data = data

    def __len__(self):
        return len(self.data)


class _Stop(Exception):
    pass


def _gap_task(cls):
    t = cls.__new__(cls)
    t.provisions = {'namd_global_config': {}}
    t.outcome = {}
    t.basename = 'b'
    t.next_basename = mock.Mock()
    state = mock.Mock()
    state.pdb.name = 'in.pdb'
    t.get_current_artifact = lambda k: state
    t.base_molecule = mock.Mock()
    t.base_molecule.obstacle_coords.return_value = OBS
    return t


def _has_obstacles(env):
    return env is not None and all(any(np.allclose(r, o) for r in env) for o in OBS)


class TestEveryModelerSeesThem(unittest.TestCase):

    def test_the_ccd_closer(self):
        t = _gap_task(LigateTask)
        t.base_molecule.protein_loop_gaps.return_value = [
            {'segname': 'A', 'loop_resids': [10, 11], 'n_anchor_resid': 9, 'c_anchor_resid': 12}]
        seen = {}

        def close_one_loop(**kw):
            seen['env'] = kw['extra_env']
            raise _Stop()
        with mock.patch('pestifer.psfutil.loop_ccd.loop_atoms_from_pdb', return_value=([], [], [])), \
             mock.patch('pestifer.psfutil.loop_ccd.close_one_loop', side_effect=close_one_loop):
            with self.assertRaises(_Stop):
                t.close_loops_ccd({})
        self.assertTrue(_has_obstacles(seen['env']))

    def test_the_tail_modeler(self):
        t = _gap_task(PsfgenTask)
        t.specs = {'source': {'sequence': {'loops': {'declash': {}}}}}
        t.base_molecule.protein_terminal_tails.return_value = [
            {'segname': 'A', 'tail_resids': [1, 2], 'anchor_resid': 3, 'end': 'N'}]
        t.base_molecule.protein_loop_gaps.return_value = []
        t._crotation_targeted_resids = lambda: {}
        seen = {}

        def model_one_tail(*a, **kw):
            seen['env'] = kw['extra_env']
            raise _Stop()
        with mock.patch('pestifer.psfutil.loop_ccd.loop_atoms_from_pdb', return_value=([], [], [])), \
             mock.patch('pestifer.psfutil.tail_model.model_one_tail', side_effect=model_one_tail):
            with self.assertRaises(_Stop):
                t.model_terminal_tails()
        self.assertTrue(_has_obstacles(seen['env']))

    def test_the_loop_declash_scores_against_them_and_never_moves_them(self):
        class _CM:
            def __init__(self, *a):
                self.coords = np.zeros((2, 3))
                self._mass = np.array([12.0, 12.0])
                self._names = np.array(['N', 'CA'])
                self.atoms = mock.Mock(data=[mock.Mock(segname='A', resid='1')] * 2)
                written['cm'] = self

            def _residue_index(self):
                return np.zeros(2, dtype=int)

            def residue_of_segname(self, *a):
                return 0

            def write_pdb(self, *a):
                pass
        written, seen = {}, {}

        def declash_loop(coords, residues, env_heavy, *a):
            seen['env'] = coords[env_heavy]
        t = _gap_task(PsfgenTask)
        with mock.patch('pestifer.molecule.coordmanip.CoordManipulator', _CM), \
             mock.patch('pestifer.psfutil.declash.declash_loop', side_effect=declash_loop):
            t._run_loop_declash(mock.Mock(), [('A', ['1'])], 1)
        self.assertTrue(_has_obstacles(seen['env']))
        self.assertEqual(written['cm'].coords.shape, (2, 3), 'obstacle rows must not be written out')


class TestOrientIsRefused(unittest.TestCase):

    def test_orient_with_obstacles_stops(self):
        t = _gap_task(PsfgenTask)
        t.objmanager = {'coord': {'orient': [object()]}}
        with self.assertRaises(PestiferBuildError):
            t.coormods()
