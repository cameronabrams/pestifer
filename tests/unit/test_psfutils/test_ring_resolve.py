# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Tests for the shared glycan-pendant ring-piercing resolver (pestifer.psfutil.ring_resolve),
used both by RingCheckTask and by the psfgen task at glycan graft time.

The rotation path itself (clearing a real glycan-bond piercing) is validated end-to-end by a real
glycosylated build; here we test the detection/partition wiring: only a glycan-*bond* piercer of a
protein/glycan ring is resolvable by rotating the glycan, everything else is returned unresolved.
"""
import os
import unittest
from pathlib import Path
from unittest import mock

from pestifer.psfutil import ring_resolve

# reuse the ring_check fixture 5: a glycan pyranose ring pierced by a *lipid* bond
_FIX5 = (Path(__file__).parent.parent / 'test_tasks' / 'test_ringcheck' / '5').resolve()


class TestIsGlycanPierce(unittest.TestCase):
    def _p(self, piercee_type, piercer_type, bond_serials=(1, 2)):
        return {'piercee': {'segtype': piercee_type, 'segname': 'A', 'resid': 1},
                'piercer': {'segtype': piercer_type, 'segname': 'B', 'resid': 2,
                            'bond_serials': list(bond_serials)}}

    def test_glycan_bond_through_protein_ring_is_resolvable(self):
        self.assertTrue(ring_resolve._is_glycan_pierce(self._p('protein', 'glycan')))

    def test_glycan_bond_through_glycan_ring_is_resolvable(self):
        self.assertTrue(ring_resolve._is_glycan_pierce(self._p('glycan', 'glycan')))

    def test_lipid_piercer_not_resolvable(self):
        self.assertFalse(ring_resolve._is_glycan_pierce(self._p('glycan', 'lipid')))

    def test_glycan_piercer_without_bond_serials_not_resolvable(self):
        p = self._p('protein', 'glycan')
        p['piercer']['bond_serials'] = []
        self.assertFalse(ring_resolve._is_glycan_pierce(p))


class TestResolveGlycanPiercings(unittest.TestCase):

    def test_no_piercings_returns_none(self):
        with mock.patch.object(ring_resolve, 'ring_check', return_value=[]):
            fixed, unresolved = ring_resolve.resolve_glycan_piercings('x.psf', 'x.pdb', None, 'out')
        self.assertIsNone(fixed)
        self.assertEqual(unresolved, [])

    def test_non_glycan_piercer_left_unresolved_no_rotation(self):
        # fixture 5 is a glycan ring pierced by a LIPID bond -> the glycan-pendant path cannot
        # (and must not) touch it; it comes back unresolved with no rotation attempted
        psf, pdb, xsc = (str(_FIX5 / 'test.psf'), str(_FIX5 / 'test.pdb'), str(_FIX5 / 'test.xsc'))
        with mock.patch.object(ring_resolve, 'try_glycan_pendant') as rot:
            fixed, unresolved = ring_resolve.resolve_glycan_piercings(
                psf, pdb, xsc, 'out', segtypes=['lipid', 'glycan'])
            rot.assert_not_called()
        self.assertIsNone(fixed)
        self.assertEqual(len(unresolved), 1)
        self.assertEqual(unresolved[0]['piercer']['segtype'], 'lipid')

    def test_resolvable_piercing_is_rotated(self):
        piercing = {'piercee': {'segtype': 'protein', 'segname': 'A', 'resid': 5, 'resname': 'PRO'},
                    'piercer': {'segtype': 'glycan', 'segname': 'G', 'resid': 3, 'resname': 'BGNA',
                                'bond_serials': [10, 11]}}
        with mock.patch.object(ring_resolve, 'ring_check', return_value=[piercing]), \
             mock.patch.object(ring_resolve, 'RingChecker') as RC, \
             mock.patch.object(ring_resolve, 'cell_from_xsc', return_value=(None, None)), \
             mock.patch.object(ring_resolve, 'try_glycan_pendant', return_value='fixed.pdb') as rot:
            RC.return_value.load_coords.return_value = mock.Mock()
            fixed, unresolved = ring_resolve.resolve_glycan_piercings('x.psf', 'x.pdb', None, 'out')
            rot.assert_called_once()
        self.assertEqual(fixed, 'fixed.pdb')
        self.assertEqual(unresolved, [])

    def test_resolvable_but_rotation_fails_is_unresolved(self):
        piercing = {'piercee': {'segtype': 'glycan', 'segname': 'A', 'resid': 5, 'resname': 'BMAN'},
                    'piercer': {'segtype': 'glycan', 'segname': 'G', 'resid': 3, 'resname': 'BGNA',
                                'bond_serials': [10, 11]}}
        with mock.patch.object(ring_resolve, 'ring_check', return_value=[piercing]), \
             mock.patch.object(ring_resolve, 'RingChecker') as RC, \
             mock.patch.object(ring_resolve, 'cell_from_xsc', return_value=(None, None)), \
             mock.patch.object(ring_resolve, 'try_glycan_pendant', return_value=None):
            RC.return_value.load_coords.return_value = mock.Mock()
            fixed, unresolved = ring_resolve.resolve_glycan_piercings('x.psf', 'x.pdb', None, 'out')
        self.assertIsNone(fixed)
        self.assertEqual(len(unresolved), 1)


class TestPsfgenGating(unittest.TestCase):
    """PsfgenTask.resolve_glycan_piercings gating: toggle off / no glycans -> no work."""

    def _task(self, *, check=True, nglycans=1):
        from pestifer.tasks.psfgen import PsfgenTask
        t = PsfgenTask.__new__(PsfgenTask)
        t.specs = {'source': {'sequence': {'glycans': {'declash': {'check_piercings': check}}}}}
        t.declash_counts = {'glycan': nglycans}
        return t

    def test_toggle_off_skips(self):
        t = self._task(check=False)
        with mock.patch('pestifer.psfutil.ring_resolve.resolve_glycan_piercings') as m:
            t.resolve_glycan_piercings()
            m.assert_not_called()

    def test_no_glycans_skips(self):
        t = self._task(check=True, nglycans=0)
        with mock.patch('pestifer.psfutil.ring_resolve.resolve_glycan_piercings') as m:
            t.resolve_glycan_piercings()
            m.assert_not_called()


class TestPiercerWithin(unittest.TestCase):
    """_piercer_within: a piercing is mover-caused iff its bond lies wholly in the mover set."""

    def test_bond_within_movers(self):
        self.assertTrue(ring_resolve._piercer_within({'piercer': {'bond_serials': [10, 11]}}, {10, 11, 12}))

    def test_bond_partly_outside_movers(self):
        self.assertFalse(ring_resolve._piercer_within({'piercer': {'bond_serials': [10, 99]}}, {10, 11}))

    def test_no_bond_serials(self):
        self.assertFalse(ring_resolve._piercer_within({'piercer': {}}, {1, 2}))


class TestResolveModelBuiltPiercings(unittest.TestCase):

    def test_no_piercings_returns_none(self):
        with mock.patch.object(ring_resolve, 'ring_check', return_value=[]):
            fixed, unresolved = ring_resolve.resolve_model_built_piercings(
                'x.psf', 'x.pdb', None, 'out', {1, 2})
        self.assertIsNone(fixed)
        self.assertEqual(unresolved, [])

    def test_non_mover_piercing_is_ignored(self):
        # a piercing whose bond is NOT within the mover set is another pass's concern:
        # it must NOT be rotated here and must NOT be reported as unresolved.
        p = {'piercee': {'segtype': 'protein', 'segname': 'A', 'resid': 5},
             'piercer': {'segtype': 'glycan', 'segname': 'G', 'resid': 3, 'bond_serials': [100, 101]}}
        with mock.patch.object(ring_resolve, 'ring_check', return_value=[p]), \
             mock.patch.object(ring_resolve, 'try_glycan_pendant') as rot:
            fixed, unresolved = ring_resolve.resolve_model_built_piercings(
                'x.psf', 'x.pdb', None, 'out', {10, 11})
            rot.assert_not_called()
        self.assertIsNone(fixed)
        self.assertEqual(unresolved, [])

    def test_mover_piercing_is_rotated(self):
        # a tail backbone bond (serials within the mover set) spearing a proline ring is rotated
        p = {'piercee': {'segtype': 'protein', 'segname': 'A', 'resid': 5, 'resname': 'PRO'},
             'piercer': {'segtype': 'protein', 'segname': 'A', 'resid': 60, 'resname': 'GLY',
                         'bond_serials': [10, 11]}}
        with mock.patch.object(ring_resolve, 'ring_check', return_value=[p]), \
             mock.patch.object(ring_resolve, 'RingChecker') as RC, \
             mock.patch.object(ring_resolve, 'cell_from_xsc', return_value=(None, None)), \
             mock.patch.object(ring_resolve, 'try_glycan_pendant', return_value='fixed.pdb') as rot:
            RC.return_value.load_coords.return_value = mock.Mock()
            fixed, unresolved = ring_resolve.resolve_model_built_piercings(
                'x.psf', 'x.pdb', None, 'out', {10, 11}, mover_kind='tail')
            rot.assert_called_once()
            # the mover set is forwarded so only tail atoms may move
            self.assertEqual(rot.call_args.kwargs.get('mover_serials'), {10, 11})
        self.assertEqual(fixed, 'fixed.pdb')
        self.assertEqual(unresolved, [])


if __name__ == '__main__':
    unittest.main()


class TestPierceeTranslate(unittest.TestCase):
    """Sliding the pierced molecule -- the move that is not a rotation.

    Both existing remedies rotate: `try_glycan_pendant` bends the piercing tail,
    `try_sterol_spin` turns the speared ring.  Neither can free a molecule skewered at two
    points.  Example 17's CHL1 L1-197 was threaded through its A-ring (6-membered, the 3-OH
    end) by one acyl chain and its D-ring (5-membered, the tail end) by another, each crossing
    within 0.2 A of the ring center; all 903 and 917 candidate tail rotations left it threaded,
    and a sterol's long axis runs near both ring centers so no spin moves them.  Of 1022 rigid
    translations, 636 un-threaded it and 25 did so with no hard overlap.
    """

    def _checker(self, n_atoms=8, outside_bonds=False):
        """A stand-in RingChecker: `n_atoms` atoms of residue L1/197, optionally with a bond
        leaving the residue."""
        import networkx as nx
        import numpy as np
        serials = list(range(1, n_atoms + 1))
        atoms = []
        for s in serials:
            a = mock.Mock()
            a.serial = s
            a.segname = 'L1'
            a.resid.resid = 197
            atoms.append(a)
        G = nx.Graph()
        G.add_nodes_from(serials)
        G.add_edges_from(zip(serials, serials[1:]))
        if outside_bonds:
            G.add_edge(serials[-1], 999)      # a bond to another residue
        ck = mock.Mock()
        ck._atoms = atoms
        ck._row_of_serial = {s: i for i, s in enumerate(serials)}
        ck.topol.G = G
        return ck

    def test_a_residue_bonded_outside_itself_is_refused(self):
        """The guard that makes this safe: translating one residue of a covalent chain would
        stretch the bonds either side of it into nonsense.  A lipid is its own molecule and
        passes; a protein or glycan residue must not."""
        ck = self._checker(outside_bonds=True)
        out = ring_resolve.try_piercee_translate(
            ck, 'in.pdb', mock.Mock(), None,
            {'segname': 'L1', 'resid': 197, 'resname': 'CHL1'}, ('L1', 197), 'out.pdb')
        self.assertIsNone(out)
        ck.check_coords.assert_not_called()      # refused before any geometry was tried

    def test_a_free_molecule_is_translated(self):
        import numpy as np
        ck = self._checker()
        ck.check_coords.return_value = False     # every pose un-threads
        ck.clash_count.return_value = 0          # and is clean
        base = np.zeros((8, 3))
        with mock.patch.object(ring_resolve, 'pdb_replace_coords') as wr:
            out = ring_resolve.try_piercee_translate(
                ck, 'in.pdb', base, None,
                {'segname': 'L1', 'resid': 197, 'resname': 'CHL1'}, ('L1', 197), 'out.pdb')
        self.assertEqual(out, 'out.pdb')
        wr.assert_called_once()

    def test_the_smallest_displacement_that_works_wins(self):
        """A piercing is cleared by the least motion that clears it: the search must stop at
        the first magnitude that yields an eligible pose, not sweep them all and rank."""
        import numpy as np
        ck = self._checker()
        base = np.zeros((8, 3))
        seen = []

        def check_coords(scratch, box, targets):
            seen.append(float(np.linalg.norm(scratch[0])))
            return False                          # everything un-threads
        ck.check_coords.side_effect = check_coords
        ck.clash_count.return_value = 0
        with mock.patch.object(ring_resolve, 'pdb_replace_coords'):
            ring_resolve.try_piercee_translate(
                ck, 'in.pdb', base, None,
                {'segname': 'L1', 'resid': 197, 'resname': 'CHL1'}, ('L1', 197), 'out.pdb')
        self.assertTrue(seen)
        self.assertAlmostEqual(max(seen), ring_resolve._TRANSLATE_MAGNITUDES[0], places=6,
                               msg='the search went past the first magnitude that worked')

    def test_a_pose_inside_the_steric_gate_is_rejected(self):
        """Same eligibility rule as the rotations: burying the molecule in a neighbor is not
        a resolution, and the gate downstream refuses it anyway."""
        import numpy as np
        ck = self._checker()
        ck.check_coords.return_value = False
        ck.clash_count.side_effect = lambda *a, **k: 1 if k.get('cutoff') else 0
        base = np.zeros((8, 3))
        with mock.patch.object(ring_resolve, 'pdb_replace_coords') as wr:
            out = ring_resolve.try_piercee_translate(
                ck, 'in.pdb', base, None,
                {'segname': 'L1', 'resid': 197, 'resname': 'CHL1'}, ('L1', 197), 'out.pdb')
        self.assertIsNone(out)
        wr.assert_not_called()

    def test_the_direction_set_is_deterministic(self):
        """Not an RNG draw.  A resolver that gives a different answer on a rerun cannot be
        debugged, and this repository has already paid for that once in ligate's loop closure."""
        import numpy as np
        a, b = ring_resolve._unit_sphere(98), ring_resolve._unit_sphere(98)
        self.assertTrue(np.array_equal(a, b))
        self.assertEqual(a.shape, (98, 3))
        self.assertTrue(np.allclose(np.linalg.norm(a, axis=1), 1.0))


class TestAFailedSpinSaysSo(unittest.TestCase):
    """`try_sterol_spin` used to return None without a word at any level, so a failed spin left
    no trace -- on 3.29.5 its absence read as "the spin was never tried"."""

    def test_failure_is_logged(self):
        import numpy as np
        ck = mock.Mock()
        atoms = []
        for s in range(1, 9):
            a = mock.Mock(); a.serial = s; a.segname = 'L1'; a.resid.resid = 199; atoms.append(a)
        ck._atoms = atoms
        ck._row_of_serial = {s: s - 1 for s in range(1, 9)}
        ck.check_coords.return_value = True            # every angle leaves it pierced
        base = np.random.default_rng(0).normal(size=(8, 3))
        with self.assertLogs('pestifer.psfutil.ring_resolve', level='DEBUG') as cm:
            out = ring_resolve.try_sterol_spin(ck, 'in.pdb', base, None,
                                               {'segname': 'L1', 'resid': 199, 'resname': 'CHL1'},
                                               ('L1', 199), 'out.pdb', 'SOPE L1-109')
        self.assertIsNone(out)
        self.assertTrue(any('no spin' in l for l in cm.output), cm.output)
