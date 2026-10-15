# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""The pre-flight steric check that keeps a RATTLE failure from reaching NAMD.

Synthetic structures only -- no VMD, no NAMD, no fixture files -- so this runs on CI, which is
the only place that sees a fresh checkout.
"""
import unittest
from types import SimpleNamespace

import numpy as np

from pestifer.util.contact_gate import (DEFAULT_MIN_CONTACT, describe_atom,
                                        worst_nonbonded_contact)


def _atom(serial, name, resid=1, resname='LIG', segname='A'):
    return SimpleNamespace(serial=serial, atomname=name, resid=resid, resname=resname,
                           segname=segname)


def _struct(atoms, bonds):
    """Duck-type of PSFContents: the gate uses only `.atoms` and `.bonds`."""
    return SimpleNamespace(atoms=atoms,
                           bonds=[SimpleNamespace(serial1=a, serial2=b) for a, b in bonds])


class TestWorstNonbondedContact(unittest.TestCase):

    def test_a_bonded_pair_is_not_an_overlap(self):
        """A C-C bond is 1.5 A.  Without the 1-2 exclusion every molecule reads as clashing --
        which is how the first draft of this measurement reported 21749 overlaps in a system
        that was fine."""
        s = _struct([_atom(1, 'C1'), _atom(2, 'C2')], [(1, 2)])
        xyz = np.array([[0.0, 0, 0], [1.5, 0, 0]])
        d, i, j = worst_nonbonded_contact(s, xyz)
        self.assertEqual(d, float('inf'), f'bonded pair reported as a contact at {d}')

    def test_a_1_3_pair_is_not_an_overlap(self):
        """Across a peptide bond, CA...N is 1-3 at ~2.4 A -- geometry, not a clash."""
        s = _struct([_atom(1, 'CA'), _atom(2, 'C'), _atom(3, 'N')], [(1, 2), (2, 3)])
        xyz = np.array([[0.0, 0, 0], [1.5, 0, 0], [2.4, 0, 0]])
        d, _, _ = worst_nonbonded_contact(s, xyz)
        self.assertEqual(d, float('inf'))

    def test_a_1_4_pair_IS_examined(self):
        """1-4 is where real steric strain lives, so it must not be excluded."""
        s = _struct([_atom(1, 'C1'), _atom(2, 'C2'), _atom(3, 'C3'), _atom(4, 'C4')],
                    [(1, 2), (2, 3), (3, 4)])
        xyz = np.array([[0.0, 0, 0], [1.5, 0, 0], [2.4, 0, 0], [0.6, 0, 0]])
        d, i, j = worst_nonbonded_contact(s, xyz)
        self.assertAlmostEqual(d, 0.6, places=6)
        self.assertEqual({i, j}, {0, 3})

    def test_two_unbonded_molecules_overlapping_are_found(self):
        s = _struct([_atom(1, 'C1', segname='A'), _atom(2, 'C2', segname='B')], [])
        xyz = np.array([[0.0, 0, 0], [0.35, 0, 0]])
        d, i, j = worst_nonbonded_contact(s, xyz)
        self.assertAlmostEqual(d, 0.35, places=6)
        self.assertLess(d, DEFAULT_MIN_CONTACT)

    def test_hydrogens_are_ignored(self):
        """RATTLE constrains X-H, but a legitimate H...H contact reaches ~1.8 A, so including
        hydrogens would force a choice between false positives and a useless threshold."""
        s = _struct([_atom(1, 'HA', segname='A'), _atom(2, 'HB', segname='B')], [])
        xyz = np.array([[0.0, 0, 0], [0.3, 0, 0]])
        d, _, _ = worst_nonbonded_contact(s, xyz)
        self.assertEqual(d, float('inf'))

    def test_the_minimum_image_convention_is_applied(self):
        """Two atoms either side of a box face are neighbours, not 19 A apart.  This is the
        defect that made `_check_fast` blind across periodic boundaries."""
        s = _struct([_atom(1, 'C1', segname='A'), _atom(2, 'C2', segname='B')], [])
        xyz = np.array([[0.2, 5.0, 5.0], [19.8, 5.0, 5.0]])
        box = np.array([20.0, 20.0, 20.0])
        d, _, _ = worst_nonbonded_contact(s, xyz, box=box)
        self.assertAlmostEqual(d, 0.4, places=6)
        d_free, _, _ = worst_nonbonded_contact(s, xyz, box=None)
        self.assertEqual(d_free, float('inf'), 'non-periodic scan should not wrap')

    def test_a_clean_structure_reports_infinity(self):
        s = _struct([_atom(i + 1, 'C', segname=chr(65 + i)) for i in range(4)], [])
        xyz = np.array([[0, 0, 0], [10, 0, 0], [0, 10, 0], [0, 0, 10]], dtype=float)
        d, i, j = worst_nonbonded_contact(s, xyz)
        self.assertEqual((d, i, j), (float('inf'), None, None))

    def test_coordinate_count_mismatch_raises(self):
        s = _struct([_atom(1, 'C1'), _atom(2, 'C2')], [])
        with self.assertRaises(ValueError):
            worst_nonbonded_contact(s, np.zeros((3, 3)))

    def test_describe_atom_names_the_residue(self):
        s = _struct([_atom(7, 'C18', resid=287, resname='CHL1', segname='L1')], [])
        self.assertEqual(describe_atom(s, 0), 'L1 CHL1 287 C18')


class TestThresholdDoesNotDriftFromTheSchema(unittest.TestCase):
    """The default lives in two places -- the Python constant and `base.yaml` -- and a reader
    who changes one will not think to change the other."""

    def test_schema_default_matches_the_module_constant(self):
        from importlib.resources import files as pkg_files
        import yaml
        spec = yaml.safe_load(
            (pkg_files('pestifer.schema') / 'base.yaml').read_text())

        def find(node):
            if isinstance(node, dict):
                if node.get('name') == 'contact_gate':
                    for a in node.get('attributes', []):
                        if a.get('name') == 'min_distance':
                            return a.get('default')
                for v in node.values():
                    r = find(v)
                    if r is not None:
                        return r
            elif isinstance(node, list):
                for v in node:
                    r = find(v)
                    if r is not None:
                        return r
            return None

        schema_default = find(spec)
        self.assertIsNotNone(schema_default, 'contact_gate.min_distance not found in base.yaml')
        self.assertAlmostEqual(float(schema_default), DEFAULT_MIN_CONTACT, places=6)


class TestTheGateIsReachableFromAnMdSpec(unittest.TestCase):
    """A schema key in the wrong block resolves to nothing and the gate silently never
    configures.  The first draft of this key landed under `steer` -- whose `ensemble` entry is
    worded identically to `md`'s -- and resolved to None."""

    def test_an_md_task_resolves_the_gate_with_its_defaults(self):
        from importlib.resources import files as pkg_files
        from pathlib import Path
        import yaml
        from pestifer.core.config import Config
        ex = (Path(str(pkg_files('pestifer.resources'))) / 'examples' / '19' / 'inputs'
              / 'sperm-whale-myoglobin.yaml')
        # userdict, not userfile: a `userfile` makes Config verify charmrun/namd3/vmd/catdcd
        # on PATH, and the clean-export gate hides them on purpose
        c = Config(userdict=yaml.safe_load(ex.read_text()) or {}, quiet=True).configure_new()
        md = [t['md'] for t in c['user']['tasks'] if 'md' in t]
        self.assertTrue(md, 'example 19 has no md task')
        gate = md[0].get('contact_gate')
        self.assertIsInstance(gate, dict, 'contact_gate did not resolve under the md task')
        self.assertTrue(gate.get('enabled'))
        self.assertAlmostEqual(float(gate['min_distance']), DEFAULT_MIN_CONTACT, places=6)


class TestTheGateRunsOncePerTaskNotOncePerChunk(unittest.TestCase):
    """``ChunkedEquilibrateTask`` calls ``namdrun()`` in a loop -- a real example-17 patch
    equilibration ran **1,105** chunks.  At ~10 s per check that is ~3 hours of overhead on one
    task, and it buys nothing: chunk N+1 starts from chunk N's output, and if that were
    overlapped NAMD would have died producing it.  The gate is for the hand-off INTO a task.

    The first version of this feature had no such guard.  It was caught by asking how often the
    call site fires, not by any test going red, which is why there is a test now.
    """

    def _task(self):
        from pestifer.tasks.mdtask import MDTask
        calls = []

        class _Probe(MDTask):
            def __init__(self):            # bypass BaseTask.__init__; this exercises one method
                self.taskname = 'md'

        t = _Probe()
        t._checked = calls
        return t, calls

    def test_second_and_later_calls_short_circuit(self):
        task, calls = self._task()
        state = SimpleNamespace(psf=None, pdb=None, coor=None, xsc=None)
        specs = {'contact_gate': {'enabled': True, 'min_distance': 1.8}}
        # state has no psf/coor, so the body returns early *after* setting the flag; what this
        # asserts is the flag's lifecycle, which is what bounds the cost
        task._gate_steric_overlaps(state, 'NVT', specs)
        self.assertTrue(getattr(task, '_contact_gate_done', False),
                        'first dynamics call did not mark the task as gated')
        task._gate_steric_overlaps(state, 'NVT', specs)   # must be a no-op
        self.assertTrue(task._contact_gate_done)

    def test_a_minimize_does_not_consume_the_one_check(self):
        """A minimize is skipped before the flag is set, so the dynamics stage after it is still
        gated -- otherwise a `minimize -> dynamics` pair would never be checked at all."""
        task, _ = self._task()
        state = SimpleNamespace(psf=None, pdb=None, coor=None, xsc=None)
        specs = {'contact_gate': {'enabled': True}}
        task._gate_steric_overlaps(state, 'minimize', specs)
        self.assertFalse(getattr(task, '_contact_gate_done', False),
                         'a minimize consumed the task\'s one steric check')

    def test_disabled_does_not_consume_the_one_check_either(self):
        task, _ = self._task()
        state = SimpleNamespace(psf=None, pdb=None, coor=None, xsc=None)
        task._gate_steric_overlaps(state, 'NVT', {'contact_gate': {'enabled': False}})
        self.assertFalse(getattr(task, '_contact_gate_done', False))
