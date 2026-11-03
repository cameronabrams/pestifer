# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""`merge` checks the interfaces it creates, which no input was built or checked against.

Example 25 merges an Env trimer built with its 17b Fabs excluded and three Fabs built on their
own.  The Env build models the V3 tip, which 5vn3 does not resolve, into space a Fab occupies:
every replica of the 3.29.8 and 3.30.0 sweeps interpenetrated there, and in one a Tyr ring and
a Phe ring ended up threaded through each other.  Nothing looked until `terminate`, hours later.

No VMD, no NAMD: the PSF, the coordinates and `ring_check` are stubbed, because what is under
test is which findings are cross-system and what the task does about them.  The geometry
engines have their own tests; `test_merge.py` runs a real merge through the call site.
"""
import logging
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from pestifer.core.errors import PestiferBuildError
from pestifer.tasks.merge import MergeTask


def _atom(seg, resname, resid, name):
    return SimpleNamespace(segname=seg, resname=resname, resid=resid, atomname=name)


# two segments of system 0 (A, B) and one of system 1 (H); A:1 CA sits on H:9 CB
ATOMS = [_atom('A', 'ALA', 1, 'CA'), _atom('A', 'ALA', 1, 'HA'),
         _atom('B', 'GLY', 5, 'CA'), _atom('H', 'TYR', 9, 'CB'), _atom('H', 'TYR', 9, 'HB1')]
XYZ = np.array([[0.0, 0, 0], [0.1, 0, 0], [-1.0, 0, 0], [0.3, 0, 0], [50.0, 0, 0]])
SYSTEM_OF = {'A': 0, 'B': 0, 'H': 1}


def _piercing(pe_seg, pr_seg):
    return {'piercee': {'segname': pe_seg, 'resname': 'TYR', 'resid': 100},
            'piercer': {'segname': pr_seg, 'resname': 'PHE', 'resid': 316}}


def _check(piercings, mode='error', atoms=ATOMS, xyz=XYZ):
    t = MergeTask.__new__(MergeTask)
    with patch('pestifer.tasks.merge.PSFContents', return_value=SimpleNamespace(atoms=atoms)), \
         patch('pestifer.util.density_convergence._read_coor_xyz', return_value=xyz), \
         patch('pestifer.psfutil.psfring.ring_check', return_value=piercings) as rc:
        t._check_interfaces('m.psf', 'm.pdb', SYSTEM_OF, mode)
    return rc


class TestCrossSystemPiercing(unittest.TestCase):

    def test_a_ring_pierced_from_another_system_stops_the_build(self):
        with self.assertRaises(PestiferBuildError) as cm:
            _check([_piercing('H', 'A')])
        msg = str(cm.exception)
        self.assertIn('RING_PIERCING_FOUND', msg)
        for token in ('H TYR 100', 'A PHE 316', 'interface_check: warn'):
            self.assertIn(token, msg)

    def test_a_piercing_inside_one_system_is_not_merges_to_report(self):
        """A and B are both system 0: that piercing came in with the input, and the input's
        own checks own it.  Raising here would blame the merge for it."""
        _check([_piercing('A', 'B')])   # must not raise

    def test_warn_reports_and_continues(self):
        with self.assertLogs('pestifer.tasks.merge', level='WARNING') as cm:
            _check([_piercing('H', 'A')], mode='warn')
        self.assertTrue(any('RING_PIERCING_FOUND' in m for m in cm.output))

    def test_every_crossing_is_listed(self):
        with self.assertRaises(PestiferBuildError) as cm:
            _check([_piercing('H', 'A'), _piercing('A', 'H'), _piercing('A', 'B')])
        self.assertIn('2 ring(s)', str(cm.exception))

    def test_protein_rings_are_checked(self):
        rc = _check([])
        self.assertEqual(rc.call_args.kwargs['segtypes'], ['lipid', 'glycan', 'protein'])


class TestCrossSystemOverlap(unittest.TestCase):

    def test_an_overlap_between_systems_is_named(self):
        with self.assertLogs('pestifer.tasks.merge', level='WARNING') as cm:
            _check([])
        out = '\n'.join(cm.output)
        self.assertIn('A ALA 1 / H TYR 9: 1', out)

    def test_an_overlap_inside_one_system_is_not_reported(self):
        """A:1 CA and B:5 CA are 1.0 A apart but both in system 0; hydrogens never count."""
        with self.assertLogs('pestifer.tasks.merge', level='WARNING') as cm:
            _check([])
        out = '\n'.join(cm.output)
        self.assertNotIn('GLY 5', out)
        self.assertIn('1 heavy-atom pair(s)', out)

    def test_no_overlap_no_warning(self):
        far = XYZ.copy()
        far[3] = [20.0, 0, 0]
        with self.assertNoLogs('pestifer.tasks.merge', level='WARNING'):
            _check([], xyz=far)
