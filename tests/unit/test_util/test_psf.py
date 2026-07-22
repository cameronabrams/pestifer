# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Unit tests for pestifer.util.psf.dedupe_psf_topology.
"""

import os
import shutil
import tempfile
import unittest

from pestifer.util.psf import dedupe_psf_topology


def _section(name, k, records, width=10):
    """Render a counted PSF integer section (header + data) the way psfgen does."""
    gpl = {'NBOND': 4, 'NTHETA': 3, 'NPHI': 2, 'NIMPHI': 2}[name]
    label = {'NBOND': 'NBOND: bonds', 'NTHETA': 'NTHETA: angles',
             'NPHI': 'NPHI: dihedrals', 'NIMPHI': 'NIMPHI: impropers'}[name]
    lines = [f'{len(records):>{width}} !{label}']
    buf = []
    for i, rec in enumerate(records):
        buf.extend(rec)
        if (i + 1) % gpl == 0:
            lines.append(''.join(f'{v:>{width}}' for v in buf))
            buf = []
    if buf:
        lines.append(''.join(f'{v:>{width}}' for v in buf))
    return lines


def _make_psf(bonds, angles, dihedrals=None, natom=4):
    dihedrals = dihedrals or []
    lines = ['PSF EXT CMAP', '']
    lines += [f'{2:>10} !NTITLE', '* test psf', '*', '']
    lines += [f'{natom:>10} !NATOM']
    for a in range(1, natom + 1):
        # atom records are copied verbatim by the cleaner, so exact columns are unimportant
        lines.append(f'{a:>10} A        1    ALA  X    C      0.000000       12.0000           0')
    lines.append('')
    lines += _section('NBOND', 2, bonds) + ['']
    lines += _section('NTHETA', 3, angles) + ['']
    lines += _section('NPHI', 4, dihedrals) + ['']
    return '\n'.join(lines) + '\n'


def _read_section(path, marker, k):
    lines = open(path).read().splitlines()
    for i, l in enumerate(lines):
        if marker in l:
            count = int(l.split()[0])
            toks, j = [], i + 1
            while len(toks) < count * k:
                toks += [int(x) for x in lines[j].split()]
                j += 1
            recs = [tuple(toks[p:p + k]) for p in range(0, count * k, k)]
            return count, recs
    return None, []


class TestDedupePsfTopology(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_removes_duplicate_bond_and_degenerate_angle(self):
        path = os.path.join(self.tmp, 'bad.psf')
        # bond (1,2) appears twice; angle (1,2,1) is degenerate; the rest are valid
        bonds = [(1, 2), (2, 3), (3, 4), (1, 2)]
        angles = [(1, 2, 3), (1, 2, 1), (2, 3, 4)]
        # a dihedral duplicated in reverse order, plus a degenerate one
        dihedrals = [(1, 2, 3, 4), (4, 3, 2, 1), (1, 2, 3, 1)]
        with open(path, 'w') as fh:
            fh.write(_make_psf(bonds, angles, dihedrals))

        removed = dedupe_psf_topology(path)
        self.assertEqual(removed.get('NBOND'), (0, 1))    # one duplicate bond
        self.assertEqual(removed.get('NTHETA'), (1, 0))   # one degenerate angle
        self.assertEqual(removed.get('NPHI'), (1, 1))     # one degenerate + one reverse-duplicate

        nb, bset = _read_section(path, '!NBOND', 2)
        self.assertEqual(nb, 3)
        self.assertEqual(len(bset), 3)
        na, aset = _read_section(path, '!NTHETA', 3)
        self.assertEqual(na, 2)
        self.assertTrue(all(len(set(a)) == 3 for a in aset))
        nd, dset = _read_section(path, '!NPHI', 4)
        self.assertEqual(nd, 1)

    def test_clean_psf_is_untouched(self):
        path = os.path.join(self.tmp, 'clean.psf')
        content = _make_psf([(1, 2), (2, 3), (3, 4)],
                            [(1, 2, 3), (2, 3, 4)],
                            [(1, 2, 3, 4)])
        with open(path, 'w') as fh:
            fh.write(content)
        removed = dedupe_psf_topology(path)
        self.assertEqual(removed, {})
        self.assertEqual(open(path).read(), content)   # byte-for-byte unchanged

    def test_real_ext_psf_fixture_is_clean(self):
        fixture = os.path.join(os.path.dirname(__file__), '..', '..', 'inputs', 'existing.psf')
        if not os.path.exists(fixture):
            self.skipTest('existing.psf fixture not present')
        path = os.path.join(self.tmp, 'existing.psf')
        shutil.copy(fixture, path)
        before = open(path).read()
        removed = dedupe_psf_topology(path)
        self.assertEqual(removed, {})                  # a well-formed PSF has nothing to remove
        self.assertEqual(open(path).read(), before)

    def test_nonexistent_and_non_psf(self):
        self.assertEqual(dedupe_psf_topology(os.path.join(self.tmp, 'nope.psf')), {})
        junk = os.path.join(self.tmp, 'junk.txt')
        with open(junk, 'w') as fh:
            fh.write('not a psf\n')
        self.assertEqual(dedupe_psf_topology(junk), {})


if __name__ == '__main__':
    unittest.main()
