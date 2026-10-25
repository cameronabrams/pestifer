# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""PDB atom serials past 99999 are hexadecimal, and every reader that maps a serial to an atom
must decode them in order.

VMD and psfgen write 100000 as ``186a0``.  ``int(line[6:11])`` read ``18700`` (atom 100096) as
atom 18700 and raised on ``186a0``, so every coordinate write-back through
`pdb_replace_coords` -- each rotation, spin and slide the piercing resolver makes, tail
modeling, CCD loop closure -- landed some atoms on other atoms' coordinates in any system past
99,999 atoms.  Measured 2026-10-25 on a 138,488-atom example-17 package: an identity rewrite
moved 3,300 atoms by up to 142 A.

The synthetic file below is written exactly as VMD writes one, so these tests need no tools.
"""
import os
import tempfile
import unittest

import numpy as np

from pestifer.util.coord import PDBSerialDecoder, pdb_replace_coords


def _vmd_serial(n):
    if n <= 99999:
        return f'{n:5d}'
    if n <= 0xfffff:
        return f'{n:5x}'
    return '*****'


def _write_pdb(path, n, first=1):
    with open(path, 'w') as fh:
        for i, n_ser in enumerate(range(first, first + n)):
            # columns as VMD writes them: name 13-16, resName 18-21, chain 22, resSeq 23-26
            fh.write(f'ATOM  {_vmd_serial(n_ser)}  OH2 TIP3W{(1 + i // 3) % 10000:4d}    '
                     f'{i * 0.001:8.3f}{i * 0.002:8.3f}{-i * 0.001:8.3f}  1.00  0.00      WT1  O\n')
        fh.write('END\n')


def _coords(path):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in open(path) if l.startswith('ATOM')])


class TestDecoder(unittest.TestCase):

    def test_decimal_then_hex(self):
        dec = PDBSerialDecoder()
        got = [dec(f'ATOM  {_vmd_serial(n)}') for n in (99998, 99999, 100000, 100096, 100106)]
        self.assertEqual(got, [99998, 99999, 100000, 100096, 100106])

    def test_an_all_digit_hex_serial_is_not_read_as_decimal(self):
        dec = PDBSerialDecoder()
        dec('ATOM  99999')
        self.assertEqual(dec('ATOM  18700'), 0x18700)

    def test_overflow_stars_continue_the_count(self):
        dec = PDBSerialDecoder()
        dec(f'ATOM  {_vmd_serial(0xfffff)}')
        self.assertEqual(dec('ATOM  *****'), 0xfffff + 1)
        self.assertEqual(dec('ATOM  *****'), 0xfffff + 2)


class TestRewriteAcross99999(unittest.TestCase):
    N = 100200    # crosses 99999, and 100096 = 0x18700 is all digits

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.src = os.path.join(cls.tmp, 'big.pdb')
        _write_pdb(cls.src, cls.N)
        cls.xyz = _coords(cls.src)
        cls.row_of = {s: s - 1 for s in range(1, cls.N + 1)}

    def test_an_identity_rewrite_is_byte_identical(self):
        out = os.path.join(self.tmp, 'same.pdb')
        pdb_replace_coords(self.src, out, self.xyz, self.row_of)
        self.assertEqual(open(out).read(), open(self.src).read())

    def test_moving_one_atom_past_99999_moves_only_that_atom(self):
        out = os.path.join(self.tmp, 'one.pdb')
        moved = self.xyz.copy()
        moved[100096 - 1] += 5.0
        pdb_replace_coords(self.src, out, moved, self.row_of)
        changed = np.nonzero(np.abs(_coords(out) - self.xyz).sum(axis=1) > 1e-6)[0]
        self.assertEqual(changed.tolist(), [100096 - 1])

    def test_loop_atoms_from_pdb_reports_the_true_serial(self):
        from pestifer.psfutil.loop_ccd import loop_atoms_from_pdb
        # atoms 100096-100098 are residue 1 + 100095 // 3, written modulo 10000 as VMD does
        resid = (1 + (100096 - 1) // 3) % 10000
        _order, _xyz, serials = loop_atoms_from_pdb(self.src, [resid], segname='WT1')
        self.assertIn(100096, serials)


if __name__ == '__main__':
    unittest.main()
