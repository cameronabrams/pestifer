import os
import shutil
import tempfile
import types
import unittest

import numpy as np

from pestifer.charmmff.make_solvent_box import (
    box_edge_for_density, cubic_lattice_sites, pack_cubic,
    write_box_pdb, atom_lines_of, coords_of,
    _resi_to_molblock, single_molecule_from_graph,
)


def _fake_meoh_topo():
    """A minimal stand-in for a CharmmResi: methanol's atom+bond graph, no coordinates."""
    A = lambda name, element: types.SimpleNamespace(name=name, element=element)
    B = lambda n1, n2, d=1: types.SimpleNamespace(name1=n1, name2=n2, degree=d)
    return types.SimpleNamespace(
        resname='MEOH',
        atoms=[A('CB', 'C'), A('OG', 'O'), A('HG1', 'H'), A('HB1', 'H'), A('HB2', 'H'), A('HB3', 'H')],
        bonds=[B('CB', 'OG'), B('OG', 'HG1'), B('CB', 'HB1'), B('CB', 'HB2'), B('CB', 'HB3')],
    )


class TestSolventBoxGeometry(unittest.TestCase):
    def test_edge_matches_known_water_box(self):
        # 216 TIP3P waters at 1.0 g/cc -> ~18.6 A (VMD's built-in 216-water box is ~18.77)
        e = box_edge_for_density(216, 18.015, 1.0)
        self.assertAlmostEqual(e, 18.63, places=1)

    def test_edge_scaling_and_validation(self):
        # 8x the molecules at the same density -> 2x the edge
        e1 = box_edge_for_density(216, 18.0, 1.0)
        e8 = box_edge_for_density(216 * 8, 18.0, 1.0)
        self.assertAlmostEqual(e8 / e1, 2.0, places=6)
        for bad in [(0, 18, 1), (216, 0, 1), (216, 18, 0)]:
            with self.assertRaises(ValueError):
                box_edge_for_density(*bad)

    def test_lattice_sites_fill_and_bound(self):
        e = 20.0
        s = cubic_lattice_sites(216, e)         # perfect cube -> 6x6x6
        self.assertEqual(s.shape, (216, 3))
        self.assertTrue((s >= 0).all() and (s < e).all())
        self.assertEqual(cubic_lattice_sites(100, e).shape, (100, 3))   # non-cube -> next cube, truncated
        self.assertEqual(cubic_lattice_sites(1, e).shape, (1, 3))

    def test_pack_is_rigid_and_deterministic(self):
        mol = np.array([[0.0, 0, 0], [1.5, 0, 0], [0, 1.0, 0]])
        e = box_edge_for_density(64, 40.0, 0.8)
        p = pack_cubic(mol, 64, e, seed=7)
        self.assertEqual(p.shape, (64, 3, 3))
        # random rotation preserves internal geometry for every copy
        for i in range(64):
            self.assertAlmostEqual(np.linalg.norm(p[i, 1] - p[i, 0]), 1.5, places=5)
            self.assertAlmostEqual(np.linalg.norm(p[i, 2] - p[i, 0]), 1.0, places=5)
        self.assertTrue(np.allclose(p, pack_cubic(mol, 64, e, seed=7)))   # deterministic

    def test_write_box_pdb_roundtrips(self):
        template = [
            "ATOM      1  O   GRL X   1       0.000   0.000   0.000  1.00  0.00      X    O\n",
            "ATOM      2  C1  GRL X   1       1.400   0.000   0.000  1.00  0.00      X    C\n",
        ]
        mol = coords_of(template)
        self.assertEqual(mol.shape, (2, 3))
        placed = pack_cubic(mol, 8, box_edge_for_density(8, 32.0, 0.79), seed=0)
        d = tempfile.mkdtemp()
        out = os.path.join(d, 'box.pdb')
        # 8 molecules become one segment 'QQQ' of 8 residues -- VMD solvate's custom-box path
        # hardwires that segid, so a box with any other segid fails to tile
        segments = write_box_pdb(placed, template, out)
        self.assertEqual(segments, [('QQQ', 8)])
        lines = atom_lines_of(out)
        self.assertEqual(len(lines), 16)                       # 8 molecules x 2 atoms
        self.assertEqual(lines[0][72:76].strip(), 'QQQ')       # segid written in the right column
        self.assertEqual(lines[0][21], 'A')                    # chainID stamped (in manager's pool)
        self.assertEqual(lines[0][22:26].strip(), '1')         # molecule 1 -> resid 1
        self.assertEqual(lines[2][22:26].strip(), '2')         # molecule 2 -> resid 2
        # coordinates round-trip back to the placed positions
        self.assertTrue(np.allclose(coords_of(lines[:2]), placed[0], atol=1e-3))
        # more than 9999 molecules cannot fit one segment's PDB resid field
        with self.assertRaises(ValueError):
            write_box_pdb(np.zeros((10001, 2, 3)), template, os.path.join(d, 'big.pdb'))


class TestGraphToMolecule(unittest.TestCase):
    def test_resi_to_molblock(self):
        block = _resi_to_molblock(_fake_meoh_topo())
        lines = block.splitlines()
        # title, program, comment, counts line, then 6 atoms + 5 bonds + M END
        self.assertEqual(lines[0], 'MEOH')
        self.assertTrue(lines[3].startswith('  6  5'))          # 6 atoms, 5 bonds
        self.assertEqual(lines[3].rstrip().endswith('V2000'), True)
        atom_lines = lines[4:10]
        self.assertEqual(len(atom_lines), 6)
        self.assertTrue(atom_lines[0].strip().endswith('C'.ljust(1)) or ' C ' in atom_lines[0])
        # coordinates are all zero (obabel --gen3d fills them in)
        for al in atom_lines:
            self.assertIn('0.0000', al)
        # the V2000 valence field (vvv, cols 49-51) is pinned to each atom's bond-order sum,
        # so obabel adds no implicit H (CB has 4 bonds, OG has 2)
        self.assertEqual(int(atom_lines[0][48:51]), 4)
        self.assertEqual(int(atom_lines[1][48:51]), 2)
        # bond block references 1-based atom indices in RESI order (CB=1, OG=2, HG1=3, ...)
        bond_lines = lines[10:15]
        self.assertEqual(bond_lines[0].split()[:3], ['1', '2', '1'])   # CB-OG single
        self.assertEqual(bond_lines[1].split()[:3], ['2', '3', '1'])   # OG-HG1 single
        self.assertEqual(lines[-1], 'M  END')

    @unittest.skipUnless(shutil.which('obabel'), 'obabel not on PATH')
    def test_single_molecule_from_graph_gen3d(self):
        d = tempfile.mkdtemp()
        cwd = os.getcwd()
        os.chdir(d)
        try:
            topo = _fake_meoh_topo()
            out = single_molecule_from_graph(topo, 'MEOH-seed.pdb')
            self.assertTrue(os.path.isfile(out))
            lines = atom_lines_of(out)
            # one atom per RESI atom, in RESI order, stamped with the CHARMM names
            self.assertEqual(len(lines), 6)
            names = [ln[12:16].strip() for ln in lines]
            self.assertEqual(names, ['CB', 'OG', 'HG1', 'HB1', 'HB2', 'HB3'])
            self.assertEqual(lines[0][17:21].strip(), 'MEOH')   # resName cols 18-21
            # obabel produced real 3D coordinates (not all at the origin)
            xyz = coords_of(lines)
            self.assertGreater(np.ptp(xyz, axis=0).max(), 0.5)
            # bonded C-O distance is chemically sane (~1.4 A)
            self.assertLess(np.linalg.norm(xyz[0] - xyz[1]), 1.8)
        finally:
            os.chdir(cwd)


if __name__ == '__main__':
    unittest.main()
