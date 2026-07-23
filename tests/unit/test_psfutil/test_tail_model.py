import unittest
from types import SimpleNamespace as NS

import numpy as np

from pestifer.psfutil.loop_ccd import place_atom_nerf
from pestifer.psfutil.tail_model import model_one_tail, downstream_anchor_target, _selection_key
from pestifer.molecule.stateinterval import StateInterval, StateIntervalList
from pestifer.molecule.molecule import Molecule

# Ideal internal coordinates for building a synthetic all-ALA test peptide.
_IC = dict(C_N=1.33, N_CA=1.45, CA_C=1.52, C_O=1.23, CA_CB=1.53,
           CA_C_N=116.6, C_N_CA=121.7, N_CA_C=111.0, N_C_O=122.7, N_CA_CB=110.5)


def _build_peptide(nres, phi=-120.0, psi=140.0):
    """A synthetic all-ALA peptide (N, CA, C, O, CB per residue) grown by NeRF; resids 1..nres."""
    atoms, prevN, prevCA, prevC = [], None, None, None
    for r in range(1, nres + 1):
        if r == 1:
            N = np.array([0.0, 0.0, 0.0]); CA = np.array([1.45, 0.0, 0.0])
            C = place_atom_nerf(np.array([-1.0, 1.0, 0.0]), N, CA, _IC['CA_C'], _IC['N_CA_C'], phi)
        else:
            N = place_atom_nerf(prevN, prevCA, prevC, _IC['C_N'], _IC['CA_C_N'], psi)
            CA = place_atom_nerf(prevCA, prevC, N, _IC['N_CA'], _IC['C_N_CA'], 180.0)
            C = place_atom_nerf(prevC, N, CA, _IC['CA_C'], _IC['N_CA_C'], phi)
        O = place_atom_nerf(N, CA, C, _IC['C_O'], _IC['N_C_O'], 0.0)
        CB = place_atom_nerf(N, C, CA, _IC['CA_CB'], _IC['N_CA_CB'], 120.0)
        atoms += [(r, 'N', N), (r, 'CA', CA), (r, 'C', C), (r, 'O', O), (r, 'CB', CB)]
        prevN, prevCA, prevC = N, CA, C
    return atoms


def _write_pdb(atoms, path, seg='A'):
    with open(path, 'w') as f:
        for i, (r, n, x) in enumerate(atoms, 1):
            f.write(f"ATOM  {i:>5d} {n:<4s} ALA {seg[0]}{r:>4d}    "
                    f"{x[0]:>8.3f}{x[1]:>8.3f}{x[2]:>8.3f}  1.00  0.00      {seg:<4s}\n")


def _atom_from_file(path, resid, name):
    with open(path) as f:
        for line in f:
            if line[22:26].strip() == str(resid) and line[12:16].strip() == name:
                return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return None


class TestDownstreamAnchorTarget(unittest.TestCase):
    """The C-tail junction target is the mirror of loop_ccd.anchor_closure_target."""

    def test_peptide_bond_geometry(self):
        A_N = np.array([0.0, 0.0, 0.0]); A_CA = np.array([1.46, 0.0, 0.0]); A_C = np.array([2.0, 1.42, 0.0])
        N, CA, C = downstream_anchor_target(A_N, A_CA, A_C)
        self.assertAlmostEqual(np.linalg.norm(N - A_C), 1.33, places=2)   # anchor C -> tail N
        self.assertAlmostEqual(np.linalg.norm(CA - N), 1.45, places=2)
        self.assertAlmostEqual(np.linalg.norm(C - CA), 1.52, places=2)


class TestModelOneTail(unittest.TestCase):

    def setUp(self):
        import tempfile, os
        self.d = tempfile.mkdtemp()
        self.pdb = os.path.join(self.d, 'pep.pdb')
        atoms = _build_peptide(8)
        _write_pdb(atoms, self.pdb)
        self.serials = {}
        for i, (r, _n, _x) in enumerate(atoms, 1):
            self.serials.setdefault(r, []).append(i)

    def _modeled(self, tail):
        return [s for r in tail for s in self.serials[r]]

    def test_c_tail_junction_and_atoms(self):
        # C-tail: anchor = res 1, tail = res 2..8; first tail residue N must bond anchor C (~1.33 A)
        tail = list(range(2, 9))
        res = model_one_tail(self.pdb, 'A', tail, 1, 'C', self._modeled(tail), seed=7, ensemble=6)
        self.assertEqual(len(res['serials']), 5 * len(tail))      # all atoms carried
        row = {k: i for i, k in enumerate(res['order'])}
        tailN = res['modeled'][row[(2, 'N')]]
        anchorC = _atom_from_file(self.pdb, 1, 'C')
        self.assertAlmostEqual(float(np.linalg.norm(tailN - anchorC)), 1.33, places=1)

    def test_n_tail_junction(self):
        # N-tail: anchor = res 8, tail = res 1..7; last tail residue C must bond anchor N (~1.33 A)
        tail = list(range(1, 8))
        res = model_one_tail(self.pdb, 'A', tail, 8, 'N', self._modeled(tail), seed=7, ensemble=6)
        row = {k: i for i, k in enumerate(res['order'])}
        tailC = res['modeled'][row[(7, 'C')]]
        anchorN = _atom_from_file(self.pdb, 8, 'N')
        self.assertAlmostEqual(float(np.linalg.norm(tailC - anchorN)), 1.33, places=1)

    def test_deterministic_for_fixed_seed(self):
        tail = list(range(2, 9))
        a = model_one_tail(self.pdb, 'A', tail, 1, 'C', self._modeled(tail), seed=7, ensemble=6)
        b = model_one_tail(self.pdb, 'A', tail, 1, 'C', self._modeled(tail), seed=7, ensemble=6)
        np.testing.assert_allclose(a['modeled'], b['modeled'])

    def test_bond_lengths_preserved(self):
        # a rigid-placed Ramachandran backbone keeps ideal peptide bond lengths within the tail
        tail = list(range(2, 9))
        res = model_one_tail(self.pdb, 'A', tail, 1, 'C', self._modeled(tail), seed=3, ensemble=6)
        row = {k: i for i, k in enumerate(res['order'])}
        X = res['modeled']
        for r in tail:
            self.assertAlmostEqual(float(np.linalg.norm(X[row[(r, 'N')]] - X[row[(r, 'CA')]])), 1.45, places=1)
            self.assertAlmostEqual(float(np.linalg.norm(X[row[(r, 'CA')]] - X[row[(r, 'C')]])), 1.52, places=1)


class TestSelectionKey(unittest.TestCase):
    """_selection_key selects by segid when present, else falls back to the chain column."""

    def _write(self, path, segid):
        with open(path, 'w') as f:
            f.write(f"ATOM      1 N    ALA A   1    "
                    f"{0.0:>8.3f}{0.0:>8.3f}{0.0:>8.3f}  1.00  0.00      {segid:<4s}\n")

    def test_prefers_segid_when_present(self):
        import tempfile, os
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'a.pdb'); self._write(p, 'A')
            self.assertEqual(_selection_key(p, 'A'), {'segname': 'A'})

    def test_falls_back_to_chain_when_segid_blank(self):
        import tempfile, os
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'a.pdb'); self._write(p, '')     # blank segid (standard dialect)
            self.assertEqual(_selection_key(p, 'A'), {'chainID': 'A'})


class TestProteinTerminalTails(unittest.TestCase):
    """Molecule.protein_terminal_tails enumeration logic, exercised on a lightweight stub."""

    @staticmethod
    def _stub(subsegs, resids, segname='A'):
        # segment with the given subsegments over residues numbered by `resids`
        residues = [NS(resid=NS(resid=n), resname='ALA') for n in resids]
        seg = NS(segtype='protein', subsegments=StateIntervalList(subsegs),
                 residues=residues, segname=segname)
        transform = NS(chainIDmap={})
        return NS(asymmetric_unit=NS(segments=[seg]),
                  active_biological_assembly=NS(transforms=[transform]))

    def test_enumerates_built_n_and_c_tails(self):
        # [MISSING built 1-2][RESOLVED 3-6][MISSING built 7-8]
        subsegs = [StateInterval(state='MISSING', bounds=[0, 1], build=True),
                   StateInterval(state='RESOLVED', bounds=[2, 5]),
                   StateInterval(state='MISSING', bounds=[6, 7], build=True)]
        stub = self._stub(subsegs, [1, 2, 3, 4, 5, 6, 7, 8])
        tails = Molecule.protein_terminal_tails(stub)
        by_end = {t['end']: t for t in tails}
        self.assertEqual(set(by_end), {'N', 'C'})
        self.assertEqual(by_end['N']['tail_resids'], [1, 2])
        self.assertEqual(by_end['N']['anchor_resid'], 3)      # first resolved residue
        self.assertEqual(by_end['C']['tail_resids'], [7, 8])
        self.assertEqual(by_end['C']['anchor_resid'], 6)      # last resolved residue

    def test_unbuilt_tail_excluded(self):
        # an unbuilt terminal MISSING run must not be enumerated
        subsegs = [StateInterval(state='MISSING', bounds=[0, 1], build=False),
                   StateInterval(state='RESOLVED', bounds=[2, 5])]
        stub = self._stub(subsegs, [1, 2, 3, 4, 5, 6])
        self.assertEqual(Molecule.protein_terminal_tails(stub), [])

    def test_interior_gap_excluded(self):
        # an interior MISSING run is handled by protein_loop_gaps, not here
        subsegs = [StateInterval(state='RESOLVED', bounds=[0, 1]),
                   StateInterval(state='MISSING', bounds=[2, 3], build=True),
                   StateInterval(state='RESOLVED', bounds=[4, 5])]
        stub = self._stub(subsegs, [1, 2, 3, 4, 5, 6])
        self.assertEqual(Molecule.protein_terminal_tails(stub), [])


if __name__ == '__main__':
    unittest.main()
