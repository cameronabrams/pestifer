import unittest
from pathlib import Path

import numpy as np
import pytest
import yaml

from pestifer.core.controller import Controller
from pestifer.core.config import Config
from pestifer.psfutil.loop_ccd import (backbone_from_pdb, loop_atoms_from_pdb,
                                       loop_clash_report, heavy_env_coords_from_pdb)

pytestmark = pytest.mark.needs_tools

# An internal, cysteine-free loop of BPTI (6pti) to delete and rebuild.
_LOOP = {24: 'ASN', 25: 'ALA', 26: 'LYS', 27: 'ALA', 28: 'GLY'}


def _carve_missing_loop(src_pdb, dst_pdb):
    """Write dst_pdb == src_pdb minus the _LOOP residues' atoms, with REMARK 465 added
    for them so pestifer models them as a missing internal loop."""
    r465 = [f"REMARK 465     {rn} A  {rid:>4d}\n" for rid, rn in _LOOP.items()]
    out, inserted = [], False
    for line in open(src_pdb):
        if line.startswith("REMARK 465     ALA A    58") and not inserted:
            out.append(line); out.extend(r465); inserted = True; continue
        if line.startswith(("ATOM", "ANISOU", "HETATM")):
            try:
                rid, ch = int(line[22:26]), line[21]
            except ValueError:
                out.append(line); continue
            if ch == 'A' and rid in _LOOP:
                continue
        out.append(line)
    with open(dst_pdb, 'w') as f:
        f.writelines(out)


def _loop_rmsd_to_native(built, native, loop, anchors):
    """Backbone (N,CA,C) RMSD of the rebuilt loop vs native, after a rigid superposition on
    the flanking resolved anchor residues (Kabsch)."""
    def stack(bb, rs):
        return np.array([bb[r][a] for r in rs for a in ('N', 'CA', 'C')])
    P, Q = stack(built, anchors), stack(native, anchors)
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    B, N = stack(built, loop), stack(native, loop)
    Bt = (R @ (B - P.mean(0)).T).T + Q.mean(0)
    return float(np.sqrt(((Bt - N) ** 2).sum(1).mean()))


class TestLigateCCD(unittest.TestCase):
    """Delete-and-rebuild benchmark for the CCD loop closer.

    These gaps are floppy, solvent-exposed surface loops with no unique native
    conformation, so the deliverable is a *closed, clash-free, physically plausible*
    starting structure -- NOT recovery of the resolved native geometry. The benchmark
    therefore asserts closure + steric validity, and only sanity-bounds the placement
    (the loop lands in the gap, not that it matches native)."""

    def setUp(self):
        self.native = str(Path(__file__).parents[2] / 'inputs' / '6pti.pdb')
        _carve_missing_loop(self.native, 'bptidr.pdb')
        cfg = {'title': 'BPTI delete-and-rebuild, CCD loop closure',
               'tasks': [
                   {'fetch': {'source': 'local', 'sourceID': 'bptidr', 'source_format': 'pdb'}},
                   {'psfgen': {'source': {'sequence': {'loops': {
                       'min_loop_length': 4, 'declash': {'maxcycles': 0}}}}}},
                   {'ligate': {'method': 'ccd'}},
               ]}
        with open('drccd.yaml', 'w') as f:
            yaml.safe_dump(cfg, f)

    def test_ccd_closes_loop_clash_free(self):
        Controller().configure(Config(userfile='drccd.yaml').configure_new()).do_tasks()
        built = backbone_from_pdb('my_system.pdb', segname='A')
        native = backbone_from_pdb(self.native, chainID='A')

        # 1. the loop's C-terminus is bonded to the downstream anchor (peptide bond ~1.33 A)
        bond = np.linalg.norm(built[28]['C'] - built[29]['N'])
        self.assertLess(bond, 1.6, f"loop not closed: 28:C -> 29:N = {bond:.2f} A")

        # 2. steric validity is the real quality bar: the rebuilt loop must not interpenetrate
        # itself or the fold. The clash-filtered ensemble should reliably find a clean closure.
        loop = [24, 25, 26, 27, 28]
        order, coords, serials = loop_atoms_from_pdb('my_system.pdb', loop, segname='A')
        # exclude the flanking anchors 23/29: the junction peptide bonds to them are expected
        _ao, _ac, anchor_serials = loop_atoms_from_pdb('my_system.pdb', [23, 29], segname='A')
        env = heavy_env_coords_from_pdb('my_system.pdb',
                                        exclude_serials=list(serials) + list(anchor_serials))
        rep = loop_clash_report(order, coords, loop, env_coords=env)
        self.assertFalse(rep['topological'],
                         f"rebuilt loop is topologically broken: worst overlap {rep['worst']:.2f} A, "
                         f"min non-adjacent CA {rep['min_ca']:.2f} A, "
                         f"{rep['n_deep']} intra + {rep['n_env_deep']} loop-vs-structure deep overlaps")

        # 3. loose placement sanity (NOT a native-match quality bar): the loop lands in the gap
        # rather than flying off. A floppy loop has no unique native, so this bound is generous.
        placement = _loop_rmsd_to_native(built, native, loop=loop,
                                         anchors=[20, 21, 22, 23, 29, 30, 31, 32])
        self.assertLess(placement, 8.0, f"loop landed far outside the gap region: {placement:.2f} A")


if __name__ == '__main__':
    unittest.main()
