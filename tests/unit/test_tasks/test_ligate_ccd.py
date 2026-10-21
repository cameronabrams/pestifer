import unittest
from unittest import mock
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



class TestLoopPiercingResolution(unittest.TestCase):
    """Closing a loop is an unphysical displacement, so it owes a piercing check.

    CCD drives backbone dihedrals to meet a target and `steer` pulls termini together with a
    colvars bias; both supply the energy to cross a van der Waals barrier that free dynamics at
    a kT of about 0.6 cannot.  And a threading has no atomic clash, so the clash-filtered
    closure ensemble is blind to it -- the same reason psfgen needs `resolve_tail_piercings`
    beside its declash.  Nine shipped examples run `ligate` and then `md` with nothing between.
    """

    def _task(self, **specs):
        from pestifer.tasks.ligate import LigateTask
        t = LigateTask.__new__(LigateTask)
        t.taskname, t.basename, t.specs, t.outcome = 'ligate', 'test', dict(specs), {}
        t.next_basename = mock.Mock()
        t.register = mock.Mock()
        state = mock.Mock()
        state.psf.name, state.pdb.name = 'in.psf', 'in.pdb'
        state.xsc = None
        t.get_current_artifact = lambda k: state
        t.base_molecule = mock.Mock()
        t.base_molecule.protein_loop_gaps.return_value = [
            {'segname': 'A', 'loop_resids': ['10', '11', '12']}]
        t.get_current_artifact_data = lambda k: t.base_molecule
        t.log_message = mock.Mock()
        return t

    def test_the_loop_serials_reach_the_engine(self):
        t = self._task()
        with mock.patch('pestifer.psfutil.loop_ccd.loop_atoms_from_pdb',
                        return_value=(None, None, [7, 8, 9])), \
             mock.patch('pestifer.psfutil.ring_resolve.resolve_model_built_piercings',
                        return_value=(None, [])) as engine:
            t.resolve_loop_piercings()
        engine.assert_called_once()
        self.assertEqual(engine.call_args.args[4], [7, 8, 9],
                         'the closed loop was not handed to the resolver as the mover set')
        self.assertEqual(engine.call_args.kwargs['mover_kind'], 'loop')

    def test_serials_come_from_the_current_state_pdb(self):
        t = self._task()
        with mock.patch('pestifer.psfutil.loop_ccd.loop_atoms_from_pdb',
                        return_value=(None, None, [1])) as lookup, \
             mock.patch('pestifer.psfutil.ring_resolve.resolve_model_built_piercings',
                        return_value=(None, [])):
            t.resolve_loop_piercings()
        self.assertEqual(lookup.call_args.args[0], 'in.pdb')

    def test_check_piercings_false_skips_it(self):
        from pestifer.tasks.ligate import LigateTask
        t = self._task(check_piercings=False)
        with mock.patch.object(LigateTask, 'resolve_loop_piercings') as r, \
             mock.patch.object(LigateTask, 'connect', return_value=0), \
             mock.patch.object(LigateTask, 'close_loops_ccd', return_value=0):
            t.base_molecule.has_protein_loops = True
            LigateTask.do(t)
        r.assert_not_called()

    def test_the_default_runs_it(self):
        """Negative control for the test above, aimed at the CALL SITE: without this, a
        `resolve_loop_piercings` never invoked from do() would satisfy it."""
        from pestifer.tasks.ligate import LigateTask
        t = self._task()
        with mock.patch.object(LigateTask, 'resolve_loop_piercings') as r, \
             mock.patch.object(LigateTask, 'connect', return_value=0), \
             mock.patch.object(LigateTask, 'close_loops_ccd', return_value=0):
            t.base_molecule.has_protein_loops = True
            LigateTask.do(t)
        r.assert_called_once()

    def test_the_resolution_runs_BEFORE_connect(self):
        """Ordering is the whole feature.

        While the loop is unpatched its C-terminus is free, so the branch carrying a threading
        bond is a small pendant on a backbone bridge -- rotatable.  `connect` joins the loop at
        BOTH ends, and from then on a hinge inside it would swing the entire rest of the chain,
        so `RingChecker.pendant_axes` skips every axis (`max_pendant` is 180 atoms).

        Measured on example 11 seed 27021974 with the call placed after `connect` (as 3.29.1
        shipped it): the resolver found `HSD D-53 <- ALA D-50`, logged "no rotatable loop hinge
        found for ALA D-50", warned, and the build ran ~30 more minutes before `terminate`
        refused to package it.  A resolver that cannot succeed is worse than none -- it reports
        the problem in a voice that sounds like it handled it.
        """
        from pestifer.tasks.ligate import LigateTask
        t = self._task()
        calls = []
        with mock.patch.object(LigateTask, 'resolve_loop_piercings',
                               side_effect=lambda: calls.append('resolve')), \
             mock.patch.object(LigateTask, 'connect',
                               side_effect=lambda *a, **k: (calls.append('connect'), 0)[1]), \
             mock.patch.object(LigateTask, 'close_loops_ccd', return_value=0):
            t.base_molecule.has_protein_loops = True
            LigateTask.do(t)
        self.assertEqual(calls, ['resolve', 'connect'],
                         'the piercing resolution must run while the loop is still a pendant')


if __name__ == '__main__':
    unittest.main()
