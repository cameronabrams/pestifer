import copy
import pytest
import os
import shutil
import unittest
from unittest import mock

import numpy as np
import pandas as pd

from pestifer.core.artifacts import DataArtifact, StateArtifacts
from pestifer.core.config import Config
from pestifer.core.controller import Controller

from pestifer.scripters import PsfgenScripter

from pestifer.scripters.vmd import VMDScripter
from pestifer.tasks.make_membrane_system import (
    MakeMembraneSystemTask, _AREA_CONVERGENCE_TOL, _per_leaflet_tension,
    _MIN_GRID_HALF_MID_ZGAP)

from pestifer.tasks.validate import ValidateTask
from pestifer.util.util import protect_str_arg

pytestmark = pytest.mark.needs_tools

class _StubMDPlot:
    """Minimal stand-in for the trailing mdplot task: just carries a ``dataframes`` dict."""
    def __init__(self, xst=None):
        self.dataframes = {} if xst is None else {'xst': xst}


def _stage_orient_oracle(vm):
    """Make the retired ``bilayer_orient.tcl`` runnable from the current working directory.

    The script and the third-party ``Orient``/``La`` packages it needs are test fixtures, not
    package resources -- no pestifer workflow loads them, and they are kept only so these tests
    can compare the Python orientation path against the Tcl path it replaced.  Because they no
    longer live under ``${PESTIFER_TCLROOT}/pkg/``, ``vmdrc.tcl`` does not glob them onto
    ``auto_path``, so we add that entry ourselves before the script is ingested.

    Call after ``vm.newscript(...)`` and before ``vm.usescript('bilayer_orient', local=True)``.
    """
    oracle = os.path.abspath(os.path.join('..', '..', 'fixtures', 'tcl_oracle'))
    shutil.copy(os.path.join(oracle, 'bilayer_orient.tcl'), '.')
    vm.addline(f'lappend auto_path {os.path.join(oracle, "pkg")}')


def _xst_from_area(area):
    """Build an xst-style dataframe whose a_x*b_y equals the given area trajectory."""
    area = np.asarray(area, dtype=float)
    return pd.DataFrame({'a_x': area, 'b_y': np.ones_like(area)})


class TestPhaseCompositionGuard(unittest.TestCase):
    """The phase/sterol mismatch warning (no NAMD required)."""

    _warn = staticmethod(MakeMembraneSystemTask._warn_phase_composition_mismatch)

    def _warns(self, leaflet, phase, species):
        with mock.patch('pestifer.tasks.make_membrane_system.logger.warning') as m:
            self._warn(leaflet, phase, species)
            return m.called

    def test_ld_on_cholesterol_rich_warns(self):
        raft = [{'name': 'PSM', 'frac': 0.36}, {'name': 'POPC', 'frac': 0.17},
                {'name': 'CHL1', 'frac': 0.47}]
        self.assertTrue(self._warns('upper_leaflet', 'Ld', raft))

    def test_lo_on_cholesterol_rich_is_fine(self):
        raft = [{'name': 'PSM', 'frac': 0.36}, {'name': 'CHL1', 'frac': 0.47}]
        self.assertFalse(self._warns('upper_leaflet', 'Lo', raft))

    def test_lo_on_sterol_poor_warns(self):
        fluid = [{'name': 'POPC', 'frac': 0.9}, {'name': 'CHL1', 'frac': 0.1}]
        self.assertTrue(self._warns('lower_leaflet', 'Lo', fluid))

    def test_ld_on_fluid_is_fine(self):
        fluid = [{'name': 'POPC', 'frac': 1.0}]
        self.assertFalse(self._warns('lower_leaflet', 'Ld', fluid))


class TestAreaConvergence(unittest.TestCase):
    """Unit tests for the calibration-cell area-convergence guard (no NAMD required)."""

    # _area_convergence does not touch self, so an unbound call with None is fine
    _check = staticmethod(lambda task, name: MakeMembraneSystemTask._area_convergence(None, task, name))

    def test_plateaued_area_is_converged(self):
        # a flat tail (already equilibrated) -> drift well under tolerance
        area = np.concatenate([np.linspace(8000, 6000, 60), np.full(60, 6000.0)])
        drift = self._check(_StubMDPlot(_xst_from_area(area)), 'flat')
        self.assertIsNotNone(drift)
        self.assertLess(abs(drift), _AREA_CONVERGENCE_TOL)

    def test_still_shrinking_area_flagged(self):
        # area decreasing monotonically right to the end -> drift exceeds tolerance, negative
        area = np.linspace(7000, 6000, 120)
        drift = self._check(_StubMDPlot(_xst_from_area(area)), 'shrinking')
        self.assertIsNotNone(drift)
        self.assertGreater(abs(drift), _AREA_CONVERGENCE_TOL)
        self.assertLess(drift, 0.0)

    def test_missing_trajectory_returns_none(self):
        self.assertIsNone(self._check(_StubMDPlot(None), 'nodata'))
        self.assertIsNone(self._check(None, 'notask'))

    def test_too_short_trajectory_returns_none(self):
        # only 10 frames -> tail of 5 < 8, cannot judge
        drift = self._check(_StubMDPlot(_xst_from_area(np.linspace(7000, 6000, 10))), 'short')
        self.assertIsNone(drift)


class TestRelaxedAreaDrift(unittest.TestCase):
    """Prefer a converged plateau-gated membrane_equilibrate's own stage-2 measure over the
    cruder whole-trajectory _area_convergence check -- but RECOMPUTE it, never read the stored
    verdict.

    This class previously contained ``test_trusts_converged_plateau_gated``, which asserted that
    a converged stage's stored ``_last_plateau_drift`` was handed straight back. That encoded the
    tautology as intended behavior: the gate declares convergence only when that number is below
    ``area_plateau_tol``, and the downstream calibration check asks whether it EXCEEDS the same
    tolerance -- so the check could never fire on a converged patch, and was blind to premature
    convergence, which is the dangerous case. The word "trusts" in the old name was accurate and
    was the problem. Renamed and inverted 2026-10-14.
    """

    def _mms(self):
        return MakeMembraneSystemTask.__new__(MakeMembraneSystemTask)

    @staticmethod
    def _me(plateau_tol, converged, drift, recomputed=None):
        from types import SimpleNamespace
        ns = SimpleNamespace(_area_plateau_tol=plateau_tol, converged=converged,
                             _last_plateau_drift=drift)
        ns._area_plateau_drift = lambda: (drift if recomputed is None else recomputed)
        return ns

    def test_uses_the_stage_2_measure_not_the_whole_trajectory(self):
        from types import SimpleNamespace
        # a converged plateau-gated ME -> its own stage-2 drift, not the still-descending whole run
        tasks = [SimpleNamespace(), self._me(0.005, True, 0.0012)]
        last = _StubMDPlot(_xst_from_area(np.linspace(7000, 6000, 120)))
        self.assertEqual(self._mms()._relaxed_area_drift(tasks, last, 'patchA'), 0.0012)

    def test_it_RECOMPUTES_rather_than_reading_the_stored_verdict(self):
        """The regression. If the stored number and a fresh computation disagree, the fresh one
        must win -- otherwise the audit is quoting the verdict it audits."""
        from types import SimpleNamespace
        tasks = [SimpleNamespace(), self._me(0.005, True, 0.0012, recomputed=-0.0400)]
        last = _StubMDPlot(_xst_from_area(np.linspace(7000, 6000, 120)))
        self.assertEqual(self._mms()._relaxed_area_drift(tasks, last, 'patchA'), -0.0400)

    def test_falls_back_when_ceilinged(self):
        # a plateau-gated ME that hit its ceiling (converged False) -> whole-trajectory _area_convergence
        tasks = [self._me(0.005, False, 0.0012)]
        last = _StubMDPlot(_xst_from_area(np.linspace(7000, 6000, 120)))
        drift = self._mms()._relaxed_area_drift(tasks, last, 'patchB')
        self.assertIsNotNone(drift)
        self.assertLess(drift, 0.0)                       # sees the descent, would warn (as intended)

    def test_falls_back_for_md_terminated(self):
        from types import SimpleNamespace
        # no plateau-gated task (legacy fixed md ladder) -> _area_convergence on the flat tail
        tasks = [SimpleNamespace(), SimpleNamespace()]
        last = _StubMDPlot(_xst_from_area(np.concatenate([np.linspace(8000, 6000, 60),
                                                          np.full(60, 6000.0)])))
        drift = self._mms()._relaxed_area_drift(tasks, last, 'quilt')
        self.assertIsNotNone(drift)
        self.assertLess(abs(drift), _AREA_CONVERGENCE_TOL)


def _pp_df_from_dp(dp, nframes=4, c_z_unused=None):
    """Build a pressureprofile-style dataframe whose Pzz-(Pxx+Pyy)/2 equals dp per slab.

    Pxx=Pyy=0, Pzz=dp, repeated over nframes dynamics frames plus a junk TS=0 frame
    (which the parser must discard)."""
    nslabs = len(dp)
    cols = {'TS': [0] + [1000 * (k + 1) for k in range(nframes)]}
    for i in range(nslabs):
        cols[f'x_{i}'] = [0.0] * (nframes + 1)
        cols[f'y_{i}'] = [0.0] * (nframes + 1)
        # TS=0 frame gets a wildly wrong value to prove it is dropped
        cols[f'z_{i}'] = [1e6] + [float(dp[i])] * nframes
    return pd.DataFrame(cols)


class TestPerLeafletTension(unittest.TestCase):
    """Unit tests for the differential-stress (per-leaflet tension) parser (no NAMD)."""

    def test_symmetric_profile_zero_dgamma(self):
        # profile symmetric about the midplane -> upper and lower tensions equal
        nslabs = 40
        idx = np.arange(nslabs)
        dp = (idx - (nslabs - 1) / 2.0) ** 2          # symmetric about the central boundary
        res = _per_leaflet_tension(_pp_df_from_dp(dp), c_z=80.0)
        self.assertIsNotNone(res)
        self.assertEqual(res['nslabs'], nslabs)
        self.assertAlmostEqual(res['upper'], res['lower'], places=6)
        self.assertAlmostEqual(res['dgamma'], 0.0, places=6)

    def test_asymmetric_profile_signed_dgamma(self):
        # make the upper half more tensile than the lower -> nonzero, positive dgamma
        nslabs = 40
        dp = np.ones(nslabs)
        dp[nslabs // 2:] *= 3.0
        res = _per_leaflet_tension(_pp_df_from_dp(dp), c_z=80.0)
        self.assertIsNotNone(res)
        self.assertGreater(res['dgamma'], 0.0)
        # quantitative: gamma = 0.01 * dp * dz, dz = 80/40 = 2
        dz = 80.0 / nslabs
        self.assertAlmostEqual(res['lower'], 0.01 * (nslabs // 2) * 1.0 * dz, places=6)
        self.assertAlmostEqual(res['upper'], 0.01 * (nslabs // 2) * 3.0 * dz, places=6)

    def test_drops_pre_dynamics_frame(self):
        # the TS=0 frame holds 1e6; if it were not dropped the tensions would explode
        dp = np.ones(20)
        res = _per_leaflet_tension(_pp_df_from_dp(dp, nframes=6), c_z=60.0)
        self.assertIsNotNone(res)
        self.assertLess(abs(res['total']), 10.0)      # sane magnitude, not ~1e6

    def test_degenerate_inputs_return_none(self):
        self.assertIsNone(_per_leaflet_tension(None, c_z=80.0))
        self.assertIsNone(_per_leaflet_tension(pd.DataFrame(), c_z=80.0))
        self.assertIsNone(_per_leaflet_tension(pd.DataFrame({'TS': [1000]}), c_z=80.0))


class TestGridHalfMidZgap(unittest.TestCase):
    """Unit tests for the grid-placement mid-plane-gap clamp (guards the segfault)."""

    # _grid_half_mid_zgap only logs + returns; call unbound with a throwaway self
    _clamp = staticmethod(lambda v: MakeMembraneSystemTask._grid_half_mid_zgap(None, v))

    def test_zero_gap_clamped_up(self):
        self.assertEqual(self._clamp(0.0), _MIN_GRID_HALF_MID_ZGAP)

    def test_small_gap_clamped_up(self):
        self.assertEqual(self._clamp(0.25), _MIN_GRID_HALF_MID_ZGAP)

    def test_adequate_gap_preserved(self):
        self.assertEqual(self._clamp(1.5), 1.5)
        self.assertEqual(self._clamp(_MIN_GRID_HALF_MID_ZGAP), _MIN_GRID_HALF_MID_ZGAP)


class TestNpatch(unittest.TestCase):
    """Normalization of the npatch tiling; guards an IndexError on a non-embedding
    grid build when npatch is absent or the schema's empty-list default."""

    @staticmethod
    def _npatch(specs):
        from types import SimpleNamespace
        return MakeMembraneSystemTask._npatch(SimpleNamespace(bilayer_specs=specs))

    def test_unset_defaults_to_1x1(self):
        self.assertEqual(self._npatch({}), [1, 1])

    def test_empty_list_defaults_to_1x1(self):
        self.assertEqual(self._npatch({'npatch': []}), [1, 1])

    def test_none_defaults_to_1x1(self):
        self.assertEqual(self._npatch({'npatch': None}), [1, 1])

    def test_malformed_defaults_to_1x1(self):
        self.assertEqual(self._npatch({'npatch': [5]}), [1, 1])

    def test_valid_preserved(self):
        self.assertEqual(self._npatch({'npatch': [2, 3]}), [2, 3])


class TestOrientationAlign(unittest.TestCase):
    """Unit tests for MakeMembraneSystemTask._orientation_align (no VMD)."""

    def _align(self, embed_specs):
        task = MakeMembraneSystemTask.__new__(MakeMembraneSystemTask)
        task.embed_specs = embed_specs
        return task._orientation_align()

    def test_none_when_no_orientation(self):
        self.assertIsNone(self._align({}))

    def test_z_head_tail_shorthand_maps_to_align_onto_z(self):
        a = self._align({'z_head_group': 'protein and resid 1',
                         'z_tail_group': 'protein and resid 100'})
        self.assertIsNotNone(a)
        self.assertEqual(a.movetype, 'ALIGN')
        # tail -> head axis carried onto +z
        self.assertEqual(a.source, ['protein and resid 100', 'protein and resid 1'])
        self.assertEqual(a.target, [0.0, 0.0, 1.0])

    def test_explicit_orient_spec_used(self):
        a = self._align({'orient': {'source': [1, 0, 0], 'target': [0, 0, 1]}})
        self.assertEqual(a.movetype, 'ALIGN')
        self.assertEqual(a.source, [1, 0, 0])
        self.assertEqual(a.target, [0, 0, 1])

    def test_orient_spec_takes_precedence_over_shorthand(self):
        a = self._align({'orient': {'source': ['a', 'b'], 'target': [0, 0, 1]},
                         'z_head_group': 'protein and resid 1',
                         'z_tail_group': 'protein and resid 100'})
        self.assertEqual(a.source, ['a', 'b'])

    def test_partial_shorthand_is_ignored(self):
        # only one of the pair given -> no orientation
        self.assertIsNone(self._align({'z_head_group': 'protein and resid 1'}))


class TestMakeMembraneSystem(unittest.TestCase):

    # A CONTROLLER PER TEST, not per class.  `Bilayer.__init__` calls
    # `charmmffcontent.provision(resnames=...)`, and provisioning with an explicit list REPLACES
    # what the content object exposes rather than adding to it -- right for a build, which has one
    # composition, and wrong for a shared fixture.  With one controller for the class, the first
    # membrane test to run narrowed its CHARMM content to its own lipids and the next one died
    # with `Cannot find residue POPC in CHARMMFFContent`.  Ordering decides which test that is
    # (unittest runs them alphabetically), so the failure moves whenever a test is added or
    # renamed; it surfaced here only once the two `5e8w_*` tests, which sort first, stopped
    # failing early on a dead `z_dist` key and started provisioning for PSM/POPE/CHL1.
    def setUp(self):
        self.controller = Controller().configure(Config().configure_new(), terminate=False)
        self.scripters = self.controller.config.scripters # shortcut

    @classmethod
    def setUpClass(cls):
        cls._patch_relaxation_protocols = [
            {'md': {'ensemble': 'minimize', 'nsteps': 1000}},
            {'md': {'ensemble': 'NVT', 'nsteps': 1000}},
            {'md': {'ensemble': 'NPT', 'nsteps': 1000}},
            # {'md': {'ensemble': 'NPT', 'nsteps': 2000}},
            # {'md': {'ensemble': 'NPT', 'nsteps': 4000}}
        ]
        cls._quilt_relaxation_protocols = [
            {'md': {'ensemble': 'minimize', 'nsteps': 1000}},
            {'md': {'ensemble': 'NVT', 'nsteps': 1000}},
            {'md': {'ensemble': 'NPT', 'nsteps': 1000}}
        ]

    # Handed to five tests, and ycleptic's `dwalk` fills in every schema default IN PLACE on the
    # dict it is given -- so the class-level list above comes back from the first test with
    # `temperature`, `pressure`, `addl_paramfiles` and `other_parameters` added and `NVT`/`NPT`
    # lower-cased, and MakeMembraneSystemTask then restructures it further (it already deepcopies
    # the protocol between patchA and patchB for exactly this reason).  Shared, that mutation
    # reaches the next test as config it never wrote: the four membrane tests after the first
    # died in `reconfigure_tasks` with "Attribute 'minimize' invalid ... under 'md'".
    #
    # Invisible until today because these are `slow` AND `needs_tools`: CI has neither VMD nor
    # NAMD so it skips the module, and a developer run without `--runslow` deselects it.  Each
    # test gets its own copy.
    @property
    def common_patch_relaxation_protocols(self):
        return copy.deepcopy(self._patch_relaxation_protocols)

    @property
    def common_quilt_relaxation_protocols(self):
        return copy.deepcopy(self._quilt_relaxation_protocols)

    def _residue_counts(self, psf_artifact):
        """``{resname: residue count}`` for a built system."""
        from collections import Counter
        from pestifer.psfutil.psfcontents import PSFContents
        return Counter(r.resname for s in PSFContents(str(psf_artifact.path)).segments
                       for r in s.residues)

    def _assert_leaflet_sizing(self, task, upper_fracs, lower_fracs, n_upper_requested):
        """The built system must hold the lipid counts the task DECIDED on.

        Do NOT assert either leaflet against a literal taken from ``patch_nlipids``.  Only the
        UPPER count is the user's request: ``make_membrane_system`` sizes the lower leaflet from
        the stress-free ratio ``n_upper * apl_upper / apl_lower`` (and, when embedding, sizes
        BOTH from the protein footprint), so ``patch_nlipids['lower']`` never reaches the quilt.
        These tests used to assert ``resname POPE <= 49`` against a 49/49 request, which held
        only while the two leaflets' APLs happened to round to the same count; POPC and POPE
        differ enough that the build correctly makes 50 POPE and the assertion correctly failed.

        What IS an invariant is that the membrane built matches the sizing the task recorded.
        ``upper_fracs``/``lower_fracs`` are the configured ``{resname: frac}`` of each leaflet,
        so a species in BOTH leaflets (cholesterol, usually) is expected at the sum of its two
        shares -- counting by resname alone cannot tell the two apart, since the quilt puts every
        lipid in one segment.

        Counts come back ``<=`` expected because ``ring_check`` deletes lipids it cannot unthread,
        and ``>=`` 90% of it because a leaflet that quietly lost a tenth of its lipids is a build
        failure rather than a rounding difference -- without that bound the ``<=`` passes on an
        empty leaflet.
        """
        counts = self._residue_counts(task.get_current_artifact('state').psf)
        n_upper = task.outcome.get('n_upper')
        n_lower = task.outcome.get('n_lower')
        self.assertIsNotNone(n_upper, 'task recorded no n_upper')
        self.assertIsNotNone(n_lower, 'task recorded no n_lower')
        # the upper leaflet IS the request, exactly; the lower one is derived, so the only thing
        # it can be checked against is the task's own recorded decision.  When EMBEDDING, neither
        # is the request -- both come from the protein footprint -- so the caller passes None.
        if n_upper_requested is not None:
            self.assertEqual(n_upper, n_upper_requested,
                             'upper leaflet was not sized from patch_nlipids')
        for resname in set(upper_fracs) | set(lower_fracs):
            # Bilayer truncates each species' share (int(frac * n)), so this is exact
            expected = (int(upper_fracs.get(resname, 0.0) * n_upper)
                        + int(lower_fracs.get(resname, 0.0) * n_lower))
            built = counts[resname]
            self.assertLessEqual(built, expected,
                                 f'{resname}: built {built} > sized {expected}; all {counts}')
            self.assertGreaterEqual(built, int(0.9 * expected),
                                    f'{resname}: built {built} << sized {expected}; all {counts}')

    @pytest.mark.slow
    def test_makemembranesystem_popc(self):
        test_dir = '__test_makemembranesystem_popc'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)

        task_list = [
            {'make_membrane_system': {
                'bilayer': {
                    'SAPL': 75,
                    'npatch': [1, 1],
                    'composition': {
                        'upper_leaflet': [{'name': 'POPC', 'frac': 1.0, 'conf': 0}],
                        'lower_leaflet': [{'name': 'POPC', 'frac': 1.0, 'conf': 0}]
                    },
                    'relaxation_protocols':{
                        'patch': self.common_patch_relaxation_protocols,
                        'quilt': self.common_quilt_relaxation_protocols
                    }
                }
            }},
            {'validate': {
                'tests': [
                    {'residue_test': {
                        'name': 'Count of lipid residues',
                        'selection': 'resname POPC',
                        'measure': 'residue_count',
                        'value': 200
                    }}
                ]
            }}]
        self.controller.reconfigure_tasks(task_list)
        self.assertIsInstance(self.controller.tasks[0], MakeMembraneSystemTask)
        self.assertIsInstance(self.controller.tasks[1], ValidateTask)
        self.assertEqual(len(self.controller.tasks), 2)
        result = self.controller.do_tasks()
        task = self.controller.tasks[0]
        validation = self.controller.tasks[1]
        state: StateArtifacts = task.get_current_artifact('state')
        self.assertTrue(state is not None)
        self.assertTrue(state.pdb.exists())
        self.assertTrue(state.psf.exists())
        self.assertTrue(state.coor.exists())
        self.assertTrue(state.xsc.exists())

        validation_results = validation.get_current_artifact('validation_results')
        self.assertIsInstance(validation_results, DataArtifact)
        self.assertIsInstance(validation_results.data, dict)
        self.assertEqual(validation_results.data['npass'], 1)
        self.assertEqual(validation_results.data['nfail'], 0)

        self.assertIsInstance(result, dict)
        os.chdir('..')
        assert result[0]['result'] == 0

    @pytest.mark.slow
    def test_makemembranesystem_popc_pope(self):
        test_dir='__test_makemembranesystem_popc_pope'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)
        task_list = [
            {'make_membrane_system': {
                'bilayer': {
                    'SAPL': 75,
                    'patch_nlipids': {
                        'upper': 49,
                        'lower': 49
                    },
                    'npatch': [1, 1],
                    'composition': {
                        'upper_leaflet': [{'name': 'POPC', 'frac': 1.0, 'conf': 0}],
                        'lower_leaflet': [{'name': 'POPE', 'frac': 1.0, 'conf': 0}]
                    },
                    'relaxation_protocols':{
                        'patch': self.common_patch_relaxation_protocols,
                        'quilt': self.common_quilt_relaxation_protocols
                    }
                }
            }},
            {'validate': {
                'tests': [
                    {
                        'residue_test': {
                            'name': 'Count of POPC residues',
                            'selection': 'resname POPC',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 49
                        }
                    },
                    {
                        # NOT 49.  POPE is the LOWER leaflet, whose count the task derives from
                        # the stress-free ratio, not from patch_nlipids -- see
                        # _assert_leaflet_sizing, which checks it against what the task recorded.
                        # This one only asserts the leaflet exists and is not wildly oversized.
                        'residue_test': {
                            'name': 'Count of POPE residues',
                            'selection': 'resname POPE',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 60
                        }
                    }
                ]
            }}]
        self.controller.reconfigure_tasks(task_list)
        self.assertIsInstance(self.controller.tasks[0], MakeMembraneSystemTask)
        self.assertIsInstance(self.controller.tasks[1], ValidateTask)
        self.assertEqual(len(self.controller.tasks), 2)
        result = self.controller.do_tasks()
        task = self.controller.tasks[0]
        validation = self.controller.tasks[1]
        state: StateArtifacts = task.get_current_artifact('state')
        self.assertTrue(state is not None)
        self.assertTrue(state.pdb.exists())
        self.assertTrue(state.psf.exists())
        self.assertTrue(state.coor.exists())
        self.assertTrue(state.xsc.exists())

        validation_results = validation.get_current_artifact('validation_results')
        self.assertIsInstance(validation_results, DataArtifact)
        self.assertIsInstance(validation_results.data, dict)
        self.assertEqual(validation_results.data['npass'], 2)
        self.assertEqual(validation_results.data['nfail'], 0)

        self._assert_leaflet_sizing(task, {'POPC': 1.0}, {'POPE': 1.0}, 49)

        # every relaxation stage ran through the subcontroller, and each must reach the run
        # record: two calibration patches and the quilt, each running its protocol (a pierced-lipid
        # guard may splice in an extra minimize ahead of the first dynamics stage)
        from pestifer.core.run_record import protocol_from_tasks
        protocol = [p for p in protocol_from_tasks(self.controller.tasks)
                    if p.get('within') == 'make_membrane_system']
        for bilayer in ('patchA', 'patchB', 'quilt'):
            stages = [p for p in protocol if p['task'] == f'make_membrane_system-md-{bilayer}']
            ensembles = [p['ensemble'] for p in stages]
            self.assertEqual(ensembles[-2:], ['nvt', 'npt'], bilayer)
            self.assertEqual(set(ensembles[:-2]), {'minimize'}, bilayer)
            self.assertEqual(sum(p['steps'] for p in stages[-2:]), 2000, bilayer)

        self.assertIsInstance(result, dict)
        os.chdir('..')
        assert result[0]['result'] == 0

    @pytest.mark.slow
    @pytest.mark.xfail(
        strict=False,
        reason="Builds a 50% cholesterol patch at SAPL 75 -- far above where a PSM/CHL1 leaflet "
               "sits -- and patchB's NPT stage dies with NAMD's 'Periodic cell has become too "
               "small for original patch grid!' (2 of 3 runs, 2026-10-14). The tool says so "
               "before the build: 'phase Ld requested but sterol fraction is 50% "
               "(cholesterol-rich) -- expect the leaflet to order ... consider phase: Lo', for "
               "BOTH leaflets, in every run. The remedy is a composition-appropriate SAPL (or "
               "phase, or a longer patch protocol), which is a decision about what this test "
               "should build, not a defect in the code under test. Not marked strict because it "
               "does sometimes survive.")
    def test_makemembranesystem_popc_chl1_psm_chl1(self):
        test_dir = '__test_makemembranesystem_popc_chl1_psm_chl1'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)

        task_list = [
            {'make_membrane_system': {
                'bilayer': {
                    'SAPL': 75,
                    'patch_nlipids': {
                        'upper': 64,
                        'lower': 64
                    },
                    'npatch': [1, 1],
                'composition':{
                    'upper_leaflet': [
                        {'name':'POPC','frac':0.5,'conf':0},
                        {'name':'CHL1','frac':0.5,'conf':0}],
                    'lower_leaflet': [
                        {'name':'PSM','frac':0.5,'conf':0},
                        {'name':'CHL1','frac':0.5,'conf':0}]},
                    'relaxation_protocols':{
                        'patch': self.common_patch_relaxation_protocols,
                        'quilt': self.common_quilt_relaxation_protocols
                    }
                }
            }},
            {'validate': {
                'tests': [
                    {
                        'residue_test': {
                            'name': 'Count of POPC residues',
                            'selection': 'resname POPC',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 32
                        }
                    },
                    {
                        # NOT 32, and CHL1 NOT 64: PSM and half the cholesterol are in the LOWER
                        # leaflet, whose count the task derives from the stress-free ratio rather
                        # than from patch_nlipids.  _assert_leaflet_sizing below checks both
                        # leaflets against what the task recorded; these are loose ceilings.
                        'residue_test': {
                            'name': 'Count of PSM residues',
                            'selection': 'resname PSM',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 40
                        }
                    },
                    {
                        'residue_test': {
                            'name': 'Count of CHL1 residues',
                            'selection': 'resname CHL1',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 80
                        }
                    },
                ]
            }}]
        self.controller.reconfigure_tasks(task_list)
        self.assertIsInstance(self.controller.tasks[0], MakeMembraneSystemTask)
        self.assertIsInstance(self.controller.tasks[1], ValidateTask)
        self.assertEqual(len(self.controller.tasks), 2)
        result = self.controller.do_tasks()
        task = self.controller.tasks[0]
        validation = self.controller.tasks[1]
        state: StateArtifacts = task.get_current_artifact('state')
        self.assertTrue(state is not None)
        self.assertTrue(state.pdb.exists())
        self.assertTrue(state.psf.exists())
        self.assertTrue(state.coor.exists())
        self.assertTrue(state.xsc.exists())

        validation_results = validation.get_current_artifact('validation_results')
        self.assertIsInstance(validation_results, DataArtifact)
        self.assertIsInstance(validation_results.data, dict)
        self.assertEqual(validation_results.data['npass'], 3)
        self.assertEqual(validation_results.data['nfail'], 0)

        # CHL1 is in both leaflets, so neither leaflet is a single resname; the helper sums the
        # species that belong to each.
        self._assert_leaflet_sizing(task, {'POPC': 0.5, 'CHL1': 0.5},
                                    {'PSM': 0.5, 'CHL1': 0.5}, 64)

        self.assertIsInstance(result, dict)
        os.chdir('..')
        assert result[0]['result'] == 0

    def test_makemembranesystem_embed_with_orient(self):
        test_dir = '__test_makemembranesystem_embed_with_orient'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)
        # psf = '5e8w-proteinonly.psf'
        # pdb = '5e8w-proteinonly.pdb'
        psf = 'van3.psf'
        pdb = 'van3.pdb'
        bilayer_psf = 'equilibrate.psf'
        bilayer_pdb = 'equilibrate.pdb'
        bilayer_xsc = 'equilibrate.xsc'
        input_data_dir = '../../fixtures/embed_inputs'
        for ftype in [psf, pdb, bilayer_psf, bilayer_pdb, bilayer_xsc]:
            shutil.copy(os.path.join(input_data_dir, ftype), '.')

        vm: VMDScripter = self.scripters['vmd']
        vm.newscript('orient')
        _stage_orient_oracle(vm)
        vm.usescript('bilayer_orient', local=True)
        vm.writescript('orient')
        result = vm.runscript(psf=psf,
                              pdb=pdb,
                              z_head_group=protect_str_arg("protein and resid 126"),
                              z_tail_group=protect_str_arg("protein and resid 158"),
                              o='oriented')
        opdb = 'oriented.pdb'
        pg: PsfgenScripter = self.scripters['psfgen']
        pg.newscript('embedded')
        pg.usescript('bilayer_embed')
        pg.writescript('embedded', guesscoord=False, regenerate=True, force_exit=True, writepsf=False, writepdb=False)
        result = pg.runscript(psf=psf,
                              pdb=opdb,
                              bilayer_psf=bilayer_psf,
                              bilayer_pdb=bilayer_pdb,
                              bilayer_xsc=bilayer_xsc,
                              z_ref_group=protect_str_arg("protein and resid 142"),
                              z_value=0.0,
                              o='embedded')
        os.chdir('..')
        assert result == 0

    def test_orient_align_matches_bilayer_orient(self):
        """Regression: the new transrot-ALIGN orientation path must reproduce the coordinates of
        the former Orient::orient-based bilayer_orient script (both compute the minimal roll-free
        rotation of the tail->head axis onto z about the protein COM)."""
        from pestifer.objs.rottrans import RotTrans
        test_dir = '__test_orient_align_matches'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)
        try:
            psf, pdb = 'van3.psf', 'van3.pdb'
            for f in (psf, pdb):
                shutil.copy(os.path.join('../../fixtures/embed_inputs', f), '.')
            head, tail = 'protein and resid 126', 'protein and resid 158'

            vm: VMDScripter = self.scripters['vmd']
            # old path: bilayer_orient (Orient::orient)
            vm.newscript('orient_old')
            _stage_orient_oracle(vm)
            vm.usescript('bilayer_orient', local=True)
            vm.writescript('orient_old')
            assert vm.runscript(psf=psf, pdb=pdb,
                                z_head_group=protect_str_arg(head),
                                z_tail_group=protect_str_arg(tail),
                                o='oriented_old') == 0

            # new path: transrot ALIGN emitted through the shared scripter (mirrors orient_protein)
            align = RotTrans(movetype='ALIGN', source=[tail, head], target=[0.0, 0.0, 1.0])
            vm.newscript('orient_new')
            vm.load_psf_pdb(psf, pdb, new_molid_varname='mOR')
            vm.write_rottrans(align, molid='mOR')
            vm.write_pdb('oriented_new', 'mOR')
            vm.writescript('orient_new')
            assert vm.runscript() == 0

            def coords(path):
                out = []
                with open(path) as fh:
                    for line in fh:
                        if line[:6] in ('ATOM  ', 'HETATM'):
                            out.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                return np.array(out)

            a, b = coords('oriented_old.pdb'), coords('oriented_new.pdb')
            self.assertEqual(a.shape, b.shape)
            rmsd = np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1)))
            self.assertLess(rmsd, 1.0e-2,
                            f'ALIGN orientation differs from bilayer_orient by RMSD {rmsd:.4f} A')
        finally:
            os.chdir('..')

    @pytest.mark.slow
    def test_makemembranesystem_5e8w_psm_chl1_pope_chl1(self):
        test_dir='__test_makemembranesystem_5e8w_psm_chl1_pope_chl1'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)
        psf='5e8w-proteinonly.psf'
        pdb='5e8w-proteinonly.pdb'
        task_list = [
            {'continuation': {'psf': '5e8w-proteinonly.psf', 'pdb': '5e8w-proteinonly.pdb'}},
            {'make_membrane_system': {
                'bilayer': {
                    'SAPL': 70.0,
                    'composition': {
                        'upper_leaflet': [
                            {'name': 'PSM', 'frac': 0.5},
                            {'name': 'CHL1', 'frac': 0.5}
                        ],
                        'lower_leaflet': [
                            {'name': 'POPE', 'frac': 0.5},
                            {'name': 'CHL1', 'frac': 0.5}
                        ]
                    },
                    'relaxation_protocols': {
                        'patch': self.common_patch_relaxation_protocols,
                        'quilt': self.common_quilt_relaxation_protocols
                    }
                },
                'embed': {
                    'z_head_group': 'protein and resid 667',
                    'z_tail_group': 'protein and resid 710',
                    'z_ref_group': {
                        'text': 'protein and resid 696',
                        'z_value': 0.0
                    }
                }
            }},
            {'md': {
                'ensemble': 'minimize',
                'nsteps': 1000
            }},
            {'validate' : {
                'tests': [
                    {
                        'residue_test': {
                            'name': 'Count of PSM residues',
                            'selection': 'resname PSM',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 200
                        }
                    },
                    {
                        'residue_test': {
                            'name': 'Count of POPE residues',
                            'selection': 'resname POPE',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 200
                        }
                    },
                    {
                        'residue_test': {
                            'name': 'Count of CHL1 residues',
                            'selection': 'resname CHL1',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 400
                        }
                    },
                    {
                        'residue_test' : {
                            'name': 'Test of protein presence',
                            'selection': 'protein',
                            'measure': 'residue_count',
                            'relation': '>=',
                            'value': 1
                        }
                    }
                ]
            }
            }
        ]
        input_data_dir = '../../fixtures/embed_inputs'
        for ftype in [psf, pdb]:
            shutil.copy(os.path.join(input_data_dir,ftype),'.')
        self.controller.reconfigure_tasks(task_list)
        result = self.controller.do_tasks()
        os.chdir('..')
        assert result[0]['result'] == 0

    @pytest.mark.slow
    def test_makemembranesystem_5e8w_grid_asymmetric_embed(self):
        # asymmetric (PSM/CHL1 over POPE/CHL1) bilayer embedded around the 5e8w fragment,
        # built with the packmol-free grid packer: calibrate two symmetric patches, size the
        # box to the protein, grid at stress-free per-leaflet counts, stage the relaxation,
        # then embed. Exercises the asymmetric + embedding grid path end to end.
        test_dir = '__test_makemembranesystem_5e8w_grid_asymmetric_embed'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)
        psf = '5e8w-proteinonly.psf'
        pdb = '5e8w-proteinonly.pdb'
        task_list = [
            {'continuation': {'psf': psf, 'pdb': pdb}},
            {'make_membrane_system': {
                'bilayer': {
                    'packer': 'grid',
                    'SAPL': 70.0,
                    'patch_nlipids': {'upper': 36, 'lower': 36},
                    'composition': {
                        'upper_leaflet': [
                            {'name': 'PSM', 'frac': 0.5},
                            {'name': 'CHL1', 'frac': 0.5}
                        ],
                        'lower_leaflet': [
                            {'name': 'POPE', 'frac': 0.5},
                            {'name': 'CHL1', 'frac': 0.5}
                        ]
                    },
                    'relaxation_protocols': {
                        'patch': self.common_patch_relaxation_protocols,
                        'quilt': self.common_quilt_relaxation_protocols
                    }
                },
                'embed': {
                    'z_head_group': 'protein and resid 667',
                    'z_tail_group': 'protein and resid 710',
                    'z_ref_group': {
                        'text': 'protein and resid 696',
                        'z_value': 0.0
                    }
                }
            }},
            {'validate': {
                'tests': [
                    {'residue_test': {
                        'name': 'Count of PSM residues', 'selection': 'resname PSM',
                        'measure': 'residue_count', 'relation': '<=', 'value': 200}},
                    {'residue_test': {
                        'name': 'Count of POPE residues', 'selection': 'resname POPE',
                        'measure': 'residue_count', 'relation': '<=', 'value': 200}},
                    {'residue_test': {
                        'name': 'Count of CHL1 residues', 'selection': 'resname CHL1',
                        'measure': 'residue_count', 'relation': '<=', 'value': 400}},
                    {'residue_test': {
                        'name': 'Test of protein presence', 'selection': 'protein',
                        'measure': 'residue_count', 'relation': '>=', 'value': 1}}
                ]
            }}
        ]
        input_data_dir = '../../fixtures/embed_inputs'
        for ftype in [psf, pdb]:
            shutil.copy(os.path.join(input_data_dir, ftype), '.')
        self.controller.reconfigure_tasks(task_list)
        result = self.controller.do_tasks()
        task = self.controller.tasks[1]
        validation = self.controller.tasks[2]
        state: StateArtifacts = task.get_current_artifact('state')
        self.assertTrue(state is not None)
        self.assertTrue(state.pdb.exists())
        self.assertTrue(state.psf.exists())
        validation_results = validation.get_current_artifact('validation_results')
        self.assertEqual(validation_results.data['nfail'], 0)
        os.chdir('..')
        assert result[0]['result'] == 0

