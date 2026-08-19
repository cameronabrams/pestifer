# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Unit tests for the chunked-equilibration loop (``ChunkedEquilibrateTask.do``).

The convergence *criterion* is well covered elsewhere (``test_density_convergence.py``); what was
not covered at all is the machinery around it -- the loop that decides how long each chunk should
be, recovers when NAMD outruns its patch grid, decides why it stopped, and reports what it
actually ran.  That code only ran under a live NAMD, so it was exercised end to end and asserted
on nowhere.

It is testable without NAMD.  ``ChunkedEquilibrateTask`` already isolates everything
physics-dependent behind hooks a subclass implements (``_setup``, ``_ingest_chunk``,
``_converged_stop_reason``, ...), so a stub subclass can script the verdicts and drive the loop
directly.  These tests are therefore about control flow and bookkeeping: chunk sizing, retry
behaviour, stop reasons, and the numbers that end up in ``run-record.json``.
"""

import os
import unittest
from unittest import mock

from pestifer.core.errors import PestiferBuildError
from pestifer.tasks import equilibrate_base as EB
from pestifer.tasks.equilibrate_base import ChunkedEquilibrateTask, _MIN_RETRY_CHUNK


_DEFAULT_SPECS = dict(margin=4.0, shrink_safety=0.5, chunk_min=500, chunk_max=5000,
                      max_steps=6000, min_steps=2000, chunk_growth=1.5, max_shrink_retries=5)


class _StubTask(ChunkedEquilibrateTask):
    """A ChunkedEquilibrateTask whose physics is scripted rather than simulated.

    ``verdicts`` is consumed one entry per chunk as ``(blowup, converged)``; when it runs out the
    task simply keeps going, which is how a run reaches its ceiling.  ``crash_on`` names chunk
    numbers (1-based) on which ``namdrun`` should raise as NAMD does when the shrinking cell
    outruns its patch grid.
    """

    _ensemble = 'NPT'

    def __init__(self, tmpdir, specs=None, verdicts=(), crash_on=(), crash_text=None,
                 firsttimestep=0):
        self.taskname = 'stub_equilibrate'
        self.basename = os.path.join(tmpdir, 'stub')
        self.specs = dict(_DEFAULT_SPECS, **(specs or {}))
        self.namd_global_config = {'solvated': {'stepspercycle': 10}}
        self.outcome = {}
        self._rows = []
        self._verdicts = list(verdicts)
        self._crash_on = set(crash_on)
        # NAMD's actual abort; is_patch_grid_crash matches the 'too small' + 'patch grid'
        # fragments, so the stub must carry both or the retry path is never reached
        self._crash_text = (crash_text if crash_text is not None else
                            'FATAL ERROR: Periodic cell has become too small for original '
                            'patch grid!  Possible solutions are to restart from a recent '
                            'checkpoint, increase margin, or disable useFlexibleCell.')
        self._firsttimestep = firsttimestep
        # observations for the assertions
        self.ran_chunks = []          # nsteps passed to each namdrun call
        self.written = None           # stop reason handed to _write_outputs
        self.namdrun_calls = 0

    # --- the bits the loop pulls from task/pipeline state -------------------------------------
    def get_current_artifact(self, key):
        return mock.Mock(psf=mock.Mock(name='stub.psf'))

    def get_current_artifact_data(self, key):
        return self._firsttimestep if key == 'firsttimestep' else None

    def record_outcome(self, **facts):
        self.outcome.update({k: v for k, v in facts.items() if v is not None})

    # --- NAMD, scripted -------------------------------------------------------------------------
    def namdrun(self, *a, **kw):
        # `crash_on` names namdrun *call* numbers, not chunk numbers: a retried chunk keeps its
        # chunk number, so keying off that would crash the same chunk forever.
        self.namdrun_calls += 1
        if self.namdrun_calls in self._crash_on:
            with open(f'{self.basename}.log', 'w') as f:
                f.write(self._crash_text)
            raise PestiferBuildError('namd died')
        self.ran_chunks.append(self.specs['nsteps'])   # only a completed chunk contributes steps
        # a successful chunk leaves an .xst behind
        with open(f'{self.basename}.xst', 'w') as f:
            f.write('# stub\n')
        return 0

    # --- hooks --------------------------------------------------------------------------------
    def _setup(self, state, min_steps):
        self.setup_min_steps = min_steps

    def _ingest_chunk(self, xst, total_steps, this_chunk, n_chunk):
        return self._verdicts.pop(0) if self._verdicts else (False, False)

    def _converged_stop_reason(self, total_steps):
        return f'CONVERGED: stub at {total_steps}'

    def _ceiling_stop_reason(self, total_steps, max_steps):
        return f'CEILING: reached max_steps ({max_steps}) at {total_steps}'

    def _blowup_stop_reason(self, total_steps):
        return f'BLOWUP at {total_steps}'

    def _blowup_error(self, total_steps):
        return f'blew up at {total_steps}'

    def _write_outputs(self, stop_reason):
        self.written = stop_reason

    def _log_patch_grid(self, log_path):
        pass


class _LoopTestCase(unittest.TestCase):
    """Shared scaffolding: a temp dir, and a fixed shrink rate so chunk sizing is deterministic."""

    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.tmpdir = self._tmp.name
        # the loop sizes each next chunk from the observed shrink rate; pin it so the only
        # variable under test is the loop's own arithmetic
        patcher = mock.patch.object(EB, 'xst_max_shrink_rate', return_value=0.0)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.addCleanup(self._tmp.cleanup)


class TestStoppingConditions(_LoopTestCase):

    def test_converged_verdict_stops_the_loop(self):
        t = _StubTask(self.tmpdir, verdicts=[(False, False), (False, True)])
        self.assertEqual(t.do(), 0)
        self.assertEqual(len(t.ran_chunks), 2, 'loop did not stop on the converged verdict')
        self.assertTrue(t.converged)
        self.assertTrue(t.written.startswith('CONVERGED'))

    def test_ceiling_reached_when_never_converging(self):
        t = _StubTask(self.tmpdir)          # no verdicts -> never converges
        self.assertEqual(t.do(), 0)
        self.assertFalse(t.converged)
        self.assertTrue(t.written.startswith('CEILING'))

    def test_blowup_writes_outputs_then_raises(self):
        # the outputs matter: a blown-up run must still leave its diagnostics behind
        t = _StubTask(self.tmpdir, verdicts=[(True, False)])
        with self.assertRaises(PestiferBuildError):
            t.do()
        self.assertEqual(t.written, 'BLOWUP at 500')

    def test_missing_xst_is_an_error(self):
        t = _StubTask(self.tmpdir)
        orig = t.namdrun

        def namdrun_without_xst(*a, **kw):
            orig(*a, **kw)
            os.remove(f'{t.basename}.xst')
        t.namdrun = namdrun_without_xst
        with self.assertRaises(PestiferBuildError) as cm:
            t.do()
        self.assertIn('no .xst', str(cm.exception))


class TestBudgetAccounting(_LoopTestCase):

    def test_total_steps_respects_the_budget(self):
        t = _StubTask(self.tmpdir, specs=dict(max_steps=6000))
        t.do()
        self.assertLessEqual(sum(t.ran_chunks), 6000 + t.specs['chunk_max'])
        self.assertEqual(t.outcome['steps'], sum(t.ran_chunks))

    def test_max_steps_is_a_budget_not_an_absolute_ceiling(self):
        """A run inheriting a large firsttimestep must still get its full budget.

        Regression: an asymmetric build's second calibration patch runs on an accumulating step
        counter, and treating max_steps as absolute starved it after a few thousand steps and
        reported an unrelaxed area.
        """
        t = _StubTask(self.tmpdir, specs=dict(max_steps=6000), firsttimestep=50_000)
        t.do()
        self.assertGreater(sum(t.ran_chunks), 5000,
                           'inherited firsttimestep starved the run of its budget')

    def test_outcome_records_what_actually_ran(self):
        t = _StubTask(self.tmpdir, verdicts=[(False, False), (False, False), (False, True)])
        t.do()
        self.assertEqual(t.outcome['chunks'], 3)
        self.assertEqual(t.outcome['steps'], sum(t.ran_chunks))
        self.assertTrue(t.outcome['adaptive'])
        self.assertTrue(t.outcome['converged'])
        self.assertIn('CONVERGED', t.outcome['stopped_because'])

    def test_every_chunk_is_a_whole_number_of_cycles(self):
        t = _StubTask(self.tmpdir, specs=dict(max_steps=9000))
        t.do()
        spc = t.namd_global_config['solvated']['stepspercycle']
        self.assertTrue(all(c % spc == 0 for c in t.ran_chunks),
                        f'NAMD requires whole cycles; got {t.ran_chunks}')

    def test_first_chunk_is_the_conservative_minimum(self):
        t = _StubTask(self.tmpdir, specs=dict(chunk_min=500))
        t.do()
        self.assertEqual(t.ran_chunks[0], 500)

    def test_chunk_growth_is_capped(self):
        t = _StubTask(self.tmpdir, specs=dict(chunk_growth=1.5, max_steps=20_000))
        t.do()
        for prev, nxt in zip(t.ran_chunks, t.ran_chunks[1:]):
            self.assertLessEqual(nxt, int(prev * 1.5) + 1,
                                 f'chunk grew faster than chunk_growth: {t.ran_chunks}')


class TestPatchGridRecovery(_LoopTestCase):

    def test_crash_halves_the_chunk_and_retries(self):
        t = _StubTask(self.tmpdir, crash_on=[1])
        t.do()
        self.assertGreaterEqual(t.namdrun_calls, 2, 'no retry after a patch-grid crash')
        self.assertEqual(t.ran_chunks[0], 250, 'retry did not halve the chunk')

    def test_retries_are_bounded_and_the_error_is_actionable(self):
        # crash on every attempt: the loop must give up rather than spin
        t = _StubTask(self.tmpdir, specs=dict(max_shrink_retries=3), crash_on=range(1, 50))
        with self.assertRaises(PestiferBuildError) as cm:
            t.do()
        msg = str(cm.exception)
        self.assertIn('margin', msg, 'the error should name the knob that fixes it')
        self.assertLessEqual(t.namdrun_calls, 10, 'retries are not bounded')

    def test_an_unrelated_namd_failure_is_re_raised_unchanged(self):
        t = _StubTask(self.tmpdir, crash_on=[1], crash_text='FATAL ERROR: something else entirely')
        with self.assertRaises(PestiferBuildError) as cm:
            t.do()
        self.assertEqual(str(cm.exception), 'namd died',
                         'a non-patch-grid failure must not be reinterpreted as a grid crash')

    def test_a_successful_chunk_clears_the_retry_counter(self):
        # crash once, succeed, then crash again: the second crash must still get its own retries
        t = _StubTask(self.tmpdir, specs=dict(max_shrink_retries=1), crash_on=[1, 3])
        t.do()
        self.assertGreaterEqual(len(t.ran_chunks), 2,
                                'a later crash was not retried, so the counter never reset')

    def test_retry_never_goes_below_the_hard_floor(self):
        t = _StubTask(self.tmpdir, specs=dict(chunk_min=40, max_shrink_retries=5),
                      crash_on=[1, 2, 3])
        try:
            t.do()
        except PestiferBuildError:
            pass
        self.assertTrue(all(c >= _MIN_RETRY_CHUNK for c in t.ran_chunks),
                        f'chunk shrank below the {_MIN_RETRY_CHUNK}-step floor: {t.ran_chunks}')
