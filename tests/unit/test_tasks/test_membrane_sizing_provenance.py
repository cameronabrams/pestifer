# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""The numbers that size a membrane are recorded, and they say where they came from.

``make_membrane_system`` called ``record_outcome`` **zero times**. ``protocol_from_tasks``
contributes an entry only for a task that recorded one, so the step that decides the whole
membrane's size was *absent entirely* from the protocol section of every ``run-record.json``.
Those numbers are decided at run time and appear nowhere in the config, so they survived only in
the build log — one ``tail`` away from gone. CLAUDE.md already records what that costs: a
phase/composition warning sat correct and unread in every ex17 build log for two releases.

The second half is provenance. The APL summary line read

    Calibrated preferred APL: upper X, lower Y (from relaxed symmetric patches)

unconditionally — true only because nothing can currently supply an APL. A per-leaflet override
was built and reverted (``ea2fb6d1`` / ``705a56d3``) and would sit *upstream* of that line, so
the summary a reader greps for would have labeled a hand-supplied number "calibrated". That is
worse than silence: it is a wrong provenance claim in the one place someone looks. The source is
now carried per leaflet and printed, so if an override returns the label follows it rather than
having to be remembered.

Recorded at ``_grid_membrane`` — the single point BOTH sizing paths pass through
(``build_grid_membrane`` symmetric and ``build_grid_membrane_asymmetric``) — so a third path
cannot be added later and quietly leave no trace. That is the mistake this fixes, not repeated.
"""
import inspect
import unittest
from unittest import mock

from pestifer.core.run_record import protocol_from_tasks
from pestifer.tasks.make_membrane_system import (_AREA_CONVERGENCE_TOL, MakeMembraneSystemTask,
                                                 _as_float)


class TestAsFloat(unittest.TestCase):

    def test_none_stays_none_so_record_outcome_drops_it(self):
        """An unmeasured drift must be ABSENT from the record, not recorded as a null."""
        self.assertIsNone(_as_float(None))

    def test_a_number_becomes_a_plain_float(self):
        import numpy as np
        self.assertIsInstance(_as_float(np.float64(1.5)), float)
        self.assertEqual(_as_float(2), 2.0)


class TestTheChokePointRecords(unittest.TestCase):
    """Both sizing paths go through `_grid_membrane`, so the common facts are recorded there."""

    def setUp(self):
        self.src = inspect.getsource(MakeMembraneSystemTask._grid_membrane)

    def test_it_records(self):
        self.assertIn('self.record_outcome(', self.src)

    def test_it_records_the_counts_and_the_box(self):
        for field in ('n_upper', 'n_lower', 'box_SAPL', 'aspect', 'embedding'):
            self.assertIn(field, self.src, f'{field} is not in the run record')

    def test_BOTH_sizing_paths_reach_it(self):
        """NEGATIVE CONTROL for recording in one branch only, which is how this started."""
        for meth in (MakeMembraneSystemTask.build_grid_membrane,
                     MakeMembraneSystemTask.build_grid_membrane_asymmetric):
            self.assertIn('_grid_membrane(', inspect.getsource(meth), meth.__name__)


class TestTheAsymmetricPathRecordsProvenance(unittest.TestCase):

    @staticmethod
    def _code(text):
        """Executable lines only.

        The first version of the negative control below matched the COMMENT that quotes the old
        string, not the code -- a guard aimed one level off its target, which is the failure
        CLAUDE.md names twice. Comments are stripped so it can only pass or fail on what runs.
        """
        return '\n'.join(ln.split('#', 1)[0] for ln in text.splitlines()
                         if ln.split('#', 1)[0].strip())

    def setUp(self):
        self.src = self._code(inspect.getsource(
            MakeMembraneSystemTask.build_grid_membrane_asymmetric))

    def test_the_apls_and_their_sources_are_recorded(self):
        for field in ('apl_upper', 'apl_lower', 'apl_upper_source', 'apl_lower_source'):
            self.assertIn(field, self.src, f'{field} is not in the run record')

    def test_the_drift_and_the_tolerance_it_was_judged_against_are_recorded(self):
        """A drift without its tolerance cannot be read later without guessing the tolerance."""
        self.assertIn('apl_upper_drift', self.src)
        self.assertIn('apl_lower_drift', self.src)
        self.assertIn('area_convergence_tol', self.src)

    def test_the_log_line_names_the_source_of_each_apl(self):
        self.assertIn("apl_source = {'upper': 'calibrated', 'lower': 'calibrated'}", self.src)
        self.assertIn("apl_source['upper']", self.src)
        self.assertIn("apl_source['lower']", self.src)

    def test_the_old_unconditional_CALIBRATED_claim_is_gone(self):
        """NEGATIVE CONTROL. The string is the defect: it asserts provenance it cannot know."""
        self.assertNotIn('Calibrated preferred APL', self.src)
        self.assertNotIn('from relaxed symmetric patches', self.src)

    def test_the_source_is_read_from_a_variable_not_hardcoded_in_the_message(self):
        """So an override, if it returns, changes the label by changing the value -- rather than
        requiring whoever adds it to remember this line exists."""
        line = next(l for l in self.src.splitlines() if 'Preferred APL' in l)
        self.assertIn('_src', line)


class TestItReachesTheRunRecord(unittest.TestCase):
    """`record_outcome` only populates `task.outcome`; the record is built elsewhere. Pin the
    whole path, because a fact recorded into a dict nobody serializes is the same as no fact."""

    def test_a_recorded_outcome_becomes_a_protocol_entry(self):
        t = mock.Mock(taskname='make_membrane_system', substage_outcomes=[])
        t.outcome = {'apl_upper': 41.58, 'apl_upper_source': 'calibrated', 'n_upper': 300}
        entries = protocol_from_tasks([t])
        self.assertEqual(len(entries), 1)
        self.assertEqual(entries[0]['apl_upper'], 41.58)
        self.assertEqual(entries[0]['apl_upper_source'], 'calibrated')

    def test_a_task_that_records_NOTHING_contributes_nothing(self):
        """The old state, stated as a test: this is why the step was missing from every record."""
        t = mock.Mock(taskname='make_membrane_system', substage_outcomes=[])
        t.outcome = {}
        self.assertEqual(protocol_from_tasks([t]), [])


if __name__ == '__main__':
    unittest.main()
