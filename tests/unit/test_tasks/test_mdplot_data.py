# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Unit tests for mdplot's data preparation and figure content.

Plotting code is usually left untested because the obvious approach -- comparing rendered PNGs
against baselines -- is hostage to matplotlib, FreeType and font versions, and answers "something
changed" rather than "what changed".  Worse, a baseline generated from buggy output blesses the
bug.

So nothing here compares pixels.  Two things are tested instead:

* the **data preparation**, which is where unit conversion, series selection and smoothing
  decisions live.  These are pure functions of data and need no figure at all.
* the **content of the figure object** -- what was plotted, what it was labelled, what the axes
  say -- via matplotlib's own inspectable model.

That second kind is not hypothetical.  A build in this repository emitted

    mdplot.py:806: UserWarning: No artists with labels found to put in legend.

which meant a multi-trace plot had been produced with a legend containing nothing, because every
requested trace was missing from the data -- and a blank PNG was saved and registered as a build
artifact anyway, explained only at debug level.  A pixel comparison would not have caught that; an
assertion on the axes does, in one line.
"""

import logging
import os
import tempfile
import unittest
from unittest import mock

import matplotlib
matplotlib.use('Agg')                  # no display, no file writing during these tests
import matplotlib.pyplot as plt
import pandas as pd

from pestifer.tasks.mdplot import (
    MDPlotTask, auto_block_average, select_commensurable_runs, stage_of,
)


# --- pure data preparation ---------------------------------------------------------------------

class TestAutoBlockAverage(unittest.TestCase):
    """Automatic smoothing must touch the noisy observables and leave the rest alone."""

    def test_only_noisy_columns_are_smoothed(self):
        self.assertGreater(auto_block_average('PRESSURE', 10_000), 0)
        self.assertGreater(auto_block_average('GPRESSURE', 10_000), 0)
        self.assertEqual(auto_block_average('density', 10_000), 0,
                         'density is not noisy enough to warrant automatic smoothing')
        self.assertEqual(auto_block_average('TEMP', 10_000), 0)

    def test_column_matching_is_case_insensitive(self):
        self.assertEqual(auto_block_average('pressure', 10_000),
                         auto_block_average('PRESSURE', 10_000))

    def test_a_short_series_is_not_smoothed_away(self):
        # smoothing a warm-up shorter than twice the window would flatten the thing being looked at
        self.assertEqual(auto_block_average('PRESSURE', 8), 0)

    def test_the_window_is_bounded_at_both_ends(self):
        small = auto_block_average('PRESSURE', 1_000)
        huge = auto_block_average('PRESSURE', 10_000_000)
        self.assertGreaterEqual(small, 5)
        self.assertLessEqual(huge, 200, 'the window must not grow without bound')


class TestStageOf(unittest.TestCase):
    """Stage boundaries are recovered from the filename, so the standalone re-plot path -- which
    has no pipeline and hence no MD history -- can still draw them."""

    def test_a_run_basename_parses_into_stage_and_label(self):
        self.assertEqual(stage_of('00-07-003_density_equilibrate-NPT'),
                         ('00-07', 'density_equilibrate-NPT'))

    def test_runs_of_one_task_share_a_stage_key(self):
        a, _ = stage_of('00-07-000_density_equilibrate-NPT')
        b, _ = stage_of('00-07-019_density_equilibrate-NPT')
        self.assertEqual(a, b, 'chunks of one equilibration are one stage')

    def test_different_tasks_are_different_stages(self):
        self.assertNotEqual(stage_of('00-06-000_md-nvt')[0],
                            stage_of('00-07-000_density_equilibrate-NPT')[0])

    def test_an_unparseable_name_is_not_a_stage(self):
        self.assertEqual(stage_of('some-other-file'), (None, None))


class TestSelectCommensurableRuns(unittest.TestCase):
    """Only dynamics on *this* system belong in the plot: a solvate between two MD stages changes
    the atom count and ends the span."""

    def test_a_change_in_atom_count_ends_the_span(self):
        entries = [{'natoms': 1000}, {'natoms': 1000}, {'natoms': 14000}, {'natoms': 14000}]
        self.assertEqual(select_commensurable_runs(entries, 14000),
                         [{'natoms': 14000}, {'natoms': 14000}])

    def test_the_trailing_span_is_kept_in_order(self):
        entries = [{'natoms': 14000, 'i': 1}, {'natoms': 14000, 'i': 2}]
        self.assertEqual([e['i'] for e in select_commensurable_runs(entries, 14000)], [1, 2])

    def test_an_unknown_count_is_not_treated_as_a_boundary(self):
        # an unrecorded count is not evidence of a different system; treating it as one silently
        # truncated the plot, which is the failure this replaces
        entries = [{'natoms': 14000}, {'natoms': None}, {'natoms': 14000}]
        self.assertEqual(len(select_commensurable_runs(entries, 14000)), 3)

    def test_everything_is_kept_when_the_system_size_is_unknown(self):
        entries = [{'natoms': 1000}, {'natoms': 14000}]
        self.assertEqual(select_commensurable_runs(entries, None), entries)


class TestUnitResolution(unittest.TestCase):
    """A histogram that ignored the unit spec would label a density axis 0.62 where the timeseries
    beside it reads 1.03 -- the same quantity in two units in one figure set."""

    @staticmethod
    def _task(units=None):
        t = MDPlotTask.__new__(MDPlotTask)
        t.specs = {'units': units or {}}
        return t

    def test_density_is_converted_to_g_per_cc(self):
        df = pd.DataFrame({'density': [0.62]})
        key, got, mult, label = self._task({'density': 'g/cc'})._resolve_units(
            'density', {'density': df})
        self.assertEqual(key, 'density')
        self.assertAlmostEqual(0.62 * mult, 1.03, places=1)
        self.assertEqual(label, 'g/cc')

    def test_a_namd_quantity_gets_its_documented_unit(self):
        df = pd.DataFrame({'TEMP': [300.0]})
        _key, _df, mult, label = self._task()._resolve_units('TEMP', {'TEMP': df})
        self.assertEqual((mult, label), (1.0, 'K'))

    def test_large_energies_are_rescaled_and_relabelled_together(self):
        # the multiplier and the label must move as a pair, or the axis lies
        df = pd.DataFrame({'TOTAL': [-250_000.0, -251_000.0]})
        _key, _df, mult, label = self._task()._resolve_units('TOTAL', {'TOTAL': df})
        self.assertAlmostEqual(mult, 1e-3)
        self.assertEqual(label, '1000 kcal/mol')

    def test_a_missing_quantity_resolves_to_no_dataframe(self):
        _key, df, _mult, _label = self._task()._resolve_units('nonesuch', {})
        self.assertIsNone(df, 'a missing quantity must be reported, not invented')

    def test_lookup_is_case_insensitive_in_both_directions(self):
        upper = pd.DataFrame({'TEMP': [300.0]})
        lower = pd.DataFrame({'density': [1.0]})
        self.assertIsNotNone(self._task()._resolve_units('temp', {'TEMP': upper})[1])
        self.assertIsNotNone(self._task()._resolve_units('DENSITY', {'density': lower})[1])


# --- figure content ------------------------------------------------------------------------------

class TestEmptyPlotHandling(unittest.TestCase):
    """Regression: when every requested trace was missing from the data, the loop drew nothing,
    called ``ax.legend()`` on empty axes (matplotlib's "No artists with labels" warning), and saved
    a blank PNG which it then registered as a build artifact -- explained only at debug level."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmpdir = self._tmp.name
        self.addCleanup(self._tmp.cleanup)
        self.addCleanup(plt.close, 'all')

    def _task(self, timeseries):
        t = MDPlotTask.__new__(MDPlotTask)
        t.taskname = 'mdplot'
        t.basename = 'test'
        t.specs = {'timeseries': timeseries, 'legend': True, 'grid': False,
                   'stage-markers': False, 'block-average': 0, 'units': {}}
        t.colormap = plt.get_cmap('tab10')
        t.colormap_direction = 1
        t.register = mock.Mock()
        return t

    def _plot(self, task, df_of_column):
        """Call the extracted timeseries plotter and hand back the figure it built."""
        captured = {}

        def capture_clf(*a, **kw):
            # capture the figure and deliberately do NOT clear it: the assertions inspect the
            # axes this call was about to discard.  tearDown closes everything.
            captured.setdefault('fig', plt.gcf())

        with mock.patch.object(plt, 'clf', capture_clf):
            task._plot_timeseries(task.specs['timeseries'], df_of_column, self.tmpdir,
                                  block_average=0, stage_markers=False, legend=True, grid=False)
        return captured.get('fig')

    def test_no_blank_png_is_written_when_every_trace_is_missing(self):
        task = self._task([['density', 'a_x']])
        with self.assertLogs('pestifer.tasks.mdplot', level='WARNING') as cm:
            self._plot(task, {})
        self.assertFalse(os.listdir(self.tmpdir), 'a blank plot was written anyway')
        self.assertIn('no data', ''.join(cm.output).lower())
        task.register.assert_not_called()

    def test_a_partial_trace_set_still_plots_and_warns(self):
        df = pd.DataFrame({'TS': [0, 100, 200], 'density': [1.0, 1.01, 1.0]})
        task = self._task([['density', 'nonesuch']])
        with self.assertLogs('pestifer.tasks.mdplot', level='WARNING') as cm:
            fig = self._plot(task, {'density': df})
        self.assertIn('nonesuch', ''.join(cm.output))
        self.assertEqual(len(fig.axes[0].get_lines()), 1, 'the available trace should still plot')

    def test_the_legend_is_never_empty(self):
        """The direct assertion matplotlib's warning was asking for."""
        df = pd.DataFrame({'TS': [0, 100, 200], 'density': [1.0, 1.01, 1.0]})
        task = self._task([['density', 'nonesuch']])
        with self.assertLogs('pestifer.tasks.mdplot', level='WARNING'):
            fig = self._plot(task, {'density': df})
        legend = fig.axes[0].get_legend()
        if legend is not None:
            self.assertTrue(legend.get_texts(), 'a legend was drawn with nothing in it')

    def test_a_fully_available_trace_set_plots_and_labels_both(self):
        df_d = pd.DataFrame({'TS': [0, 100, 200], 'density': [1.0, 1.01, 1.0]})
        df_a = pd.DataFrame({'TS': [0, 100, 200], 'a_x': [55.0, 55.1, 55.0]})
        task = self._task([['density', 'a_x']])
        fig = self._plot(task, {'density': df_d, 'a_x': df_a})
        ax = fig.axes[0]
        self.assertEqual(len(ax.get_lines()), 2)
        _handles, labels = ax.get_legend_handles_labels()
        self.assertEqual(len(labels), 2, 'both traces should be labelled for the legend')
        self.assertIn('density', ax.get_ylabel())

    def test_a_png_is_written_and_registered_when_there_is_data(self):
        df = pd.DataFrame({'TS': [0, 100, 200], 'density': [1.0, 1.01, 1.0]})
        task = self._task(['density'])
        self._plot(task, {'density': df})
        written = os.listdir(self.tmpdir)
        self.assertEqual(len(written), 1, f'expected one figure, got {written}')
        task.register.assert_called_once()
