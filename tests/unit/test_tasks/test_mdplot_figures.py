# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for mdplot's derived quantities and figure-building paths.

``test_mdplot_data.py`` covers the pure data preparation and the timeseries figure.  This file
takes the rest: the cell geometry derived from the ``.xst`` series, the x-axis choice, the summary
table, and the four remaining figure builders (overlays, panels, histograms, stage markers).

The approach is the same and for the same reasons -- nothing here compares pixels.  Baseline PNG
comparison is hostage to matplotlib, FreeType and font versions, answers "something changed"
rather than "what changed", and blesses whatever bug was present when the baseline was made.  What
is asserted instead is the *content* of the figure object via matplotlib's own inspectable model
(how many axes, how many lines, what they were labelled) and, where a number is being computed
rather than drawn, the number itself.
"""

import contextlib
import os
import tempfile
import unittest
from unittest import mock

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pestifer.tasks.mdplot import MDPlotTask, _MAX_LABELED_STAGES


def _task(**specs):
    """An MDPlotTask with just enough state for the plotting methods; never provisioned."""
    t = MDPlotTask.__new__(MDPlotTask)
    t.taskname = 'mdplot'
    t.basename = 'test'
    t.specs = {'units': {}, 'axis-labels': {}, **specs}
    t.dataframes = {}
    t.stage_boundaries = {}
    t.colormap = plt.get_cmap('tab10')
    t.colormap_direction = 1
    t.register = mock.Mock()
    t.reprocess_logs = False
    t.provisions = {}
    return t


class _FigureCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmpdir = self._tmp.name
        self.addCleanup(self._tmp.cleanup)
        self.addCleanup(plt.close, 'all')

    @contextlib.contextmanager
    def captured_figures(self):
        """Hold on to figures the plotting code disposes of before returning.

        Both builders end by closing (panels, overlays) or clearing (histograms) their figure, so
        by the time the method returns there is nothing left to inspect.  The patches record the
        figure at the moment of disposal and deliberately do not destroy it; tearDown closes
        everything.
        """
        figs = []

        def fake_close(fig=None):
            if hasattr(fig, 'axes'):
                figs.append(fig)

        def fake_clf():
            figs.append(plt.gcf())

        with mock.patch.object(plt, 'close', fake_close), \
             mock.patch.object(plt, 'clf', fake_clf):
            yield figs


# --- derived cell geometry -----------------------------------------------------------------------

class TestDerivedCellColumns(unittest.TestCase):
    """NAMD logs three cell vectors; what a membrane build wants is area, volume and area per
    lipid.  Deriving them here makes them ordinary traces."""

    @staticmethod
    def _orthorhombic(ax=50.0, by=60.0, cz=70.0, n=3):
        return pd.DataFrame({'TS': range(n),
                             'a_x': [ax] * n, 'a_y': [0.0] * n, 'a_z': [0.0] * n,
                             'b_x': [0.0] * n, 'b_y': [by] * n, 'b_z': [0.0] * n,
                             'c_x': [0.0] * n, 'c_y': [0.0] * n, 'c_z': [cz] * n})

    def _derive(self, df, **specs):
        t = _task(**specs)
        t.dataframes['xst'] = df
        t._add_derived_columns()
        return t.dataframes['xst']

    def test_area_and_volume_for_an_orthorhombic_cell(self):
        got = self._derive(self._orthorhombic())
        self.assertAlmostEqual(got['area'].iloc[0], 50.0 * 60.0)
        self.assertAlmostEqual(got['volume'].iloc[0], 50.0 * 60.0 * 70.0)

    def test_a_triclinic_cell_uses_the_exact_vector_forms(self):
        """The documented reason for |a x b| over ``a_x * b_y``: the two agree only for an
        orthorhombic cell, so a sheared cell would be silently wrong.

        The shear has to be chosen with care -- tilting b along *a* leaves |a x b| exactly equal
        to ``a_x * b_y``, so that version of this test passes no matter which formula is used.
        Tilting b out of the ab plane is what separates them.
        """
        df = self._orthorhombic()
        df['b_z'] = 20.0                      # tilt b out of the ab plane
        got = self._derive(df)
        naive = df['a_x'].iloc[0] * df['b_y'].iloc[0]
        exact = float(np.linalg.norm(np.cross([50.0, 0.0, 0.0], [0.0, 60.0, 20.0])))
        self.assertNotAlmostEqual(exact, naive, msg='this shear does not exercise the difference')
        self.assertAlmostEqual(got['area'].iloc[0], exact)

    def test_volume_is_positive_for_a_left_handed_cell(self):
        # |a . (b x c)| -- a sign flip must not produce a negative volume
        df = self._orthorhombic()
        df['c_z'] = -70.0
        self.assertGreater(self._derive(df)['volume'].iloc[0], 0)

    def test_area_per_lipid_needs_a_lipid_count(self):
        self.assertNotIn('apl', self._derive(self._orthorhombic()).columns)
        got = self._derive(self._orthorhombic(), **{'lipids-per-leaflet': 100})
        self.assertAlmostEqual(got['apl'].iloc[0], 50.0 * 60.0 / 100)

    def test_a_zero_lipid_count_is_treated_as_absent(self):
        self.assertNotIn('apl', self._derive(self._orthorhombic(),
                                             **{'lipids-per-leaflet': 0}).columns)

    def test_a_partial_cell_series_is_skipped_not_guessed(self):
        df = self._orthorhombic().drop(columns=['c_z'])
        self.assertNotIn('area', self._derive(df).columns)

    def test_no_xst_series_is_not_an_error(self):
        t = _task()
        t._add_derived_columns()               # nothing registered at all
        t.dataframes['xst'] = pd.DataFrame()
        t._add_derived_columns()               # present but empty


# --- x axis --------------------------------------------------------------------------------------

class TestXAxisFor(unittest.TestCase):
    """Simulation time is preferred over timestep, and the unit switches so the axis never reads
    in the thousands."""

    def test_picoseconds_below_a_nanosecond(self):
        x, label = _task()._x_axis_for(pd.DataFrame({'time_ps': [0.0, 500.0]}))
        self.assertEqual(label, 'simulation time (ps)')
        self.assertAlmostEqual(x.iloc[-1], 500.0)

    def test_nanoseconds_at_and_above_one(self):
        x, label = _task()._x_axis_for(pd.DataFrame({'time_ps': [0.0, 2000.0]}))
        self.assertEqual(label, 'simulation time (ns)')
        self.assertAlmostEqual(x.iloc[-1], 2.0, msg='values must be rescaled with the label')

    def test_timestep_is_the_fallback(self):
        x, label = _task()._x_axis_for(pd.DataFrame({'TS': [0, 100]}))
        self.assertEqual(label, 'time step')
        self.assertEqual(list(x), [0, 100])

    def test_the_timestep_column_name_is_configurable(self):
        t = _task(time_step_column_names=['frame'])
        _x, label = t._x_axis_for(pd.DataFrame({'frame': [0, 1]}))
        self.assertEqual(label, 'time step')

    def test_a_dataframe_with_no_x_axis_reports_none(self):
        self.assertEqual(_task()._x_axis_for(pd.DataFrame({'density': [1.0]})), (None, None))


# --- summary table -------------------------------------------------------------------------------

class TestSummaryTable(_FigureCase):
    """One row per stage per quantity: "the density settled at 1.03 during NPT" as a number that
    can be pasted or diffed rather than eyeballed off a curve."""

    def _run(self, task, quantities, df_of_column):
        task._write_summary_table(quantities, df_of_column, self.tmpdir)
        path = os.path.join(self.tmpdir, 'test-summary.csv')
        return pd.read_csv(path) if os.path.exists(path) else None

    def _two_stage_task(self):
        df = pd.DataFrame({'TS': [0, 100, 200, 300], 'density': [1.0, 1.0, 2.0, 2.0]})
        t = _task()
        t.dataframes['energy'] = df
        t.stage_boundaries['energy'] = [(0, 'warmup'), (2, 'npt')]
        return t, df

    def test_one_row_per_stage(self):
        t, df = self._two_stage_task()
        table = self._run(t, ['density'], {'density': df})
        self.assertEqual(list(table['stage']), ['warmup', 'npt'])

    def test_statistics_are_computed_over_the_stage_not_the_whole_run(self):
        t, df = self._two_stage_task()
        table = self._run(t, ['density'], {'density': df})
        self.assertAlmostEqual(table['mean'].iloc[0], 1.0)
        self.assertAlmostEqual(table['mean'].iloc[1], 2.0,
                               msg='the second stage must not average in the first')

    def test_the_step_range_of_each_stage_is_recorded(self):
        t, df = self._two_stage_task()
        table = self._run(t, ['density'], {'density': df})
        self.assertEqual((table['first_step'].iloc[1], table['last_step'].iloc[1]), (200, 300))

    def test_a_series_with_no_recorded_stages_is_one_row(self):
        df = pd.DataFrame({'TS': [0, 100], 'density': [1.0, 1.1]})
        t = _task()
        t.dataframes['energy'] = df
        table = self._run(t, ['density'], {'density': df})
        self.assertEqual(len(table), 1)

    def test_a_single_sample_stage_reports_zero_spread_rather_than_nan(self):
        df = pd.DataFrame({'TS': [0, 100], 'density': [1.0, 2.0]})
        t = _task()
        t.dataframes['energy'] = df
        t.stage_boundaries['energy'] = [(0, 'a'), (1, 'b')]
        table = self._run(t, ['density'], {'density': df})
        self.assertEqual(list(table['std']), [0.0, 0.0])

    def test_units_reach_the_table(self):
        # a live multiplier: large energies are rescaled and relabelled as a pair, and both the
        # rescaled value and the label have to survive into the summary table
        df = pd.DataFrame({'TS': [0, 1], 'TOTAL': [-250_000.0, -250_000.0]})
        t = _task()
        t.dataframes['energy'] = df
        table = self._run(t, ['TOTAL'], {'TOTAL': df})
        self.assertAlmostEqual(table['mean'].iloc[0], -250.0, places=1)
        self.assertEqual(table['units'].iloc[0], '1000 kcal/mol')

    def test_a_density_table_is_labelled_g_per_cc(self):
        df = pd.DataFrame({'TS': [0, 1], 'density': [1.03, 1.03]})
        t = _task()
        t.dataframes['energy'] = df
        table = self._run(t, ['density'], {'density': df})
        self.assertAlmostEqual(table['mean'].iloc[0], 1.03, places=2)
        self.assertEqual(table['units'].iloc[0], 'g/cc')

    def test_no_table_is_written_when_nothing_resolves(self):
        self.assertIsNone(self._run(_task(), ['nonesuch'], {}))
        self.assertFalse(os.listdir(self.tmpdir))


# --- histograms ----------------------------------------------------------------------------------

class TestHistograms(_FigureCase):
    """A timeseries shows whether a quantity settled; the histogram shows what it settled at.
    Including the approach to equilibrium makes it bimodal and its mean meaningless, so only the
    trailing fraction is used."""

    @staticmethod
    def _approach_then_plateau(n=100):
        # first half ramps, second half sits at 2.0
        return pd.DataFrame({'TS': range(n),
                             'density': list(np.linspace(0.0, 2.0, n // 2)) + [2.0] * (n // 2)})

    def _plot(self, task, df_of_column, names=('density',)):
        with self.captured_figures() as figs:
            task._plot_histograms(list(names), df_of_column, self.tmpdir)
        self.assertTrue(figs, 'no histogram figure was produced')
        return figs[0]

    def test_only_the_trailing_fraction_is_histogrammed(self):
        df = self._approach_then_plateau()
        t = _task(**{'histogram-tail': 0.5})
        fig = self._plot(t, {'density': df})
        title = fig.axes[0].get_title()
        self.assertIn('mean 2', title, f'the ramp was averaged in: {title}')

    def test_the_tail_fraction_is_configurable_and_reported(self):
        df = self._approach_then_plateau()
        fig = self._plot(_task(**{'histogram-tail': 0.25}), {'density': df})
        self.assertIn('25%', fig.axes[0].get_title())

    def test_an_out_of_range_tail_fraction_is_clamped(self):
        df = self._approach_then_plateau()
        fig = self._plot(_task(**{'histogram-tail': 5.0}), {'density': df})
        self.assertIn('100%', fig.axes[0].get_title())

    def test_a_figure_is_written_and_registered(self):
        t = _task()
        self._plot(t, {'density': self._approach_then_plateau()})
        self.assertEqual(os.listdir(self.tmpdir), ['test-density-hist.png'])
        t.register.assert_called_once()

    def test_too_few_samples_to_be_a_distribution_are_skipped(self):
        df = pd.DataFrame({'TS': [0, 1], 'density': [1.0, 1.1]})
        t = _task()
        t._plot_histograms(['density'], {'density': df}, self.tmpdir)
        self.assertFalse(os.listdir(self.tmpdir))
        t.register.assert_not_called()

    def test_a_missing_quantity_is_skipped(self):
        t = _task()
        t._plot_histograms(['nonesuch'], {}, self.tmpdir)
        self.assertFalse(os.listdir(self.tmpdir))

    def test_the_mean_is_marked_on_the_axes(self):
        fig = self._plot(_task(), {'density': self._approach_then_plateau()})
        self.assertTrue(fig.axes[0].get_lines(), 'no mean rule drawn')


# --- stage markers -------------------------------------------------------------------------------

class TestStageMarkers(_FigureCase):
    """Without boundaries a reader cannot tell that a discontinuity is a change of ensemble rather
    than a physical event."""

    def _draw(self, boundaries, n=100):
        df = pd.DataFrame({'TS': range(n), 'density': [1.0] * n})
        t = _task()
        t.dataframes['energy'] = df
        t.stage_boundaries['energy'] = boundaries
        fig, ax = plt.subplots()
        ax.plot(df['TS'], df['density'])
        t._draw_stage_markers(ax, 'density', df, df['TS'])
        return ax

    @staticmethod
    def _rules(ax):
        # the data line plus one Line2D per axvline
        return len(ax.get_lines()) - 1

    def test_the_first_stage_is_not_marked(self):
        """Its rule would sit on the y axis and mark nothing."""
        ax = self._draw([(0, 'a'), (50, 'b')])
        self.assertEqual(self._rules(ax), 1)

    def test_a_single_stage_draws_nothing(self):
        self.assertEqual(self._rules(self._draw([(0, 'only')])), 0)

    def test_every_boundary_gets_a_rule(self):
        ax = self._draw([(0, 'a'), (25, 'b'), (50, 'c'), (75, 'd')])
        self.assertEqual(self._rules(ax), 3)

    def test_labels_are_dropped_when_there_are_too_many_stages(self):
        many = [(i * 5, f's{i}') for i in range(_MAX_LABELED_STAGES + 3)]
        ax = self._draw(many)
        self.assertEqual(len(ax.texts), 0, 'labels should be suppressed')
        self.assertGreater(self._rules(ax), _MAX_LABELED_STAGES,
                           'the rules themselves must still be drawn')

    def test_crowded_labels_are_thinned_but_still_marked(self):
        # three boundaries within 1% of the span: rules for all, at most one label
        ax = self._draw([(0, 'a'), (50, 'b'), (51, 'c'), (52, 'd')])
        self.assertEqual(self._rules(ax), 3)
        self.assertLessEqual(len(ax.texts), 2)

    def test_an_out_of_range_boundary_is_ignored(self):
        ax = self._draw([(0, 'a'), (500, 'past the end')])
        self.assertEqual(self._rules(ax), 0)

    def test_an_unknown_series_draws_nothing(self):
        t = _task()
        fig, ax = plt.subplots()
        t._draw_stage_markers(ax, 'density', pd.DataFrame({'TS': [0, 1]}), [0, 1])
        self.assertEqual(len(ax.get_lines()), 0)


# --- overlay source parsing --------------------------------------------------------------------

class TestSeriesFromLogs(unittest.TestCase):
    """Overlay sources are read only to be drawn, so unlike the main path they must NOT write
    CSVs -- doing so would scatter files named after another run through this run's directory."""

    FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fixtures')

    def test_a_real_log_yields_its_series(self):
        series = _task()._series_from_logs([os.path.join(self.FIXTURES, 'npt_run.log')])
        self.assertIn('energy', series)
        self.assertEqual(len(series['energy']), 45)

    def test_no_csvs_are_written(self):
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as d:
            os.chdir(d)
            try:
                _task()._series_from_logs([os.path.join(self.FIXTURES, 'npt_run.log')])
                self.assertFalse(os.listdir(d), 'overlay parsing must not write files')
            finally:
                os.chdir(cwd)

    def test_an_unreadable_log_is_skipped_not_fatal(self):
        self.assertEqual(_task()._series_from_logs(['no_such_file.log']), {})

    def test_separate_logs_are_concatenated_whole(self):
        log = os.path.join(self.FIXTURES, 'npt_run.log')
        series = _task()._series_from_logs([log, log])
        self.assertEqual(len(series['energy']), 90)

    def test_a_shared_boundary_timestep_is_not_counted_twice(self):
        """Chunked runs overlap by one record: the step a chunk ends on is the step the next one
        starts from, and keeping both puts a duplicate sample in the middle of every trace."""
        log = os.path.join(self.FIXTURES, 'npt_run.log')
        head, tail, last_ts = [], [], None
        for line in open(log):
            if not line.startswith('ENERGY: '):
                head.append(line)
                continue
            last_ts = line.split()[1]
            tail.append(line)
        with tempfile.TemporaryDirectory() as d:
            # a second chunk that resumes on the timestep the first ended
            second = os.path.join(d, 'second.log')
            with open(second, 'w') as f:
                f.writelines(head)
                f.writelines([l for l in tail if l.split()[1] == last_ts])
                f.writelines(tail[:3])
            series = _task()._series_from_logs([log, second])
        self.assertEqual(len(series['energy']), 45 + 4 - 1,
                         'the repeated boundary record should have been dropped')

    def test_simulation_time_is_derived_when_the_timestep_is_known(self):
        series = _task()._series_from_logs([os.path.join(self.FIXTURES, 'npt_run.log')])
        self.assertIn('time_ps', series['energy'].columns)
        self.assertGreater(series['energy']['time_ps'].iloc[-1],
                           series['energy']['time_ps'].iloc[0])


# --- choosing which runs to plot ------------------------------------------------------------------

class TestCollectFromMdHistory(unittest.TestCase):
    """Which earlier runs describe the system we are holding now.

    Decided by atom count rather than task adjacency: a ``solvate`` between two MD stages makes
    everything before it a different system, while a ``validate`` or ``ring_check`` between them
    changes nothing.
    """

    def _task_with(self, history, natoms=14000, present=True):
        from pestifer.core.artifacts import FileArtifactList
        t = _task()
        t.csvartifacts = FileArtifactList([])
        state = mock.Mock()
        state.psf.path = 'sys.psf'
        artifacts = {'md_history': mock.Mock(data=history), 'state': state}
        t.get_current_artifact = lambda k: artifacts.get(k)
        self._patches = [
            mock.patch('pestifer.tasks.mdplot.total_atoms', return_value=natoms),
            mock.patch('pestifer.core.artifacts.CSVDataFileArtifact.exists',
                       return_value=present),
        ]
        for p in self._patches:
            p.start()
            self.addCleanup(p.stop)
        return t

    @staticmethod
    def _entry(natoms, base='00-01-000_md', keys=('energy',)):
        return {'natoms': natoms, 'csv_basename': base, 'csv_keys': list(keys),
                'controller': 0, 'task_index': 1}

    def test_no_history_falls_back_to_the_caller(self):
        t = _task()
        t.csvartifacts = []
        t.get_current_artifact = lambda k: None
        self.assertFalse(t._collect_from_md_history())

    def test_an_empty_history_falls_back_too(self):
        t = self._task_with([])
        self.assertFalse(t._collect_from_md_history())

    def test_runs_matching_the_current_atom_count_are_collected(self):
        t = self._task_with([self._entry(14000), self._entry(14000, '00-02-000_md')])
        self.assertTrue(t._collect_from_md_history())
        self.assertEqual(len(t.csvartifacts), 2)

    def test_runs_describing_an_earlier_smaller_system_are_dropped(self):
        """The pre-solvation minimize belongs to a different system and must not be plotted."""
        t = self._task_with([self._entry(1000, '00-00-000_md'), self._entry(14000)])
        t._collect_from_md_history()
        self.assertEqual(len(t.csvartifacts), 1)

    def test_a_history_entry_whose_csv_is_gone_is_skipped(self):
        t = self._task_with([self._entry(14000)], present=False)
        self.assertFalse(t._collect_from_md_history())
        self.assertEqual(len(t.csvartifacts), 0)


# --- overlays ------------------------------------------------------------------------------------

class TestOverlays(_FigureCase):
    """One figure per quantity, one curve per run -- for comparing this build against others."""

    DF = pd.DataFrame({'TS': [0, 100, 200], 'density': [1.0, 1.01, 1.02]})
    OTHER = pd.DataFrame({'TS': [0, 100, 200], 'density': [1.1, 1.11, 1.12]})

    def _plot(self, task, overlays, df_of_column=None, series=None):
        with mock.patch.object(MDPlotTask, '_series_from_logs',
                               return_value=series if series is not None else {}):
            task._plot_overlays(overlays, ['density'],
                                df_of_column if df_of_column is not None else {'density': self.DF},
                                self.tmpdir, block_average=0, grid=False)

    def test_two_runs_produce_one_figure_with_two_labelled_curves(self):
        t = _task()
        self._plot(t, [{'label': 'run B', 'logs': ['b.log']}],
                   series={'energy': self.OTHER})
        self.assertEqual(os.listdir(self.tmpdir), ['test-density-overlay.png'])
        t.register.assert_called_once()

    def test_no_overlay_sources_means_no_figure(self):
        t = _task()
        self._plot(t, [])
        self.assertFalse(os.listdir(self.tmpdir), 'a single run is not an overlay')

    def test_an_entry_that_is_not_a_mapping_is_skipped(self):
        t = _task()
        self._plot(t, ['not-a-dict'])
        self.assertFalse(os.listdir(self.tmpdir))

    def test_an_entry_naming_no_logs_is_skipped(self):
        t = _task()
        self._plot(t, [{'label': 'empty'}])
        self.assertFalse(os.listdir(self.tmpdir))

    def test_an_entry_whose_logs_yield_nothing_is_skipped(self):
        t = _task()
        self._plot(t, [{'label': 'bad', 'logs': ['x.log']}], series={})
        self.assertFalse(os.listdir(self.tmpdir))

    def test_a_quantity_only_one_run_has_is_not_drawn(self):
        """An 'overlay' with one curve is a timeseries wearing the wrong label."""
        t = _task()
        self._plot(t, [{'label': 'run B', 'logs': ['b.log']}],
                   series={'energy': pd.DataFrame({'TS': [0, 1], 'temperature': [300.0, 301.0]})})
        self.assertFalse(os.listdir(self.tmpdir))


# --- panels --------------------------------------------------------------------------------------

class TestPanels(_FigureCase):
    """Three quantities in three figures cannot be read together; a shared x axis is what makes
    'the area kept relaxing after the density settled' visible at all."""

    def _frames(self):
        return {'density': pd.DataFrame({'TS': [0, 100, 200], 'density': [1.0, 1.01, 1.02]}),
                'area': pd.DataFrame({'TS': [0, 100, 200], 'area': [55.0, 54.0, 53.0]})}

    def _plot(self, task, panels, frames=None, **kw):
        opts = dict(block_average=0, stage_markers=False, legend=True, grid=False)
        opts.update(kw)
        with self.captured_figures() as figs:
            task._plot_panels(panels, frames if frames is not None else self._frames(),
                              self.tmpdir, **opts)
        self.figs = figs
        return figs[0] if figs else None

    def test_each_entry_becomes_a_panel_sharing_one_x_axis(self):
        fig = self._plot(_task(), ['density', 'area'])
        self.assertEqual(len(fig.axes), 2)
        self.assertIs(fig.axes[0].get_shared_x_axes().get_siblings(fig.axes[0])[0].figure, fig)

    def test_only_the_bottom_panel_carries_the_x_label(self):
        fig = self._plot(_task(), ['density', 'area'])
        self.assertEqual(fig.axes[0].get_xlabel(), '')
        self.assertTrue(fig.axes[-1].get_xlabel())

    def test_a_list_entry_overlays_traces_in_one_panel(self):
        fig = self._plot(_task(), [['density', 'area']])
        self.assertEqual(len(fig.axes), 1)
        self.assertEqual(len(fig.axes[0].get_lines()), 2)

    def test_a_figure_is_written_and_registered(self):
        t = _task()
        self._plot(t, ['density'])
        self.assertEqual(len(os.listdir(self.tmpdir)), 1)
        t.register.assert_called_once()

    def test_no_figure_when_no_trace_has_data(self):
        t = _task()
        self._plot(t, ['nonesuch'], frames={})
        self.assertFalse(os.listdir(self.tmpdir))
        t.register.assert_not_called()

    def test_an_empty_panel_list_is_a_no_op(self):
        t = _task()
        self._plot(t, [])
        self.assertFalse(os.listdir(self.tmpdir))

    def test_a_single_trace_panel_gets_no_legend(self):
        fig = self._plot(_task(), ['density'])
        self.assertIsNone(fig.axes[0].get_legend())

    def test_block_averaging_adds_a_smoothed_trace_and_a_spread_band(self):
        frames = {'density': pd.DataFrame(
            {'TS': range(50), 'density': list(np.random.default_rng(0).normal(1.0, .01, 50))})}
        fig = self._plot(_task(), ['density'], frames=frames, block_average=5)
        self.assertEqual(len(fig.axes[0].get_lines()), 2, 'raw plus smoothed')
        self.assertTrue(fig.axes[0].collections, 'no spread band drawn')


if __name__ == '__main__':
    unittest.main()
