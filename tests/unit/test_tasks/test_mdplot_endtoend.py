# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""End-to-end tests for the ``mdplot`` task, driven the way ``pestifer mdplot`` drives it.

The other two mdplot test files exercise the pieces -- data preparation, derived quantities, each
figure builder in isolation.  This one runs ``MDPlotTask.do()`` whole, through a real ``Config``
and ``Controller``, because ``do()`` is where those pieces are wired together and none of the unit
tests can see the wiring: which series feed which figures, whether stage boundaries survive
concatenation, whether the derived cell columns are available by the time the traces are drawn.

It needs no VMD and no NAMD -- only a log file -- so it belongs in the unit suite despite being
integration-shaped.  The fixtures are a **real** NAMD log and its ``.xst``, trimmed to 45 energy
records.  Written by hand they would be a guess at NAMD's output format, and a parser regression
would sail through; taken from a real run they cannot be subtly wrong about the thing being parsed.

The log is a chunk of a ``density_equilibrate`` run from example 12, which is why it starts at
timestep 19,700 rather than zero and why ``success()`` is False for it -- it is a chunk that was
continued, not a run that ended.
"""

import os
import shutil
import tempfile
import unittest

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

from pestifer.core.config import Config
from pestifer.core.controller import Controller

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fixtures')
LOG = os.path.join(FIXTURES, 'npt_run.log')
XST = os.path.join(FIXTURES, 'npt_run.xst')


class _ReplotCase(unittest.TestCase):
    """Runs mdplot over the fixture log in a scratch directory."""

    @classmethod
    def setUpClass(cls):
        # Config provisions the force field; doing it once keeps this file well under a second
        # per test rather than a second per test.
        cls.config = Config().configure()

    def setUp(self):
        self._cwd = os.getcwd()
        self._tmp = tempfile.TemporaryDirectory()
        os.chdir(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)
        self.addCleanup(os.chdir, self._cwd)
        shutil.copy(LOG, 'run.log')
        shutil.copy(XST, 'run.xst')
        # mdplot ends its figures with plt.clf(), which clears but does not deregister them, so a
        # test that re-plots many times in one process accumulates figures and trips matplotlib's
        # max-open warning.  One figure per invocation is harmless in a build; many in one process
        # is a test artifact, so the tests clean up after themselves.
        self.addCleanup(plt.close, 'all')

    def replot(self, **specs):
        base = {'basename': 'plt', 'reprocess-logs': True, 'logs': ['run.log'],
                'legend': True, 'grid': True, 'units': {'density': 'g/cc'}}
        base.update(specs)
        C = Controller().configure(self.config, userspecs={'tasks': [{'mdplot': base}]},
                                   terminate=False, validate=False)
        C.do_tasks()
        self.task = C.tasks[0]
        return sorted(os.listdir('mdplots')) if os.path.isdir('mdplots') else []


class TestFiguresAreProduced(_ReplotCase):

    def test_a_timeseries_figure_is_written(self):
        self.assertIn('plt-density.png', self.replot(timeseries=['density']))

    def test_each_requested_output_appears(self):
        got = self.replot(timeseries=['density'], histograms=['density'],
                          panels=['density'], **{'summary-table': True})
        self.assertEqual(got, ['plt-density-hist.png', 'plt-density.png',
                               'plt-panels-density.png', 'plt-summary.csv'])

    def test_co_plotted_quantities_share_one_figure(self):
        got = self.replot(timeseries=[['TEMP', 'PRESSURE']])
        self.assertEqual(got, ['plt-TEMP-PRESSURE.png'])

    def test_nothing_is_produced_for_a_quantity_the_run_does_not_have(self):
        with self.assertLogs('pestifer.tasks.mdplot', level='WARNING'):
            got = self.replot(timeseries=['nonesuch'])
        self.assertEqual(got, [], 'a blank figure must not be written')


class TestDerivedCellColumnsReachTheFigures(_ReplotCase):
    """The wiring the unit tests cannot see: ``_add_derived_columns`` runs on the xst series
    during ``do()``, so area/volume/apl must be plottable as ordinary traces."""

    def test_cell_edges_from_the_xst_are_plottable(self):
        self.assertIn('plt-a_x-b_y-c_z.png', self.replot(timeseries=[['a_x', 'b_y', 'c_z']]))

    def test_area_and_volume_are_available_as_traces(self):
        self.assertEqual(self.replot(timeseries=['area', 'volume']),
                         ['plt-area.png', 'plt-volume.png'])

    def test_area_per_lipid_appears_only_when_a_lipid_count_is_given(self):
        self.assertEqual(self.replot(timeseries=['apl']), [])
        self.assertEqual(self.replot(timeseries=['apl'], **{'lipids-per-leaflet': 100}),
                         ['plt-apl.png'])


class TestSummaryTableContent(_ReplotCase):

    def _table(self, **specs):
        self.replot(**{'summary-table': True, **specs})
        return pd.read_csv(os.path.join('mdplots', 'plt-summary.csv'))

    def test_the_table_reports_the_quantity_and_its_units(self):
        table = self._table(timeseries=['density'])
        self.assertEqual(list(table['quantity'].unique()), ['density'])
        self.assertEqual(list(table['units'].unique()), ['g/cc'])

    def test_the_density_is_converted_not_reported_raw(self):
        # raw NAMD density here is ~0.62 amu/A^3; in g/cc that is ~1.03
        mean = self._table(timeseries=['density'])['mean'].iloc[0]
        self.assertAlmostEqual(mean, 1.03, places=1)

    def test_the_step_range_matches_the_log(self):
        table = self._table(timeseries=['density'])
        self.assertEqual(int(table['first_step'].iloc[0]), 19700)
        self.assertEqual(int(table['last_step'].iloc[-1]), 24100)

    def test_the_sample_count_matches_the_energy_records(self):
        self.assertEqual(int(self._table(timeseries=['density'])['samples'].sum()), 45)


class TestProvenanceOfReplottedFigures(_ReplotCase):
    """A re-plot is marked with the seed of the run that produced the data, read from the log --
    not the seed in the configuration driving the re-plot."""

    def test_the_mark_carries_the_logged_seed(self):
        self.replot(timeseries=['density'])
        self.assertIn('1210520238', self.task.build_stamp())

    def test_the_configurations_own_seed_is_not_claimed(self):
        self.replot(timeseries=['density'])
        self.assertNotIn('27021972', self.task.build_stamp())


class TestCsvArtifacts(_ReplotCase):
    """Re-plotting writes the parsed series out beside the figures, which is what makes the data
    behind a figure inspectable."""

    def test_csvs_are_written_for_each_series(self):
        self.replot(timeseries=['density'])
        csvs = sorted(f for f in os.listdir('.') if f.endswith('.csv'))
        self.assertTrue(any('energy' in c for c in csvs), csvs)
        self.assertTrue(any('xst' in c for c in csvs), csvs)

    def test_the_csv_matches_the_log(self):
        self.replot(timeseries=['density'])
        energy = next(f for f in os.listdir('.') if f.endswith('-energy.csv'))
        self.assertEqual(len(pd.read_csv(energy)), 45)


class TestPressureProfiles(unittest.TestCase):
    """The lateral pressure profile: a quantity plotted along the box normal rather than against
    time, and the one figure a bilayer paper usually wants.

    It needs three things to agree -- a ``pressureprofile`` series, an ``.xst`` carrying ``c_z``,
    and matching timesteps between them -- so the fixture is a real NPgT membrane run trimmed to
    twelve records rather than anything synthesized.
    """

    @classmethod
    def setUpClass(cls):
        cls.config = Config().configure()

    def setUp(self):
        self._cwd = os.getcwd()
        self._tmp = tempfile.TemporaryDirectory()
        os.chdir(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)
        self.addCleanup(os.chdir, self._cwd)
        shutil.copy(os.path.join(FIXTURES, 'npgt_profile.log'), 'run.log')
        shutil.copy(os.path.join(FIXTURES, 'npgt_profile.xst'), 'run.xst')
        self.addCleanup(plt.close, 'all')

    def replot(self, **specs):
        base = {'basename': 'plt', 'reprocess-logs': True, 'logs': ['run.log'],
                'legend': True, 'grid': True}
        base.update(specs)
        C = Controller().configure(self.config, userspecs={'tasks': [{'mdplot': base}]},
                                   terminate=False, validate=False)
        C.do_tasks()
        return sorted(os.listdir('mdplots')) if os.path.isdir('mdplots') else []

    def test_a_pressure_profile_figure_is_produced(self):
        self.assertIn('plt-pressureprofile.png', self.replot(profiles=['pressure']))

    def test_the_block_size_is_configurable(self):
        # fewer profiles per block means more curves on the figure, not a different figure
        self.assertIn('plt-pressureprofile.png',
                      self.replot(profiles=['pressure'], **{'profiles-per-block': 4}))

    def test_pressure_units_are_accepted(self):
        for unit in ('bar', 'Pa', '*'):
            with self.subTest(unit=unit):
                self.assertIn('plt-pressureprofile.png',
                              self.replot(profiles=['pressure'], units={'pressure': unit}))

    def test_an_unknown_profile_is_ignored(self):
        self.assertEqual(self.replot(profiles=['nosuchprofile']), [])

    def test_profiles_coexist_with_timeseries(self):
        got = self.replot(profiles=['pressure'], timeseries=['density'])
        self.assertIn('plt-pressureprofile.png', got)
        self.assertIn('plt-density.png', got)


class TestBadInput(_ReplotCase):

    def test_a_missing_log_is_reported_rather_than_silently_plotting_nothing(self):
        C = Controller().configure(
            self.config,
            userspecs={'tasks': [{'mdplot': {'basename': 'plt', 'reprocess-logs': True,
                                             'logs': ['no_such.log'], 'timeseries': ['density']}}]},
            terminate=False, validate=False)
        with self.assertRaises(Exception):
            C.do_tasks()


if __name__ == '__main__':
    unittest.main()
