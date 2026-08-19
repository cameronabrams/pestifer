# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the provenance mark pestifer stamps on the files it authors.

Two properties are being protected here, and they pull in opposite directions.

*Files pestifer writes* should say what produced them.  A figure lifted into a slide, a psfgen
script mailed to a colleague, a ``.dat`` convergence report attached to an issue -- all of these
routinely travel without the build directory that would otherwise explain them.  The mark carries
the version and the build-wide NAMD seed, the latter because a triplicate sweep produces three
near-identical artifacts that nothing else distinguishes.

*Files pestifer copies* must be left exactly alone.  The CHARMM ``.rtf``/``.str``/``.prm`` files
are staged into the build directory verbatim, and byte-identity with the upstream release is
precisely how a reader verifies they were not tampered with.  Stamping those would destroy the
provenance the mark exists to establish, so there is a test for it below.

The mark deliberately replaced a wall-clock ``Created <ctime>`` banner.  A timestamp says nothing
about what produced a file and makes byte-identical regeneration impossible; the build time is
recorded in the log and in ``run-record.json``, where it costs nothing.
"""

import hashlib
import os
import re
import tempfile
import unittest

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pestifer.charmmff.charmmffprm import CharmmParamFile
from pestifer.scripters.tcl import TcLScripter
from pestifer.util.provenance import stamp, stamp_figure
from pestifer.util.stringthings import __pestifer_version__


class TestTheMark(unittest.TestCase):

    def test_the_version_is_always_present(self):
        self.assertIn(__pestifer_version__, stamp())

    def test_the_seed_is_included_when_there_is_one(self):
        self.assertIn('27021973', stamp(27021973))

    def test_a_missing_seed_is_omitted_rather_than_rendered_as_none(self):
        for empty in (None, '', 0):
            self.assertNotIn('None', stamp(empty))
            self.assertEqual(stamp(empty), f'pestifer {__pestifer_version__}')


class TestGeneratedScripts(unittest.TestCase):
    """Every ``.tcl`` and ``.namd`` pestifer writes goes through TcLScripter's header."""

    def setUp(self):
        self._cwd = os.getcwd()
        self._tmp = tempfile.TemporaryDirectory()
        os.chdir(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)
        self.addCleanup(os.chdir, self._cwd)

    @staticmethod
    def _write(name, seed=27021973):
        s = TcLScripter(namd_config={'seed': seed} if seed else None)
        s.newscript(name)
        s.addline('mol new foo.pdb')
        s.writefile()
        return open(f'{name}.tcl').read()

    def test_the_header_carries_version_and_seed(self):
        header = self._write('demo')
        self.assertIn(f'pestifer {__pestifer_version__}', header)
        self.assertIn('27021973', header)

    def test_no_wall_clock_timestamp_is_written(self):
        # a ctime string, e.g. 'Tue Aug 19 11:04:07 2026'
        self.assertNotRegex(self._write('demo'),
                            r'\w{3} \w{3} ?\d{1,2} \d{2}:\d{2}:\d{2} \d{4}')

    def test_repeated_generation_is_byte_identical(self):
        """The property the timestamp used to make impossible.

        This is also what stops test runs from dirtying tracked fixture scripts in the working
        tree, which they did on every invocation while the banner held a clock reading.
        """
        digests = []
        for _ in range(2):
            digests.append(hashlib.sha256(self._write('repro').encode()).hexdigest())
        self.assertEqual(*digests)

    def test_a_scripter_without_a_seed_still_gets_a_mark(self):
        # the standalone paths (subcommands, --check) have no build config to draw a seed from
        header = self._write('noseed', seed=None)
        self.assertIn(f'pestifer {__pestifer_version__}', header)
        self.assertNotIn('None', header)

    def test_kwargs_reach_the_base_scripter(self):
        """Regression: TcLScripter forwarded only comment_char, so build_seed was always None and
        every generated script went out unstamped despite the seed being configured."""
        self.assertEqual(TcLScripter(namd_config={'seed': 42}).build_seed, 42)


class TestTaskStamp(unittest.TestCase):
    """``BaseTask.build_stamp`` is called from inside plotting and reporting code, so it has to
    survive every state a task can be in -- including never having been provisioned, which is how
    the standalone re-plot paths construct one."""

    @staticmethod
    def _task(provisions='unset'):
        from pestifer.tasks.mdplot import MDPlotTask
        t = MDPlotTask.__new__(MDPlotTask)
        if provisions != 'unset':
            t.provisions = provisions
        return t

    def test_an_unprovisioned_task_still_produces_a_mark(self):
        # regression: attribute access raised AttributeError here and took the figure down with it
        self.assertIn(__pestifer_version__, self._task().build_stamp())

    def test_empty_provisions_are_tolerated(self):
        for p in (None, {}, {'namd_global_config': None}, {'namd_global_config': {}}):
            self.assertIn(__pestifer_version__, self._task(p).build_stamp())

    def test_a_configured_seed_reaches_the_mark(self):
        got = self._task({'namd_global_config': {'seed': 27021973}}).build_stamp()
        self.assertIn('27021973', got)


class TestReplotStampsTheDataItPlots(unittest.TestCase):
    """Re-plotting is the one path where the running configuration is not the provenance.

    ``pestifer mdplot`` and ``reprocess-logs: true`` draw figures from logs some other run
    produced.  Stamping those with the current configuration's seed attributes the data to a run
    that never made it -- which is exactly what happened in practice: a figure re-plotted from a
    sweep log came out marked with the schema's default seed, while the log itself recorded
    1609943847.
    """

    @staticmethod
    def _task(log_seeds, config_seed=27021972):
        from pestifer.tasks.mdplot import MDPlotTask
        t = MDPlotTask.__new__(MDPlotTask)
        t.reprocess_logs = True
        t._log_seeds = set(log_seeds)
        t.provisions = {'namd_global_config': {'seed': config_seed}}
        return t

    def test_the_seed_comes_from_the_log_not_the_configuration(self):
        got = self._task([1609943847]).build_stamp()
        self.assertIn('1609943847', got)
        self.assertNotIn('27021972', got, "the current config's seed is not this data's provenance")

    def test_disagreeing_logs_yield_no_seed_at_all(self):
        # a build's chunks each get a derived seed; there is no single right answer, so claim none
        got = self._task([111, 222]).build_stamp()
        self.assertNotIn('seed', got)
        self.assertIn(__pestifer_version__, got)

    def test_logs_without_a_recorded_seed_yield_no_seed(self):
        self.assertNotIn('seed', self._task([None]).build_stamp())
        self.assertNotIn('seed', self._task([]).build_stamp())

    def test_a_normal_build_still_uses_the_build_seed(self):
        from pestifer.tasks.mdplot import MDPlotTask
        t = MDPlotTask.__new__(MDPlotTask)
        t.reprocess_logs = False
        t.provisions = {'namd_global_config': {'seed': 27021973}}
        self.assertIn('27021973', t.build_stamp())


class TestGeneratedParameterFile(unittest.TestCase):
    """``*_minimal.prm`` is generated by pestifer, so it is stamped -- unlike the CHARMM files
    beside it, which are copied."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.path = os.path.join(self._tmp.name, 'min.prm')

    def test_the_mark_is_written_as_a_charmm_title_line(self):
        CharmmParamFile().write(self.path, title='Minimal', stamp=stamp(27021973))
        lines = open(self.path).read().splitlines()
        marked = [l for l in lines if 'pestifer' in l]
        self.assertTrue(marked, 'no mark written')
        # CHARMM title records must begin with '*', or the file will not parse
        for line in marked:
            self.assertTrue(line.startswith('*'), f'not a legal CHARMM title record: {line!r}')

    def test_the_title_block_is_terminated(self):
        CharmmParamFile().write(self.path, title='Minimal', stamp=stamp(27021973))
        lines = [l for l in open(self.path).read().splitlines() if l.startswith('*')]
        self.assertEqual(lines[-1], '*', 'the CHARMM title block must end with a bare *')

    def test_an_unstamped_write_is_still_valid(self):
        CharmmParamFile().write(self.path, title='Minimal')
        self.assertNotIn('pestifer 3', open(self.path).read().split('\n\n')[0])


class TestCopiedForceFieldFilesAreNeverStamped(unittest.TestCase):
    """The hazard a naive 'stamp everything pestifer writes' rule would create.

    Pestifer *writes* ``.rtf``/``.str``/``.prm`` files into the build directory, but it is copying
    upstream CHARMM content.  Those must stay byte-identical to the release, because that is how a
    reader verifies they are unmodified -- so the staging path must remain a pure copy.
    """

    def test_staging_uses_a_byte_for_byte_copy(self):
        import inspect
        from pestifer.charmmff import autocache
        src = inspect.getsource(autocache)
        self.assertIn('shutil.copyfile', src,
                      'CHARMM resources must be staged by copy, not rewritten')
        self.assertNotIn('util.provenance', src,
                         'copied CHARMM files must never be stamped')

    def test_a_copied_file_is_unchanged(self):
        import shutil
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, 'top_all36_prot.rtf')
            body = '* CHARMM36 all-hydrogen topology\n*\n\nRESI ALA 0.00\n'
            open(src, 'w').write(body)
            dst = os.path.join(d, 'staged.rtf')
            shutil.copyfile(src, dst)
            self.assertEqual(open(dst).read(), body)
            self.assertNotIn('pestifer', open(dst).read())


class TestFigureFooter(unittest.TestCase):

    def setUp(self):
        self.addCleanup(plt.close, 'all')

    def test_the_mark_is_drawn_on_the_figure(self):
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        stamp_figure(fig, stamp(27021973))
        texts = [t.get_text() for t in fig.texts]
        self.assertTrue(any(__pestifer_version__ in t for t in texts), texts)
        self.assertTrue(any('27021973' in t for t in texts), texts)

    def test_the_footer_sits_outside_the_axes(self):
        """It must not land on top of the data it is annotating."""
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        stamp_figure(fig, stamp())
        x, y = fig.texts[0].get_position()
        bbox = ax.get_position()
        self.assertLess(y, bbox.y0, 'the mark overlaps the axes')

    def test_an_explicit_text_overrides_the_default(self):
        fig, _ = plt.subplots()
        stamp_figure(fig, 'custom mark')
        self.assertEqual(fig.texts[0].get_text(), 'custom mark')


if __name__ == '__main__':
    unittest.main()
