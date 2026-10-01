# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Tests for :mod:`pestifer.util.update_check`.

Nothing here touches the network: every test passes an explicit ``fetcher``, so a run on a
machine with no route out behaves exactly as CI does.  A test that reached PyPI would be
testing PyPI.

The tests that matter most are the negative ones at the bottom -- a check that cannot disturb
the command it precedes is the entire premise of the feature, so it is tested by *breaking* the
check rather than by watching it work.
"""
import io
import json
import os
import unittest

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from pestifer.util.update_check import (CHECK_INTERVAL_SECONDS, ENV_OPT_OUT, _is_newer,
                                        check_for_update, emit_update_notice, is_disabled,
                                        set_enabled)


class _Tty(io.StringIO):
    """A stream that claims to be a terminal, which is what gates the whole check."""

    def isatty(self):
        return True


class UpdateCheckTestCase(unittest.TestCase):

    def setUp(self):
        self._tmp = TemporaryDirectory()
        self.cache = Path(self._tmp.name) / 'update-check.json'
        self.addCleanup(self._tmp.cleanup)
        # A stray opt-out in the developer's own environment would otherwise turn every
        # assertion below into a vacuous pass.
        self._env = patch.dict(os.environ, {ENV_OPT_OUT: ''})
        self._env.start()
        self.addCleanup(self._env.stop)

    def _check(self, current='3.24.2', latest='3.25.0', **kwargs):
        kwargs.setdefault('from_source', False)
        kwargs.setdefault('cache_path', self.cache)
        kwargs.setdefault('fetcher', lambda: latest)
        return check_for_update(current, **kwargs)


class TestVersionComparison(UpdateCheckTestCase):
    """``latest > current``, never ``latest != current``."""

    def test_newer_release_is_an_update(self):
        self.assertTrue(_is_newer('3.25.0', '3.24.2'))
        self.assertTrue(_is_newer('3.24.10', '3.24.2'))  # not a string comparison

    def test_same_version_is_not_an_update(self):
        self.assertFalse(_is_newer('3.24.2', '3.24.2'))

    def test_running_ahead_of_pypi_is_not_an_update(self):
        """The state of this repo between a version bump and its tag."""
        self.assertFalse(_is_newer('3.24.2', '3.25.0'))

    def test_unparseable_versions_are_not_an_update(self):
        self.assertFalse(_is_newer('not-a-version', '3.24.2'))
        self.assertFalse(_is_newer('3.25.0', '0.0.0+unknown!'))
        self.assertFalse(_is_newer(None, '3.24.2'))


class TestCheckForUpdate(UpdateCheckTestCase):

    def test_reports_an_available_update(self):
        result = self._check()
        self.assertTrue(result.update_available)
        self.assertEqual(result.origin, 'pypi')
        self.assertIn('3.25.0', result.notice)
        self.assertIn('pip install -U pestifer', result.notice)

    def test_up_to_date_says_nothing(self):
        self.assertIsNone(self._check(latest='3.24.2').notice)

    def test_a_source_tree_is_not_nagged(self):
        result = self._check(from_source=True)
        self.assertEqual(result.origin, 'source-tree')
        self.assertIsNone(result.notice)

    def test_force_overrides_source_tree_suppression(self):
        """`pestifer check-update` asked explicitly, so it gets a real answer."""
        result = self._check(from_source=True, force=True)
        self.assertEqual(result.origin, 'pypi')
        self.assertTrue(result.update_available)


class TestCaching(UpdateCheckTestCase):

    def test_second_check_within_a_day_does_not_fetch(self):
        self._check(now=1000.0)
        calls = []

        def fetcher():
            calls.append(1)
            return '3.25.0'

        result = self._check(now=1000.0 + CHECK_INTERVAL_SECONDS - 1, fetcher=fetcher)
        self.assertEqual(calls, [], 'a fresh cache must not hit the network')
        self.assertEqual(result.origin, 'cache')
        self.assertTrue(result.update_available)

    def test_a_stale_cache_refetches(self):
        self._check(now=1000.0)
        result = self._check(now=1000.0 + CHECK_INTERVAL_SECONDS + 1, latest='3.26.0')
        self.assertEqual(result.origin, 'pypi')
        self.assertEqual(result.latest, '3.26.0')

    def test_force_ignores_a_fresh_cache(self):
        self._check(now=1000.0)
        result = self._check(now=1000.0, latest='3.26.0', force=True)
        self.assertEqual(result.latest, '3.26.0')

    def test_a_failed_fetch_is_cached_too(self):
        """The guard against an offline machine paying the timeout on every invocation.

        Negative caching is the whole reason a 93-build sweep on a node with no route out
        costs one timeout rather than ninety-three, so it is pinned rather than assumed.
        """
        def boom():
            raise OSError('no route to host')

        result = self._check(now=1000.0, fetcher=boom)
        self.assertEqual(result.origin, 'unavailable')
        self.assertIsNone(result.notice)
        self.assertEqual(json.loads(self.cache.read_text())['last_check'], 1000.0)

        calls = []
        self._check(now=1000.0 + 60, fetcher=lambda: calls.append(1) or '3.25.0')
        self.assertEqual(calls, [], 'a failed check must back off like a successful one')

    def test_an_unreadable_cache_is_not_fatal(self):
        self.cache.write_text('{ this is not json')
        self.assertTrue(self._check().update_available)

    def test_a_cache_stamped_in_the_future_refetches(self):
        """A clock that moved backwards must not wedge the check off for good."""
        self._check(now=10_000.0)
        calls = []
        self._check(now=5_000.0, fetcher=lambda: calls.append(1) or '3.25.0')
        self.assertEqual(len(calls), 1)


class TestOptOut(UpdateCheckTestCase):

    def test_environment_variable_disables(self):
        with patch.dict(os.environ, {ENV_OPT_OUT: '1'}):
            self.assertTrue(is_disabled(self.cache))
            self.assertEqual(self._check().origin, 'disabled')

    def test_falsey_environment_values_do_not_disable(self):
        for value in ('', '0', 'false', 'no'):
            with patch.dict(os.environ, {ENV_OPT_OUT: value}):
                self.assertFalse(is_disabled(self.cache), f'{value!r} should not opt out')

    def test_persistent_opt_out(self):
        set_enabled(False, self.cache)
        self.assertTrue(is_disabled(self.cache))
        self.assertEqual(self._check().origin, 'disabled')
        set_enabled(True, self.cache)
        self.assertFalse(is_disabled(self.cache))
        self.assertTrue(self._check().update_available)

    def test_opting_out_preserves_the_rest_of_the_cache(self):
        self._check(now=1000.0)
        set_enabled(False, self.cache)
        self.assertEqual(json.loads(self.cache.read_text())['latest'], '3.25.0')


class TestEmitUpdateNotice(UpdateCheckTestCase):
    """The startup-path guard.  These are the tests the feature exists to survive."""

    def _emit(self, stream, **kwargs):
        kwargs.setdefault('current_version', '3.24.2')
        kwargs.setdefault('from_source', False)
        kwargs.setdefault('cache_path', self.cache)
        kwargs.setdefault('fetcher', lambda: '3.25.0')
        return emit_update_notice(stream, **kwargs)

    def test_writes_to_a_terminal(self):
        stream = _Tty()
        self.assertTrue(self._emit(stream))
        self.assertIn('A newer pestifer is available: 3.25.0', stream.getvalue())

    def test_a_redirected_stream_gets_nothing_and_fetches_nothing(self):
        """A redirected build must stay byte-identical, and must not touch the network.

        `pestifer build x.yaml > run.log 2>&1` is how the sweep and every cluster job run.  A
        notice that appeared there only when the network happened to be up would make two runs
        of the same example produce different logs.
        """
        stream = io.StringIO()  # a plain StringIO.isatty() is False
        calls = []
        self.assertFalse(self._emit(stream, fetcher=lambda: calls.append(1) or '3.25.0'))
        self.assertEqual(stream.getvalue(), '')
        self.assertEqual(calls, [])
        self.assertFalse(self.cache.exists(), 'a non-tty run must not even write the cache')

    def test_disabled_emits_nothing(self):
        stream = _Tty()
        self.assertFalse(self._emit(stream, enabled=False))
        self.assertEqual(stream.getvalue(), '')

    def test_a_raising_fetcher_is_silent(self):
        stream = _Tty()

        def boom():
            raise ConnectionError('DNS went away')

        self.assertFalse(self._emit(stream, fetcher=boom))
        self.assertEqual(stream.getvalue(), '')

    def test_an_internal_error_cannot_escape(self):
        """The negative control for the one promise this module makes.

        `check_for_update` deliberately lets programming errors through -- only network failures
        are caught there -- so the guard that keeps a broken update check from failing a build
        lives in `emit_update_notice`, and that is the level this is aimed at.  Break the guard
        (narrow the `except`, or drop it) and this test goes red; break only the fetcher and it
        would still pass, which is why the two cases above are not enough on their own.
        """
        stream = _Tty()
        with patch('pestifer.util.update_check.check_for_update',
                   side_effect=RuntimeError('a bug in the update check itself')):
            self.assertFalse(emit_update_notice(stream))
        self.assertEqual(stream.getvalue(), '')

    def test_a_stream_with_no_isatty_is_tolerated(self):
        class Bare:
            pass

        self.assertFalse(self._emit(Bare()))


class TestCLIWiring(unittest.TestCase):
    """That the notice is actually reachable from the command line, not merely importable."""

    def test_cli_calls_the_guard(self):
        from pestifer.cli import pestifer as cli_module
        self.assertIs(cli_module.emit_update_notice, emit_update_notice)

    def test_check_update_is_registered(self):
        from pestifer.subcommands import _subcommands
        names = [s.name for s in _subcommands]
        self.assertIn('check-update', names)

    def test_check_update_suppresses_the_passive_notice(self):
        """Otherwise `pestifer check-update` answers the question twice."""
        import argparse as ap
        from pestifer.subcommands.check_update import CheckUpdateSubcommand
        parser = ap.ArgumentParser()
        subparsers = parser.add_subparsers(dest='command')
        CheckUpdateSubcommand().add_subparser(subparsers)
        args = parser.parse_args(['check-update'])
        self.assertTrue(getattr(args, 'suppress_update_notice', False))


if __name__ == '__main__':
    unittest.main()
