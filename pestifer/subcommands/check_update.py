# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The check-update subcommand: ask PyPI, now, whether a newer pestifer exists.

The passive check on the startup path is deliberately quiet -- it skips a redirected run and
answers from a day-old cache -- so this is the way to get a definite answer, and the thing to
tell someone who may be running a stale version.  It ignores both the cache interval and the
source-tree suppression, and unlike the passive check it reports *why* it has no answer.
"""
import argparse as ap

from dataclasses import dataclass

from . import Subcommand
from ..util.stringthings import __pestifer_version_from_source__
from ..util.update_check import CHANGELOG_URL, check_for_update, set_enabled

_ORIGIN_EXPLANATION = {
    'unavailable': 'could not reach PyPI (no network, or it is down)',
    'disabled':    'update checking is disabled',
}


@dataclass
class CheckUpdateSubcommand(Subcommand):
    name: str = 'check-update'
    group: str = 'Manage the installation'
    short_help: str = 'report whether a newer pestifer has been released'
    long_help: str = (
        'Asks PyPI for the latest released pestifer and compares it with the one you are '
        'running.  pestifer also checks this on its own, at most once a day and only when it '
        'is attached to a terminal; this command always asks.  Use --disable/--enable to turn '
        'the automatic check off or on for good (equivalently: PESTIFER_NO_UPDATE_CHECK=1, or '
        '--no-update-check for one invocation).'
    )

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        if args.disable is not None:
            set_enabled(not args.disable)
            print(f'Automatic update checking {"disabled" if args.disable else "enabled"}.')
            return True

        result = check_for_update(force=True)
        print(f'Running:  {result.current}'
              f'{" (from a source tree)" if __pestifer_version_from_source__ else ""}')
        if result.latest is None:
            print(f'Latest:   unknown -- {_ORIGIN_EXPLANATION.get(result.origin, result.origin)}')
            return True
        print(f'Latest:   {result.latest} (PyPI)')
        if result.update_available:
            print(f'\nA newer pestifer is available.\n'
                  f'  upgrade:   pip install -U pestifer\n'
                  f'  changelog: {CHANGELOG_URL}')
        elif result.latest == result.current:
            print('\nYou are up to date.')
        else:
            # Newer than the latest release: a working tree between a bump and its tag, or an
            # install from git.  Not a problem, but saying "up to date" would be wrong.
            print('\nYou are running a version newer than the latest release.')
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument(
            '--disable',
            default=None,
            action=ap.BooleanOptionalAction,
            help='turn the automatic check off (--disable) or back on (--no-disable) persistently')
        # Without this, the startup check would print its notice and then this command would
        # print the same thing again.
        self.parser.set_defaults(suppress_update_notice=True)
        return self.parser
