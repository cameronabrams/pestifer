# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The report-methods subcommand: draft a Methods section from completed builds.

Named after ``gmx report-methods``, which does the same job for a GROMACS ``.tpr``.  Reads the
``run-record.json`` a finished build leaves behind -- one, or several for replicas -- and writes
a LaTeX fragment and a BibTeX file.

Deliberately not run automatically at the end of a build.  ``terminate`` records the facts;
turning them into prose is a separate, explicit act, because a Methods paragraph that appears
unbidden invites being pasted unread.
"""
import argparse as ap
import os
import sys

from dataclasses import dataclass, field

from . import Subcommand

from ..core.run_record import RUN_RECORD_NAME, load_run_records
from ..util.methods_report import write_report


@dataclass
class ReportMethodsSubcommand(Subcommand):
    name: str = 'report-methods'
    group: str = 'After the run'
    short_help: str = 'draft a Methods section from one or more completed builds'
    long_help: str = ('Draft a Methods section and bibliography from the run-record.json files '
                      'that completed builds leave behind. Pass several run directories to '
                      'describe replicas as a set.')
    func_returns_type: type = bool

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument(
            'runs', nargs='*', default=['.'],
            help=f'run directories or {RUN_RECORD_NAME} files (default: the current directory). '
                 f'Give several to describe replicas as one set.')
        self.parser.add_argument(
            '--format', dest='formats', action='append', choices=['tex', 'bib', 'md'],
            help='output format; repeatable (default: tex and bib)')
        self.parser.add_argument(
            '--output-dir', default='.', help='where to write the files (default: %(default)s)')

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        try:
            records = load_run_records(args.runs)
        except FileNotFoundError as e:
            print(f'ERROR: {e}\n'
                  f'A build writes {RUN_RECORD_NAME} when it completes; an aborted build has '
                  f'none.', file=sys.stderr)
            return False
        except ValueError as e:
            print(f'ERROR: {e}', file=sys.stderr)
            return False
        if not records:
            print('ERROR: no run records given.', file=sys.stderr)
            return False

        formats = tuple(args.formats or ('tex', 'bib'))
        written = write_report(records, outdir=args.output_dir, formats=formats)
        n = len(records)
        print(f'Drafted a Methods section from {n} run record{"s" if n != 1 else ""}:')
        for p in written:
            print(f'  {p}')
        print('\nThis is a draft. Every number in it is what the build did, but the prose is '
              'machine-generated\nand the TODO markers are deliberate. Read and rewrite it '
              'before it goes anywhere.')
        return True
