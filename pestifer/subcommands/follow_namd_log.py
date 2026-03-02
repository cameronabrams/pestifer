# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The follow-namd-log subcommand.  Allows for real-time following and parsing of NAMD log files
actively being written to by a namd3 execution.
"""
from dataclasses import dataclass
import argparse as ap

from ..cli.subcommand import Subcommand

from ..logparsers.namdlogparser import subcommand_follow_namd_log

@dataclass
class FollowNAMDLogSubcommand(Subcommand):
    name: str = "follow-namd-log"
    short_help: str = "follow and parse an actively updating NAMD log file"
    long_help: str = "Monitor a NAMD log file for changes and display relevant information."

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        log=args.log
        basename=args.basename
        subcommand_follow_namd_log(log, basename=basename)
        print()
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('log', type=str, help='NAMD log file to follow')
        self.parser.add_argument('--basename', type=str, default=None, help='base name for output files (default: derived from log file name)')
        return self.parser
