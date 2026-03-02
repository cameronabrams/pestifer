# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The make-namd-restart subcommand.  Generates a default NAMD restart configuration based on a
current configuration and log file.
"""
from dataclasses import dataclass
import argparse as ap

from ..cli.subcommand import Subcommand

from ..util.namdrestart import make_namd_restart_subcommand

@dataclass
class MakeNAMDRestartSubcommand(Subcommand):
    name: str = "make-namd-restart"
    short_help: str = 'generate a restart NAMD config file based on current checkpoint'
    long_help: str = 'This command generates a NAMD configuration file for restarting a simulation from the current checkpoint.'

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        make_namd_restart_subcommand(args)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('--namd-log', type=str, help='name of most recent NAMD log')
        self.parser.add_argument('--config', type=str, help='name of most recent NAMD config')
        self.parser.add_argument('--new-base', type=str, help='basename of new NAMD config to create (excludes .namd extension)')
        self.parser.add_argument('--run', type=int, help='number of time steps to run')
        self.parser.add_argument('--slurm', type=str, default=None, help='name of SLURM script to update')
        return self.parser
