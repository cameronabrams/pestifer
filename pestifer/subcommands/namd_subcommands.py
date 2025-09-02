# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Subcommand that allows for real-time following and parsing of NAMD log files actively being written to by a namd3 execution.
"""
from dataclasses import dataclass
import logging
import argparse as ap
from ..logparsers.namdlogparser import subcommand_follow_namd_log
from ..util.namdrestart import make_namd_restart_subcommand

from . import Subcommand

@dataclass
class FollowNAMDLogSubcommand(Subcommand):
    name: str = "follow-namd-log"
    short_help: str = "Follow an actively updating NAMD log file"
    long_help: str = "Monitor a NAMD log file for changes and display relevant information."

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        log=args.log
        basename=args.basename
        console=logging.StreamHandler()
        loglevel_numeric=getattr(logging,args.log_level.upper())
        console.setLevel(loglevel_numeric)
        formatter=logging.Formatter('%(levelname)s> %(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)
        subcommand_follow_namd_log(log, basename=basename)
        print()
        return True
    
    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('log', type=str, help='NAMD log file to follow')
        self.parser.add_argument('--basename', type=str, default=None, help='Base name for output files (default: derived from log file name)')
        return self.parser

@dataclass
class MakeNAMDRestartSubcommand(Subcommand):
    """
    Subcommand that generates a default NAMD restart configuration based on a current configuration and log file.
    """
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