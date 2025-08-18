# Author: Cameron F. Abrams <cfa22@drexel.edu>

from dataclasses import dataclass
import logging

from ..logparsers.namdlogparser import subcommand_follow_namd_log
from ..util.namdrestart import make_namd_restart_subcommand

from . import Subcommand

@dataclass
class FollowNAMDLogSubcommand(Subcommand):
    name: str = "follow-namd-log"
    short_help: str = "Follow an actively updating NAMD log file"
    long_help: str = "Monitor a NAMD log file for changes and display relevant information."

    @staticmethod
    def func(args, **kwargs):
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
    
    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('log', type=str, help='NAMD log file to follow')
        self.parser.add_argument('--basename', type=str, default=None, help='Base name for output files (default: derived from log file name)')
        self.parser.add_argument('--diagnostic-log-file', type=str, default=None, help='diagnostic log file')
        self.parser.add_argument('--log-level', type=str, default='info', choices=['info', 'debug', 'warning'], help='Logging level (default: %(default)s)')
        self.parser.add_argument('--no-banner', default=True, action='store_true', help='turn off the banner')
        return self.parser

@dataclass
class MakeNAMDRestartSubcommand(Subcommand):
    name: str = "make-namd-restart"
    short_help: str = 'generate a restart NAMD config file based on current checkpoint'
    long_help: str = 'This command generates a NAMD configuration file for restarting a simulation from the current checkpoint.'

    @staticmethod
    def func(args, **kwargs):
        make_namd_restart_subcommand(args)

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('--namd-log', type=str, help='name of most recent NAMD log')
        self.parser.add_argument('--config', type=str, help='name of most recent NAMD config')
        self.parser.add_argument('--new-base', type=str, help='basename of new NAMD config to create (excludes .namd extension)')
        self.parser.add_argument('--run', type=int, help='number of time steps to run')
        self.parser.add_argument('--slurm', type=str, default=None, help='name of SLURM script to update')
        return self.parser