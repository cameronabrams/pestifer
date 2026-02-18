# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Implementation of the ``build`` subcommand for launching system preparations.
"""

import datetime
import logging
import os
import shutil
import time
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

from ..core.controller import Controller
from ..core.config import Config
from .subcommand import Subcommand
from ..util.stringthings import __pestifer_version__
from ..util.util import hmsf

@dataclass
class RunSubcommand(Subcommand):
    name: str = 'build'
    aliases: list = field(default_factory=lambda: ['run'])
    long_help: str = 'Prepare a system according to instructions in the config file.'
    short_help: str = 'prepare a system'
    func_returns_type: type = Controller

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('config', type=str, default=None, help='input configuration file in YAML format')
        self.parser.add_argument('--output-dir', type=str, default='./', help='name of output directory relative to CWD (default: %(default)s)')
        self.parser.add_argument('--gpu', default=False, action='store_true', help='force run on GPU')
        self.parser.add_argument('--complete-config', default=False, action='store_true', help='write complete config file')
        return self.parser
    
    @staticmethod
    def func(args, **kwargs):
        """
        Run the ``pestifer build`` command with the specified configuration file.
        """
        if args.output_dir != './':
            exec_dir = os.getcwd()
            if not os.path.exists(args.output_dir):
                os.mkdir(args.output_dir)
            if not os.path.exists(os.path.join(args.output_dir, args.config)):
                shutil.copy(args.config, args.output_dir)
            os.chdir(args.output_dir)
        # Set up the Controller and execute tasks
        begin_time = time.time()
        configname = args.config
        # include date and time in the message below
        logger.info(f'pestifer v. {__pestifer_version__} begins at {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(begin_time))} with config {configname}')

        allowed_extensions = ['.yaml', '.yml', '.y']
        cbase, cext = os.path.splitext(configname)
        if not cext:
            fil = [os.path.exists(f'{cbase}{ext}') for ext in allowed_extensions]
            if any(fil):
                iix = fil.index(True)
                configname = f'{cbase}{allowed_extensions[iix]}'

        config = Config(userfile=configname, **kwargs).configure_new()
        C = Controller().configure(config)
        logger.debug(f'NAMD global config: {C.tasks[0].provisions["namd_global_config"]}')
        if args.gpu:
            if not 'namd' in C.config['user']:
                C.config['user']['namd'] = {}
            C.config['user']['namd']['processor-type'] = 'gpu'

        if args.complete_config:
            C.write_complete_config(f'{cbase}-complete.yaml')
        report = C.do_tasks()
        end_time = time.time()
        elapsed_time_s = datetime.timedelta(seconds=(end_time - begin_time))
        logger.info(f'pestifer ends. Elapsed time {time.strftime("%H:%M:%S", time.gmtime(elapsed_time_s.seconds))}.')
        logger.info(f'Task durations:')
        maxnamelen = max([len(t.taskname) for t in C.tasks]) if len(C.tasks) > 0 else 0
        name_format = f'{{:>{maxnamelen}s}}'
        for task in report.values():
            logger.info(f' - {name_format.format(task["taskname"])}: {hmsf(task["duration"])} ({task["duration_frac"]:>5.1%})')
        if args.output_dir != './':
            os.chdir(exec_dir)
        return C
