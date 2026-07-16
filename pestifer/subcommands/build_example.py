# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The build-example subcommand.  Runs the system preparation for a specific example by its ID.
"""
import copy
import logging

from dataclasses import dataclass, field

from argparse import Namespace

from .build import RunSubcommand
from .fetch_example import FetchExampleSubcommand

from ..core.resourcemanager import ResourceManager
from ..util.util import remove_argument

logger = logging.getLogger(__name__)

@dataclass
class RunExampleSubcommand(RunSubcommand):
    name: str = 'build-example'
    aliases: list = field(default_factory=lambda: ['run-example'])
    short_help: str = "build a specific example system"
    long_help: str = "Run the system preparation for a specific example by its ID; \'pestifer show-resources examples\' to see the list."
    func_returns_type: type = RunSubcommand.func_returns_type

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        remove_argument(self.parser, 'config')  # Remove the config argument since we will fetch it
        self.parser.add_argument('example_id', type=int, help='the ID of the example to run')
        return self.parser

    @staticmethod
    def func(args: Namespace, **kwargs):
        config = FetchExampleSubcommand.func(args, **kwargs)
        # Multi-step examples ship helper configs under inputs/aux/; checkout_example
        # flattens them into the CWD alongside the main script.  Run them first, in
        # sorted order (which is their dependency order, e.g. helper-01-base ->
        # helper-02-fabdonor -> helper-03-position-fab{1,2,3}), so their outputs
        # (base.psf, fab_p*.psf, ...) exist before the main script consumes them.
        em = ResourceManager().example_manager
        example = em.examples.get_example_by_example_id(args.example_id)
        aux_dir = em.auxpath(example)
        aux_scripts = sorted(p.name for p in aux_dir.glob('*.yaml')) if aux_dir.is_dir() else []
        for aux in aux_scripts:
            logger.info(f'Example {args.example_id}: running auxiliary helper script {aux}')
            aux_args = copy.copy(args)
            aux_args.config = aux
            RunSubcommand.func(aux_args, **kwargs)
        args.config = config
        controller = RunSubcommand.func(args, **kwargs)
        # test the testable artifacts against gold standards if they exist
        return controller
