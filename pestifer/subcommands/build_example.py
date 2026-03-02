# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The build-example subcommand.  Runs the system preparation for a specific example by its ID.
"""
from dataclasses import dataclass, field

from argparse import Namespace

from .build import RunSubcommand
from .fetch_example import FetchExampleSubcommand

from ..util.util import remove_argument

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
        args.config = config
        controller = RunSubcommand.func(args, **kwargs)
        # test the testable artifacts against gold standards if they exist
        return controller
