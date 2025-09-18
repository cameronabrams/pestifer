# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Subcommands providing access to the examples included in Pestifer.
"""

from dataclasses import dataclass

import logging
logger = logging.getLogger(__name__)

# from pestifer.core.artifacts import FileArtifactList
# from pestifer.core.examplemanager import ExampleManager
from . import Subcommand, RunSubcommand
from ..core.resourcemanager import ResourceManager
# from ..core.example import Example
from ..util.util import remove_argument
from argparse import Namespace

@dataclass
class FetchExampleSubcommand(Subcommand):
    """
    Subcommand to fetch a specific example's YAML configuration file.
    """
    name: str = 'fetch-example'
    func_returns_type: type = str
    short_help: str = "copy the example\'s YAML config file to the CWD"
    long_help: str = "Fetch the YAML configuration file for a specific example by its ID. This command will copy the example's configuration file to the current working directory, allowing you to run simulations or analyses based on that configuration."

    @staticmethod
    def func(args: Namespace, **kwargs):
        example_id = args.example_id
        r = ResourceManager()
        example = r.example_manager.checkout_example(example_id)
        return example.scriptname

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('example_id', type=int, help='The ID of the example to fetch')
        return self.parser

@dataclass
class RunExampleSubcommand(RunSubcommand):
    """
    Subcommand to run a specific example's system preparation.
    """
    name: str = 'run-example'
    short_help: str = "Run a specific example system preparation"
    long_help: str = "Run the system preparation for a specific example by its ID; \'pestifer show-resources examples\' to see the list."
    func_returns_type: type = RunSubcommand.func_returns_type

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        remove_argument(self.parser, 'config')  # Remove the config argument since we will fetch it
        self.parser.add_argument('example_id', type=int, help='The ID of the example to run')
        return self.parser

    @staticmethod
    def func(args: Namespace, **kwargs):
        config = FetchExampleSubcommand.func(args, **kwargs)
        args.config = config
        controller = RunSubcommand.func(args, **kwargs)
        # test the testable artifacts against gold standards if they exist
        return controller