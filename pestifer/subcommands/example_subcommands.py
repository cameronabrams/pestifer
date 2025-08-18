# Author: Cameron F. Abrams <cfa22@drexel.edu>

from dataclasses import dataclass
from . import Subcommand, RunSubcommand
from ..core.resourcemanager import ResourceManager
from ..core.controller import Controller
from argparse import Namespace

@dataclass
class FetchExampleSubcommand(Subcommand):
    name: str = 'fetch-example'
    example_id: int = 0
    return_type: type = str
    short_help: str = "copy the example\'s YAML config file to the CWD"
    long_help: str = "Fetch the YAML configuration file for a specific example by its ID. This command will copy the example's configuration file to the current working directory, allowing you to run simulations or analyses based on that configuration."

    @staticmethod
    def func(args: Namespace, **kwargs):
        example_id = args.example_id
        r = ResourceManager()
        config = r.example_manager.checkout_example(example_id)
        return config

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('example_id', type=int, help='The ID of the example to fetch')
        return self.parser

@dataclass
class RunExampleSubcommand(Subcommand):
    name: str = 'run-example'
    example_id: int = 0
    short_help: str = "Run a specific example system preparation"
    long_help: str = "Run the system preparation for a specific example by its ID; \'pestifer show-resources examples\' to see the list."
    return_type: type = Controller

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('example_id', type=int, help='The ID of the example to run')
        return self.parser

    @staticmethod
    def func(args: Namespace, **kwargs):
        config = FetchExampleSubcommand.func(args, **kwargs)
        args.config = config
        return RunSubcommand.func(args, **kwargs)
