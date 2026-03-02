# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The fetch-example subcommand.  Copies a specific example's YAML configuration file to the
current working directory.
"""
from dataclasses import dataclass

from argparse import Namespace

from ..cli.subcommand import Subcommand

from ..core.resourcemanager import ResourceManager

@dataclass
class FetchExampleSubcommand(Subcommand):
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
        self.parser.add_argument('example_id', type=int, help='the ID of the example to fetch')
        return self.parser
