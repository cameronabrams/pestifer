# Author: Cameron F. Abrams, <cfa22@drexel.edu>

from abc import ABC, abstractmethod
from dataclasses import dataclass
from argparse import Namespace, ArgumentParser

@dataclass
class Subcommand(ABC):
    name: str
    short_help: str = ''
    long_help: str = ''

    @staticmethod
    def func(args: Namespace, **kwargs):
        pass

    @abstractmethod
    def add_subparser(self, subparsers: ArgumentParser) -> ArgumentParser:
        self.parser = subparsers.add_parser(self.name, help=self.short_help, description=self.long_help)
        self.parser.set_defaults(func=self.func)
