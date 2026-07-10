# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Base class for all pestifer subcommands.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from argparse import Namespace, ArgumentParser

@dataclass
class Subcommand(ABC):
    name: str
    """ subcommand name """
    aliases: list = field(default_factory=list)
    """ alternative names for this subcommand """
    log_file: str = 'pestifer_diagnostics.log'
    """ default log file for this subcommand """
    short_help: str = ''
    """ help string """
    long_help: str = ''
    """ description """
    func_returns_type: type = bool
    """ the return type of the func method """
    parser: ArgumentParser = None
    """ the subparser for this subcommand """

    @staticmethod
    def func(args: Namespace, **kwargs):
        """ the function that executes the subcommand.  Subclasses must implement this method. """
        return False

    def default_log_file(self, args: Namespace) -> str:
        """
        The default log-file name for this subcommand when the user gives no ``--log-file``.
        The base implementation returns the static :attr:`log_file`; subclasses may override
        to derive a name from the parsed arguments (e.g. from the input config filename).
        """
        return self.log_file

    def __post_init__(self):
        if not self.name:
            raise ValueError("Subcommand name cannot be empty")
        if not self.short_help:
            self.short_help = f"Help for {self.name} subcommand"
        if not self.long_help:
            self.long_help = f"Detailed help for {self.name} subcommand"

    @abstractmethod
    def add_subparser(self, subparsers: ArgumentParser) -> ArgumentParser:
        """
        Adds the subcommand parser to the given subparsers.  Subclasses declare
        arguments here, after a super().add_subparser() call.
        """
        self.parser = subparsers.add_parser(self.name, aliases=self.aliases, help=self.short_help, description=self.long_help)
        self.parser.set_defaults(func=self.func, func_returns_type=self.func_returns_type, default_log_file_func=self.default_log_file)
