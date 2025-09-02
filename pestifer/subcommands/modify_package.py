# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Subcommand to modify the pestifer package.  This subcommand is only exposed if the Pestifer installation includes the entire source tree (i.e., it was installed with pip -e from a git repository).
The two main ways you can modify the package are

1. Managing examples, including adding, updating, and deleting example scripts.  Updating can include name-changes, adding new auxiliary inputs, and modifying expected outputs.
2. Regenerating atomselect macros based on globals defined in :mod:`core/labels.py <pestifer.core.labels>`.
"""
from dataclasses import dataclass

import logging

import argparse as ap

from . import Subcommand

from ..core.resourcemanager import ResourceManager

logger = logging.getLogger(__name__)

@dataclass
class ModifyPackageSubcommand(Subcommand):
    name: str = 'modify-package'
    short_help: str = "Modify the pestifer package"
    long_help: str = "Modify the pestifer package by adding, updating, or deleting examples."

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        """
        Modify the pestifer source package.  Current modifications allowed involve example addition/update/deletion, and
        updating atomselect macros based on globals defined in :mod:`core/labels.py <pestifer.core.labels>`.

        .. note::
        
            The command ``pestifer modify-package`` can only be invoked if ``pestifer`` is installed as an editable source tree, i.e. the directory containing the ``pestifer`` package.  A simple pip installation of the pestifer package will not allow this command to run, as the pestifer package will not have a full source tree available.  If you want to modify the package, you must clone the repository from GitHub and then execute ``pip install -e .`` from the package root.  Before modifying, it would also be a good idea to create a new branch, rather than modifying in the main branch.
        """
        # First, verify that the user has a full source tree by verifying the existence of the docs directory.  If not, raise an error.
        logging.basicConfig(level=getattr(logging, args.log_level.upper()))
        RM = ResourceManager()
        match args.example_action:
            case 'add':
                example_id = args.example_id
                scriptname = args.example_scriptname
                author_name = args.example_author_name
                author_email = args.example_author_email
                title = args.example_title
                db_id = args.example_db_id
                auxiliary_inputs = args.example_auxiliary_inputs
                outputs = args.example_outputs
                if scriptname:
                    RM.add_example(scriptname, example_id=example_id, author_name=author_name, author_email=author_email, title=title, db_id=db_id, auxiliary_inputs=auxiliary_inputs, outputs=outputs)
                else:
                    raise ValueError(f'Invalid parameter for add example action: scriptname={scriptname}. Must be non-empty YAML file name.')
            case 'update':
                example_id = args.example_id
                name = args.example_name
                author_name = args.example_author_name
                author_email = args.example_author_email
                title = args.example_title
                db_id = args.example_db_id
                auxiliary_inputs = args.example_auxiliary_inputs
                outputs = args.example_outputs
                if example_id > 0:
                    RM.update_example(example_id, shortname=name, author_name=author_name, author_email=author_email, title=title, db_id=db_id, auxiliary_inputs=auxiliary_inputs, outputs=outputs)
                else:
                    raise ValueError(f'Invalid parameters for update example action: example_id={example_id}. Must be positive integer.')
            case 'delete':
                example_id = args.example_id
                if example_id > 0:
                    RM.delete_example(example_id)
                else:
                    raise ValueError(f'Invalid parameter for delete example action: example_id={example_id}. Must be positive integer.')
            case 'rename':
                example_id = args.example_id
                new_name = args.example_name
                if example_id > 0 and new_name:
                    RM.rename_example(example_id, new_name=new_name)
                else:
                    raise ValueError(f'Invalid parameters for rename example action: example_id={example_id}, new_name={new_name}. Must be positive integer and non-empty YAML file name.')
            case 'author':
                example_id = args.example_id
                author_name = args.example_author_name
                author_email = args.example_author_email
                if example_id > 0 and author_name and author_email:
                    RM.set_example_author(example_id, author_name, author_email)
                else:
                    raise ValueError(f'Invalid parameters for set-author example action: example_id={example_id}, author_name={author_name}, author_email={author_email}. Must be positive integer and non-empty strings.')
            case _:
                raise ValueError(f'Invalid example action: {args.example_action}.')
        if args.update_atomselect_macros:
            RM.update_atomselect_macros()
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('--update-atomselect-macros', action='store_true', help='update the resources/tcl/macros.tcl file based on content in core/labels.py; developer use only')
        self.parser.add_argument('--example-id', type=int, default=0, help='integer ID of example to modify; developer use only')
        self.parser.add_argument('--example-action', type=str, default=None, choices=[None, 'add', 'update', 'delete', 'rename', 'author'], help='action to perform on the example; choices are [add|update|delete|rename|author]; developer use only')
        self.parser.add_argument('--example-scriptname', type=str, default='', help='yaml file of example; developer use only')
        self.parser.add_argument('--example-name', type=str, default='', help='new name for the example for action \'rename\'; Default is basename of the script file (without extension); developer use only')
        self.parser.add_argument('--log-level', type=str, default='info', choices=['info', 'debug', 'warning'], help='Logging level (default: %(default)s)')
        self.parser.add_argument('--example-author-name', type=str, default='', help='Name of the author; if not given, pestifer attempts to extract it from the script header')
        self.parser.add_argument('--example-author-email', type=str, default='', help='Email of the author; if not given, pestifer attempts to extract it from the script header')
        self.parser.add_argument('--example-title', type=str, default='', help='Descriptive 1-line title of the example (default: extract from \'title\' directive in the script)')
        self.parser.add_argument('--example-db-id', type=str, default='', help='Database ID of the example; if not given, pestifer attempts to extract it from the script\'s "fetch" task')
        self.parser.add_argument('--example-auxiliary-inputs', type=str, nargs='*', default=[], help='List of auxiliary input files for the example')
        self.parser.add_argument('--example-outputs', type=str, nargs='*', default=[], help='List of output files for the example')
        return self.parser