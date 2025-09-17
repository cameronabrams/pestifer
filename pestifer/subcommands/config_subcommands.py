# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Subcommands pertaining to the pestifer's configuration.
"""
import os
import yaml

import argparse as ap

from dataclasses import dataclass

from . import Subcommand

from ..core.config import Config
from ..core.resourcemanager import ResourceManager

@dataclass
class ConfigHelpSubcommand(Subcommand):
    """
    The config-help subcommand.  This subcommand exposes Ycleptic's built in interactive help system for
    exploring Pestifer's base configuration.  The user can specify the entry point using the
    ``directives`` argument, and can disable interactivity using the ``--no-interactive`` flag.
    """
    name: str = 'config-help'
    short_help: str = "Show help for configuration options"
    long_help: str = "Display help information for all available configuration options."

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        config = Config().configure_new()
        directives = args.directives
        iprompt='pestifer-help: ' if args.interactive else ''
        config.console_help(directives, interactive_prompt=iprompt, exit=True)
        return True
    
    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('directives', type=str, nargs='*', default=[], help='specific directives at which the help system begins (default is at the root)')
        self.parser.add_argument('--interactive', default=True, action=ap.BooleanOptionalAction, help='Enable/disable interactive mode to query config options (default: %(default)s)')

@dataclass
class ConfigDefaultSubcommand(Subcommand):
    """
    Subcommand that echoes all or part of a configuration script with default values, depending on the ``directives`` entry point.
    """
    name: str = 'config-default'
    short_help: str = "Show default configuration options"
    long_help: str = "Display default values for all available configuration options."

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        config = Config().configure_new()
        directives = args.directives
        specs=config.make_default_specs(*directives)
        print(yaml.dump(specs))
        return True
    
    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('directives', type=str, nargs='*', default=[], help='specific directives to get default values for (default: all)')

@dataclass
class NewSystemSubcommand(Subcommand):
    """
    Subcommand that generates a minimal but specific system script.  The single required argument ``id`` is the database ID of the structure file
    that you want to base the system on.  It can be a PDB ID or a UniProt ID, the latter of which is interpreted as an AlphaFold entry.
    The "--full" flag can be used to include a full set of default relaxation/solvation/equilibration protocols.
    """
    name: str = 'new-system'
    short_help: str = 'create a new system script from scratch'
    long_help: str = 'Generate a new system script with a basic template.'

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        r = ResourceManager()
        build_type = 'full' if args.full else 'minimal'
        outputfilename = args.output if args.output else f'{args.id}.yaml'
        title = args.title if args.title else ''
        if os.path.exists(outputfilename):
            raise FileExistsError(f'Output file {outputfilename} already exists; please remove it or choose a different name')
        r.example_manager.new_example_yaml(db_id=args.id, build_type=build_type, outputfilename=outputfilename, title=title)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('id', type=str, default=None, help='PDB/UniProt ID  of the new system to create')
        self.parser.add_argument('--full', default=False, action=ap.BooleanOptionalAction, help='Enable/disable full set of tasks in the new system configuration (default: %(default)s)')
        self.parser.add_argument('--output', type=str, default=None, help='Output filename for the new system configuration (default: <id>.yaml)')
        self.parser.add_argument('--title', type=str, default=None, help='Title for the new system configuration (default: none)')
        return self.parser

@dataclass
class ShowResourcesSubcommand(Subcommand):
    """
    Subcommand that shows available resources for the current configuration.  The single positional argument is the resource type to show:

    1. **tcl**: shows a brief explanation of the TcL resources included in Pestifer.
    2. **examples**: shows a list of available example configurations.
    3. **charmmff**: shows selected information about the CHARMM force field.

    The ``charmff`` type has a few optional switches as well:

    * **toppar**: shows information about the CHARMM topology and parameter files.
    * **custom**: shows information about any custom CHARMM force fields being used.
    * **pdb**: shows information about the PDB files in the auxiliary PDB repository.  
      This repository is a collection of sample PDB files of small molecules and lipids useful 
      for packing.
    * **user-pdbcollection**: allows user to specify paths of their own PDB collections.

    """
    name: str = 'show-resources'
    short_help: str = 'show available resources for the current configuration'
    long_help: str = 'Display a list of all available resources for the current configuration.'

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        r = ResourceManager()
        resource_type = args.resource_type
        specs = {resource_type: []}
        if hasattr(args, resource_type):
            for subresource in getattr(args, resource_type):
                specs[resource_type].append(subresource)
        if args.user_pdbcollection:
            r.update_charmmff(user_pdbrepository_paths=args.user_pdbcollection)
        r.show(out_stream=print, components=specs, fullnames=args.fullnames, missing_fullnames=r.labels.residue_fullnames)
        return True
    
    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('resource_type', type=str, default='examples', help='Type of resource to show; [tcl|examples|charmmff]')
        self.parser.add_argument('--charmmff', type=str, nargs='+', default=[], help='show sub-resources of charmmff resources (\'toppar\', \'custom\', \'pdb\')')
        self.parser.add_argument('--fullnames', default=False, action='store_true', help='show full names of any residues shown with --charmmff pdb')
        self.parser.add_argument('--user-pdbcollection', type=str, nargs='+', default=[], help='additional collections of PDB files outside pestifer installation')
        return self.parser

@dataclass
class WhereTCLSubcommand(Subcommand):
    """
    Subcommand that exposes true locations of various TcL resources in the pestifer installation.
    This is useful for custom VMD scripts that the user may want to create.
    """
    name: str = 'wheretcl'
    short_help: str = "provides path of TcL scripts for sourcing in interactive VMD"
    long_help: str = "Display the path of TcL scripts for sourcing in interactive VMD."

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        r = ResourceManager()
        tcl_root = r.get_tcldir()
        script_dir = r.get_tcl_scriptsdir()
        pkg_dir = r.get_tcl_pkgdir()
        msg = ''
        if args.root or args.verbose:
            assert os.path.exists(tcl_root)
            if args.verbose:
                msg='Tcl root: '
            print(f'{msg}{tcl_root}')
        if args.startup_script_path or args.verbose:
            vmd_startup_script=os.path.join(tcl_root,'vmdrc.tcl')
            assert os.path.exists(vmd_startup_script)
            if args.verbose:
                msg='VMD Startup script: '
            print(f'{msg}{vmd_startup_script}')
        if args.script_dir or args.verbose:
            assert os.path.exists(script_dir)
            if args.verbose:
                msg='VMD Script directory: '
            print(f'{msg}{script_dir}')
        if args.pkg_dir or args.verbose:
            assert os.path.exists(pkg_dir)
            if args.verbose:
                msg='TcL package directory: '
            print(f'{msg}{pkg_dir}')
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('--startup-script-path', default=False, action='store_true', help='print full path of VMD startup script used by pestifer')
        self.parser.add_argument('--verbose', default=False, action='store_true', help='print full paths of all TcL resources of Pestifer')
        self.parser.add_argument('--script-dir', default=False, action='store_true', help='print full path of directory of Pestifer\'s VMD scripts')
        self.parser.add_argument('--pkg-dir', default=False, action='store_true', help='print full path of directory of Pestifer\'s VMD/TcL package library')
        self.parser.add_argument('--root', default=False, action='store_true', help='print full path of Pestifer\'s root TcL directory')
        return self.parser