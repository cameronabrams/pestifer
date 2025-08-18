# Author: Cameron F. Abrams <cfa22@drexel.edu>
import os
import yaml

from dataclasses import dataclass

from . import Subcommand

from ..core.config import Config
from ..core.resourcemanager import ResourceManager

@dataclass
class ConfigHelpSubcommand(Subcommand):
    name: str = 'config-help'
    short_help: str = "Show help for configuration options"
    long_help: str = "Display help information for all available configuration options."

    @staticmethod
    def func(args, **kwargs):
        config = Config().configure_new()
        directives = args.directives
        iprompt='pestifer-help: ' if args.interactive else ''
        config.console_help(directives=directives, interactive_prompt=iprompt, exit=True)

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('directives', type=str, nargs='*', default=[], help='specific directives to get help for (default: all)')
        self.parser.add_argument('--interactive', default=True, action='store_false', help='enter interactive mode to query config options (default: %(default)s)')

@dataclass
class ConfigDefaultSubcommand(Subcommand):
    name: str = 'config-default'
    short_help: str = "Show default configuration options"
    long_help: str = "Display default values for all available configuration options."

    @staticmethod
    def func(args, **kwargs):
        config = Config().configure_new()
        directives = args.directives
        specs=config.make_default_specs(*directives)
        print(yaml.dump(specs))

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('directives', type=str, nargs='*', default=[], help='specific directives to get default values for (default: all)')

@dataclass
class NewSystemSubcommand(Subcommand):
    name: str = 'new-system'
    short_help: str = 'create a new system script from scratch'
    long_help: str = 'Generate a new system script with a basic template.'

    @staticmethod
    def func(args, **kwargs):
        r = ResourceManager()
        build_type = 'full' if args.full else 'minimal'
        r.example_manager.new_example_yaml(id=args.id, build_type=build_type)

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('id', type=str, default=None, help='PDB/AlphaFold ID  of the new system to create')
        self.parser.add_argument('--full', default=False, action='store_true', help='Include full set of tasks in the new system configuration (default: %(default)s)')
        return self.parser

@dataclass
class ShowResourcesSubcommand(Subcommand):
    name: str = 'show-resources'
    short_help: str = 'show available resources for the current configuration'
    long_help: str = 'Display a list of all available resources for the current configuration.'

    @staticmethod
    def func(args, **kwargs):
        r = ResourceManager()
        resource_type = args.resource_type
        specs = {resource_type: []}
        if hasattr(args, resource_type):
            for subresource in getattr(args, resource_type):
                specs[resource_type].append(subresource)
        if args.user_pdbcollection:
            r.update_charmmff(user_pdbrepository_paths=args.user_pdbcollection)
        r.show(out_stream=print, components=specs, fullnames=args.fullnames, missing_fullnames=r.labels.residue_fullnames)

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('resource_type', type=str, default='examples', help='Type of resource to show; [tcl|examples|charmmff]')
        self.parser.add_argument('--charmmff', type=str, nargs='+', default=[], help='show sub-resources of charmmff resources (\'toppar\', \'custom\', \'pdb\')')
        self.parser.add_argument('--fullnames', default=False, action='store_true', help='show full names of any residues shown with --charmmff pdb')
        self.parser.add_argument('--user-pdbcollection', type=str, nargs='+', default=[], help='additional collections of PDB files outside pestifer installation')
        return self.parser

@dataclass
class WhereTCLSubcommand(Subcommand):
    name: str = 'where-tcl'
    short_help: str = "provides path of TcL scripts for sourcing in interactive VMD"
    long_help: str = "Display the path of TcL scripts for sourcing in interactive VMD."

    @staticmethod
    def func(args, **kwargs):
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

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('--startup-script-path', default=False, action='store_true', help='print full path of VMD startup script used by pestifer')
        self.parser.add_argument('--verbose', default=False, action='store_true', help='print full paths of all TcL resources of Pestifer')
        self.parser.add_argument('--script-dir', default=False, action='store_true', help='print full path of directory of Pestifer\'s VMD scripts')
        self.parser.add_argument('--pkg-dir', default=False, action='store_true', help='print full path of directory of Pestifer\'s VMD/TcL package library')
        self.parser.add_argument('--root', default=False, action='store_true', help='print full path of Pestifer\'s root TcL directory')
        return self.parser