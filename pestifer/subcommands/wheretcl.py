# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The wheretcl subcommand.  Exposes true locations of various TcL resources in the pestifer
installation, useful for custom VMD scripts.
"""
import os
import argparse as ap

from dataclasses import dataclass

from . import Subcommand

from ..core.resourcemanager import ResourceManager

@dataclass
class WhereTCLSubcommand(Subcommand):
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
