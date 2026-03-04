# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The rebuild-charmmff-cache subcommand.  Rebuilds the CHARMM force field cache from the
current topology and parameter files.
"""
import argparse as ap

from dataclasses import dataclass

from ..cli.subcommand import Subcommand

from ..charmmff.charmmffcontent import CHARMMFFContent

@dataclass
class RebuildCHARMMFFCache(Subcommand):
    name: str = 'rebuild-charmmff-cache'
    short_help: str = "rebuild the CHARMM force field cache"
    long_help: str = "Rebuild the CHARMM force field cache from the current topology and parameter files."

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        from ..core.resourcemanager import ResourceManager
        rm = ResourceManager()
        for version_dir in rm.charmmff_version_dirs():
            CC = CHARMMFFContent(version_dir, force_rebuild=True)
            CC.provision(force_rebuild=True)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        return self.parser
