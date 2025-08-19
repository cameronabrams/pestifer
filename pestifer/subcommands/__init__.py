import os
from pathlib import Path
from .subcommand import Subcommand
from .run import RunSubcommand
from .example_subcommands import FetchExampleSubcommand, RunExampleSubcommand
from .namd_subcommands import FollowNAMDLogSubcommand, MakeNAMDRestartSubcommand
from .config_subcommands import ConfigHelpSubcommand, ConfigDefaultSubcommand, ShowResourcesSubcommand, NewSystemSubcommand, WhereTCLSubcommand
from .desolvate_subcommand import DesolvateSubcommand
from .modify_package import ModifyPackageSubcommand
from .make_pdbcollection import MakePDBCollectionSubcommand
from .mdplot_subcommand import MDPlotSubcommand

package_path = Path(__file__).resolve().parent.parent.parent
is_source_package_with_git = os.path.isdir(os.path.join(package_path, '.git'))

_subcommands: list[Subcommand] = [
    RunSubcommand(),
    FetchExampleSubcommand(),
    RunExampleSubcommand(),
    ConfigHelpSubcommand(),
    ConfigDefaultSubcommand(),
    ShowResourcesSubcommand(),
    NewSystemSubcommand(),
    DesolvateSubcommand(),
    MakePDBCollectionSubcommand(),
    WhereTCLSubcommand(),
    MDPlotSubcommand(),
    MakeNAMDRestartSubcommand(),
    FollowNAMDLogSubcommand(),
    ]

if is_source_package_with_git:
    _subcommands.append(ModifyPackageSubcommand())