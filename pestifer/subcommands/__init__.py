import os
from pathlib import Path
from ..cli.subcommand import Subcommand
from .build import RunSubcommand
from .fetch_example import FetchExampleSubcommand
from .build_example import RunExampleSubcommand
from .follow_namd_log import FollowNAMDLogSubcommand
from .make_namd_restart import MakeNAMDRestartSubcommand
from .config_help import ConfigHelpSubcommand
from .config_default import ConfigDefaultSubcommand
from .new_system import NewSystemSubcommand
from .show_resources import ShowResourcesSubcommand
from .wheretcl import WhereTCLSubcommand
from .desolvate import DesolvateSubcommand
from .modify_package import ModifyPackageSubcommand
from .make_pdbcollection import MakePDBCollectionSubcommand
from .rebuild_charmmff_cache import RebuildCHARMMFFCache
from .mdplot import MDPlotSubcommand
from .setup_vmd import SetupVMDSubcommand

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
    RebuildCHARMMFFCache(),
    FollowNAMDLogSubcommand(),
    SetupVMDSubcommand(),
    ]

if is_source_package_with_git:
    _subcommands.append(ModifyPackageSubcommand())