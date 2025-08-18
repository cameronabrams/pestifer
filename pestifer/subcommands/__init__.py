from .subcommand import Subcommand
from .run import RunSubcommand
from .example_subcommands import FetchExampleSubcommand, RunExampleSubcommand
from .namd_subcommands import FollowNAMDLogSubcommand, MakeNAMDRestartSubcommand
from .config_subcommands import ConfigHelpSubcommand, ConfigDefaultSubcommand
from .desolvate_subcommand import DesolvateSubcommand
from .modify_package import ModifyPackageSubcommand
from .make_pdbcollection import MakePDBCollectionSubcommand
from .mdplot_subcommand import MDPlotSubcommand

_subcommands: list[Subcommand] = [
    RunSubcommand(),
    FetchExampleSubcommand(),
    RunExampleSubcommand(),
    FollowNAMDLogSubcommand(),
    ConfigHelpSubcommand(),
    ConfigDefaultSubcommand(),
    DesolvateSubcommand(),
    ModifyPackageSubcommand(),
    MakePDBCollectionSubcommand(),
    MDPlotSubcommand(),
    MakeNAMDRestartSubcommand(),
    ]