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
from .report_methods import ReportMethodsSubcommand
from .wheretcl import WhereTCLSubcommand
from .desolvate import DesolvateSubcommand
from .modify_package import ModifyPackageSubcommand
from .make_ligand_mol2 import MakeLigandMol2Subcommand
from .make_pdbcollection import MakePDBCollectionSubcommand
from .cache import CacheSubcommand
from .mdplot import MDPlotSubcommand
from .density_profile import DensityProfileSubcommand
from .pressure_profile_ewald import PressureProfileEwaldSubcommand
from .setup_vmd import SetupVMDSubcommand
from .setup_claude import SetupClaudeSubcommand

GROUPS: tuple[str, ...] = (
    'Build a system',
    'Configure a build',
    'Work with structures',
    'After the run',
    'Manage the installation',
)
"""The order subcommand groups are presented in, on the command line and in the docs.

Commands are grouped by *what you hand them* -- a config, a structure, a NAMD run's output, or
nothing at all -- rather than by what they do, so a reader can place an unfamiliar command without
having read the whole list.  The order follows a build's life: describe it, check it, build it,
look at what came out.

Note what the 'After the run' group contains.  Four of its five (``mdplot``, ``density-profile``,
``follow-namd-log``, ``make-namd-restart``) read NAMD output and touch nothing of pestifer's, so
they work equally well on a run pestifer never managed; only ``report-methods`` needs the
``run-record.json`` a pestifer build leaves behind.  That is the one real seam in the command set,
and naming it is the point -- see ``docs/design/namd-layer-extraction.md``.
"""

package_path = Path(__file__).resolve().parent.parent.parent
is_source_package_with_git = os.path.isdir(os.path.join(package_path, '.git'))

# Declared in presentation order within each group; the stable sort below only orders the groups
# themselves, so this list decides what comes first inside one.
_subcommands: list[Subcommand] = [
    # Build a system
    RunSubcommand(),
    RunExampleSubcommand(),
    FetchExampleSubcommand(),
    NewSystemSubcommand(),
    # Configure a build
    ConfigHelpSubcommand(),
    ConfigDefaultSubcommand(),
    ShowResourcesSubcommand(),
    # Work with structures
    DesolvateSubcommand(),
    MakeLigandMol2Subcommand(),
    MakePDBCollectionSubcommand(),
    # After the run
    MDPlotSubcommand(),
    DensityProfileSubcommand(),
    PressureProfileEwaldSubcommand(),
    FollowNAMDLogSubcommand(),
    MakeNAMDRestartSubcommand(),
    ReportMethodsSubcommand(),
    # Manage the installation
    CacheSubcommand(),
    WhereTCLSubcommand(),
    SetupVMDSubcommand(),
    SetupClaudeSubcommand(),
    ]

if is_source_package_with_git:
    # a maintainer tool; it does not exist for a pip-installed pestifer
    _subcommands.append(ModifyPackageSubcommand())

# Present in group order.  Sorting here rather than by hand keeps the declaration list above free
# to stay in import order, and makes a subcommand added without a group loud rather than silent:
# it sorts to the end instead of landing somewhere arbitrary in the middle.
_subcommands.sort(key=lambda s: (GROUPS.index(s.group) if s.group in GROUPS else len(GROUPS),))


def grouped_subcommands() -> list[tuple[str, list[Subcommand]]]:
    """``[(group title, [subcommand, ...]), ...]`` in presentation order."""
    out = []
    for g in GROUPS:
        members = [s for s in _subcommands if s.group == g]
        if members:
            out.append((g, members))
    ungrouped = [s for s in _subcommands if s.group not in GROUPS]
    if ungrouped:
        out.append(('Other', ungrouped))
    return out