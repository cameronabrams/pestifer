# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The cache subcommand.  Inspect, clear, or rebuild pestifer's on-disk caches (the parsed
CHARMM force field, the PDB repository, and the residue-name lookup index).
"""
import argparse as ap

from dataclasses import dataclass
from datetime import datetime

from . import Subcommand

from ..util.cacheable_object import CacheableObject


def _human(nbytes: int) -> str:
    size = float(nbytes)
    for unit in ('B', 'KB', 'MB', 'GB'):
        if size < 1024 or unit == 'GB':
            return f'{size:.0f} {unit}' if unit == 'B' else f'{size:.1f} {unit}'
        size /= 1024


def _cache_status(out=print):
    d = CacheableObject.cache_directory()
    files = CacheableObject.cache_files()
    out(f'pestifer cache directory: {d}')
    if not files:
        out('  (empty -- no caches have been built yet)')
        return
    total = 0
    for f in files:
        st = f.stat()
        total += st.st_size
        # cache files are named "<prefix>-<classname>-<label>-v<ver>.joblib"
        parts = f.stem.split('-')
        kind = parts[1] if len(parts) > 1 else f.stem
        when = datetime.fromtimestamp(st.st_mtime).strftime('%Y-%m-%d %H:%M')
        out(f'  {kind:<26s} {_human(st.st_size):>9s}  {when}')
    out(f'  {len(files)} file(s), {_human(total)} total')


def _cache_clear(out=print):
    removed = CacheableObject.clear_cache()
    out(f'Removed {len(removed)} cache file(s) from {CacheableObject.cache_directory()}')


def _cache_rebuild(out=print):
    from ..core.resourcemanager import ResourceManager
    from ..charmmff.charmmffcontent import CHARMMFFContent, ResnameIndex
    rm = ResourceManager()
    for version_dir in rm.charmmff_version_dirs():
        out(f'Rebuilding caches for CHARMM force field "{version_dir.name}"...')
        CC = CHARMMFFContent(version_dir, force_rebuild=True)   # content metadata
        CC.provision(force_rebuild=True)                        # residue collection + PDB repository
        ResnameIndex(version_dir, force_rebuild=True)           # residue-name lookup index
    out('Done.')


@dataclass
class CacheSubcommand(Subcommand):
    name: str = 'cache'
    short_help: str = "inspect, clear, or rebuild pestifer's on-disk caches"
    long_help: str = ("Manage pestifer's per-user caches (the parsed CHARMM force field, the PDB "
                      "repository, and the residue-name lookup index): 'status' lists them, 'clear' "
                      "deletes them, and 'rebuild' force-rebuilds them.")

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        action = args.action
        if action == 'status':
            _cache_status()
        elif action == 'clear':
            _cache_clear()
        elif action == 'rebuild':
            _cache_rebuild()
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('action', type=str, choices=['status', 'clear', 'rebuild'],
                                 help="status: list the caches; clear: delete them; "
                                      "rebuild: force-rebuild them")
        return self.parser
