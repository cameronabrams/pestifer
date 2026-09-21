# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The cache subcommand.  Inspect, clear, rebuild, or pre-populate pestifer's on-disk caches (the
parsed CHARMM force field, the PDB repository, and the residue-name lookup index).
"""
import argparse as ap
import logging

from dataclasses import dataclass
from datetime import datetime

from . import Subcommand

from ..util.cacheable_object import CacheableObject

logger = logging.getLogger(__name__)


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


def _cache_prebuild(args):
    """Generate one lipid conformer set into the per-user cache, ahead of any build.

    A build that meets an uncached ``(lipid, phase)`` pair generates the set itself, which is
    usually what you want.  This exists for the case where that is inconvenient or impossible:
    generating inside a batch job spends allocation time on a single-molecule vacuum run, and on a
    cluster whose MPI cannot be direct-launched it used to fail outright (fixed, but a pre-populated
    cache still saves the wait).  Run it on a login node, once per lipid and phase; every later
    build on that filesystem is a cache hit.

    The sampler is not a choice here: it is the one a build would use for that phase, so the cached
    entry is the entry the build would have made.
    """
    from ..charmmff.autocache import ensure_lipid_conformer, phase_entry_name
    from ..core.resourcemanager import ResourceManager

    resname = (getattr(args, 'resname', '') or '').strip()
    if not resname:
        raise ValueError("cache prebuild needs --resname (e.g. 'pestifer cache prebuild "
                         "--resname PSM --phase Lo')")
    phase = getattr(args, 'phase', 'Ld')
    release = getattr(args, 'charmmff_release', '')
    RM = ResourceManager(charmmff_config={'release': release} if release else {})
    CC = RM.charmmff_content
    CC.provision()
    if resname not in CC:
        raise ValueError(f'RESI {resname} is not defined in the CHARMM force field')
    # 'mc' for a phased ensemble, matching Bilayer's own choice: a prebuilt entry that a build
    # would not have made is worse than no entry at all
    sampler = 'mc' if phase in ('Ld', 'Lo') else 'md'
    entry = phase_entry_name(resname, phase)
    logger.info(f'prebuilding conformer set {entry} ({sampler} sampler)')
    collection_dir = ensure_lipid_conformer(resname, CC, phase=phase, sampler=sampler)
    logger.info(f'{entry} is cached at {collection_dir / entry}')
    print(collection_dir / entry)


@dataclass
class CacheSubcommand(Subcommand):
    name: str = 'cache'
    group: str = 'Manage the installation'
    short_help: str = "inspect, clear, rebuild, or pre-populate pestifer's on-disk caches"
    long_help: str = ("Manage pestifer's per-user caches (the parsed CHARMM force field, the PDB "
                      "repository, and the residue-name lookup index): 'status' lists them, 'clear' "
                      "deletes them, 'rebuild' force-rebuilds them, and 'prebuild' generates one "
                      "lipid conformer set ahead of the build that would need it.")

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        action = args.action
        if action == 'status':
            _cache_status()
        elif action == 'clear':
            _cache_clear()
        elif action == 'rebuild':
            _cache_rebuild()
        elif action == 'prebuild':
            _cache_prebuild(args)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('action', type=str,
                                 choices=['status', 'clear', 'rebuild', 'prebuild'],
                                 help="status: list the caches; clear: delete them; "
                                      "rebuild: force-rebuild them; prebuild: generate one lipid "
                                      "conformer set into the cache now")
        self.parser.add_argument('--resname', type=str, default='',
                                 help="prebuild: the lipid RESI to generate (e.g. PSM)")
        self.parser.add_argument('--phase', type=str, default='Ld', choices=['Ld', 'Lo'],
                                 help="prebuild: the bilayer phase the conformer ensemble is tuned "
                                      "for; 'Lo' caches as <RESI>__Lo (default: %(default)s)")
        self.parser.add_argument('--charmmff-release', dest='charmmff_release', type=str, default='',
                                 help='prebuild: CHARMMFF release to build against (default: the '
                                      'newest available, which is what a build uses)')
        return self.parser
