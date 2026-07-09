# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
On-demand generation cache for PDB-repository entries under ``~/.pestifer/``.

When a build references a residue that is defined in the CHARMM force field but has **no**
entry in the shipped PDB repository, pestifer can generate the needed coordinates on the fly
and cache them per-user, keyed by CHARMMFF release, so that subsequent builds reuse the result
instead of hard-erroring.  The cache lives at ``~/.pestifer/pdbrepository/<release>/<stream>/``
and mirrors the shipped repository's on-disk layout, so a cached collection is registered like
any other user PDB collection (auto-registered at provisioning; see
:meth:`CHARMMFFContent.provision_pdbrepository`).

The *kind* of artifact generated is driven by the **consumer**, not the species:

* the ``solvate`` task needs a pre-equilibrated periodic **box** (``kind: box``); a missing box
  is built with :func:`~pestifer.charmmff.make_solvent_box.make_solvent_box` and cached under
  ``.../<release>/solvent/<RESN>/`` (see :func:`ensure_solvent_box`);
* the grid membrane packer needs single-molecule **conformers** (``kind: molecule``); a missing
  lipid is built with :func:`~pestifer.charmmff.make_pdb_collection.do_resi` and cached under
  ``.../<release>/lipid/<RESN>/`` (see :func:`ensure_lipid_conformer`).

Generation is treated like a compile-cache: the first build is slow (a full pack/sample +
minimize + equilibration) and loudly logged; every build thereafter is a cache hit.  Builds are
concurrency-safe: they are serialized with an advisory ``fcntl`` lock (auto-released if a builder
dies, so no stale locks), and a finished entry is published into place with an atomic
``os.rename`` so a concurrent reader never sees a half-written entry.  Each build runs in a
**fresh** :class:`~pestifer.core.resourcemanager.ResourceManager` so the running build's
force-field state is not disturbed.
"""
import fcntl
import logging
import os
import shutil
from contextlib import contextmanager
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)

PDBCACHE_ROOT = Path('~/.pestifer/pdbrepository').expanduser()
"""Root of the per-user, release-keyed PDB-repository generation cache."""


def cache_release_root(release_key: str) -> Path:
    """Return the cache root for a CHARMMFF release (e.g. ``~/.pestifer/pdbrepository/feb26``)."""
    return PDBCACHE_ROOT / release_key


def cached_collection_dirs(release_key: str) -> list[Path]:
    """Return the collection directories present in the cache for ``release_key``.

    Used at provisioning time to auto-register previously generated collections as user PDB
    collections, so cached entries are available without any regeneration.
    """
    root = cache_release_root(release_key)
    if not root.is_dir():
        return []
    return sorted(p for p in root.iterdir() if p.is_dir() and not p.name.startswith('.'))


def _has_entry(collection_dir: Path, resname: str) -> bool:
    """True if ``collection_dir/<resname>/info.yaml`` exists (a complete cached entry)."""
    return (collection_dir / resname / 'info.yaml').is_file()


def _release_key(CC) -> str:
    """The resolved release directory name (e.g. ``feb26``) that keys the cache for ``CC``."""
    return os.path.basename(str(CC.charmmff_path))


@contextmanager
def _entry_lock(collection_dir: Path, resname: str):
    """Serialize builders of ``collection_dir/<resname>`` with an advisory ``fcntl`` lock.

    The lock is released automatically if the holding process dies (its fd closes), so there is
    no stale-lock hazard.  On contention the wait is announced, then blocks until acquired.
    """
    collection_dir.mkdir(parents=True, exist_ok=True)
    lock_path = collection_dir / f'.{resname}.lock'
    with open(lock_path, 'w') as lockf:
        try:
            fcntl.flock(lockf, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            logger.info(f'another process is building {resname}; waiting for it')
            fcntl.flock(lockf, fcntl.LOCK_EX)
        yield


def _fresh_charmmff(release_str: str):
    """Return an isolated ``(ResourceManager, CHARMMFFContent)`` for building a cache entry.

    Uses a fresh ResourceManager (seeded with the same release string) so generation does not
    mutate the running build's force-field configuration.
    """
    from ..core.resourcemanager import ResourceManager
    RM = ResourceManager(charmmff_config={'release': release_str} if release_str else {})
    CC = RM.charmmff_content
    CC.provision()
    return RM, CC


def _announce_build(resname: str, what: str, collection_dir: Path):
    logger.warning('=' * 70)
    logger.warning(f'pestifer: no PDB-repository entry for {resname!r}; building a {what} now.')
    logger.warning(f'This is a one-time cost.  The result is cached under {collection_dir} and '
                   f'reused by future builds.')
    logger.warning('=' * 70)


def _fresh_dir(path: Path) -> Path:
    shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True)
    return path


def ensure_solvent_box(resname: str, CC, *, nmol: int = 216, density: float = 1.0,
                       npt_steps: int = 50000, minimize_steps: int = 1000,
                       temperature: float = 300.0, pressure: float = 1.0,
                       seed: int = None, key_atom: str = None) -> Path:
    """Ensure a cached ``kind: box`` entry for solvent ``resname`` and return the collection dir.

    If the box is already cached this is a fast no-op.  Otherwise the box is built with
    :func:`make_solvent_box` in an isolated working directory and published atomically.

    Parameters
    ----------
    resname : str
        The solvent RESI name (must be defined in the CHARMM force field).
    CC : CHARMMFFContent
        The live force-field content; its ``charmmff_path`` keys the cache and its ``release_str``
        seeds the isolated builder.
    nmol, density, npt_steps, minimize_steps, temperature, pressure, seed, key_atom
        Forwarded to :func:`make_solvent_box`.  Defaults match the curated shipped boxes (full
        50k-step NPT equilibration); tests may shorten them.

    Returns
    -------
    Path
        The ``solvent`` collection directory, ready to register with ``PDBRepository.add_resource``.
    """
    from .make_solvent_box import make_solvent_box

    collection_dir = cache_release_root(_release_key(CC)) / 'solvent'
    if _has_entry(collection_dir, resname):
        logger.debug(f'auto-cache hit: {resname} solvent box at {collection_dir / resname}')
        return collection_dir

    with _entry_lock(collection_dir, resname):
        if _has_entry(collection_dir, resname):
            logger.info(f'{resname} solvent box was built by another process; using it')
            return collection_dir
        _announce_build(resname, 'solvent box (pack + minimize + NPT equilibration)', collection_dir)

        RM, build_CC = _fresh_charmmff(getattr(CC, 'release_str', ''))
        if resname not in build_CC:
            raise ValueError(f'RESI {resname} is not defined in the CHARMM force field; '
                             f'cannot auto-generate a solvent box')
        work_dir = _fresh_dir(collection_dir / f'.{resname}.work')
        staged_dir = collection_dir / f'.{resname}.staged'
        shutil.rmtree(staged_dir, ignore_errors=True)
        cwd = os.getcwd()
        os.chdir(work_dir)
        try:
            info = make_solvent_box(resname, build_CC, RM=RM, nmol=nmol, density=density,
                                    temperature=temperature, pressure=pressure,
                                    minimize_steps=minimize_steps, npt_steps=npt_steps,
                                    seed=seed, key_atom=key_atom)
            info['quality'] = 'auto'   # provenance: auto-generated, not hand-curated
            with open('info.yaml', 'w') as f:
                yaml.dump(info, f)
            staged_dir.mkdir()
            for fname in ('info.yaml', info['psf'], info['pdb']):
                shutil.copyfile(fname, staged_dir / fname)
        finally:
            os.chdir(cwd)
        os.rename(staged_dir, collection_dir / resname)   # atomic publish
        shutil.rmtree(work_dir, ignore_errors=True)
        logger.warning(f'pestifer: cached {resname} solvent box '
                       f"(edge {info['box_edge']} A, density {info['density']} g/cc) at "
                       f'{collection_dir / resname}')
        return collection_dir


def ensure_lipid_conformer(resname: str, CC, *, nsamples: int = 10, sample_steps: int = 5000,
                           minimize_steps: int = 500, sample_temperature: float = 300.0,
                           force_constant: float = 1.0) -> Path:
    """Ensure a cached ``kind: molecule`` conformer entry for ``resname`` and return the
    (``lipid``) collection dir.

    If the conformers are already cached this is a fast no-op.  Otherwise they are built with
    :func:`do_resi` in an isolated working directory and published atomically.  ``do_resi`` writes
    the finished entry directory (``<out>/<resname>/`` with ``info.yaml`` + conformer PDBs) itself;
    we mark it ``quality: auto`` and move it into the cache.

    Parameters
    ----------
    resname : str
        The RESI name (must be defined in the CHARMM force field).
    CC : CHARMMFFContent
        The live force-field content (see :func:`ensure_solvent_box`).
    nsamples, sample_steps, minimize_steps, sample_temperature, force_constant
        Forwarded to :func:`do_resi`.  Defaults match ``make-pdb-collection``; tests may shorten.

    Returns
    -------
    Path
        The ``lipid`` collection directory, ready to register with ``PDBRepository.add_resource``.
    """
    from .make_pdb_collection import do_resi

    collection_dir = cache_release_root(_release_key(CC)) / 'lipid'
    if _has_entry(collection_dir, resname):
        logger.debug(f'auto-cache hit: {resname} conformer at {collection_dir / resname}')
        return collection_dir

    with _entry_lock(collection_dir, resname):
        if _has_entry(collection_dir, resname):
            logger.info(f'{resname} conformer was built by another process; using it')
            return collection_dir
        _announce_build(resname, 'single-molecule conformer set (sample + minimize)', collection_dir)

        RM, build_CC = _fresh_charmmff(getattr(CC, 'release_str', ''))
        if resname not in build_CC:
            raise ValueError(f'RESI {resname} is not defined in the CHARMM force field; '
                             f'cannot auto-generate a conformer')
        work_dir = _fresh_dir(collection_dir / f'.{resname}.work')
        cwd = os.getcwd()
        os.chdir(work_dir)
        try:
            do_resi(resname, build_CC, RM=RM, outdir='out', faildir='fail', cleanup=True,
                    nsamples=nsamples, sample_steps=sample_steps, minimize_steps=minimize_steps,
                    sample_temperature=sample_temperature, force_constant=force_constant)
            produced = Path('out') / resname
            if not (produced / 'info.yaml').is_file():
                raise RuntimeError(f'conformer generation for {resname} failed (no entry produced; '
                                   f'see {work_dir / "fail" / resname} if present)')
            info = yaml.safe_load((produced / 'info.yaml').read_text())
            info['quality'] = 'auto'
            (produced / 'info.yaml').write_text(yaml.dump(info))
        finally:
            os.chdir(cwd)
        os.rename(work_dir / 'out' / resname, collection_dir / resname)   # atomic publish
        shutil.rmtree(work_dir, ignore_errors=True)
        logger.warning(f'pestifer: cached {resname} conformer set '
                       f"({len(info.get('conformers', []))} conformers) at {collection_dir / resname}")
        return collection_dir
