# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
On-demand generation cache for PDB-repository entries under ``~/.pestifer/``.

When a build references a residue that is defined in the CHARMM force field but has **no**
entry in the shipped PDB repository, pestifer can generate the needed coordinates on the fly
and cache them per-user, keyed by CHARMMFF release, so that subsequent builds reuse the result
instead of hard-erroring.  The cache lives at ``~/.pestifer/pdbrepository/<release>/<stream>/``
and mirrors the shipped repository's on-disk layout, so a cached collection is registered like
any other user PDB collection.

This module implements the **solvent-box** slice: a missing ``kind: box`` solvent entry is built
with :func:`~pestifer.charmmff.make_solvent_box.make_solvent_box` and cached under
``~/.pestifer/pdbrepository/<release>/solvent/<RESN>/``.  Generation is treated like a
compile-cache: the first build is slow (a full pack + minimize + NPT equilibration) and loudly
logged; every build thereafter is a cache hit.

Concurrency: builds are serialized with an advisory ``fcntl`` lock (auto-released if a builder
dies, so no stale locks), and a finished entry is published into place with an atomic
``os.rename`` so a concurrent reader never sees a half-written entry.
"""
import fcntl
import logging
import os
import shutil
from pathlib import Path

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


def ensure_solvent_box(resname: str, release_key: str, release_str: str = '', *,
                       nmol: int = 216, density: float = 1.0, npt_steps: int = 50000,
                       minimize_steps: int = 1000, temperature: float = 300.0,
                       pressure: float = 1.0, seed: int = None, key_atom: str = None) -> Path:
    """Ensure a cached ``kind: box`` entry for solvent ``resname`` and return the collection dir.

    If the box is already cached, this is a fast no-op returning the ``solvent`` collection
    directory.  Otherwise the box is built with
    :func:`~pestifer.charmmff.make_solvent_box.make_solvent_box` in an isolated working directory
    (using a **fresh** :class:`~pestifer.core.resourcemanager.ResourceManager` seeded with
    ``release_str`` so the running build's force-field state is not disturbed) and published into
    the cache atomically.

    Parameters
    ----------
    resname : str
        The solvent RESI name (must be defined in the CHARMM force field for ``release_str``).
    release_key : str
        The resolved release directory name (e.g. ``feb26``); keys the cache directory.
    release_str : str
        The release string as configured (e.g. ``February2026``, or ``''`` for the newest);
        used to construct the isolated ResourceManager.
    nmol, density, npt_steps, minimize_steps, temperature, pressure, seed, key_atom
        Forwarded to :func:`make_solvent_box`.  Defaults match the curated shipped boxes
        (full 50k-step NPT equilibration); tests may shorten them.

    Returns
    -------
    Path
        The ``solvent`` collection directory (``~/.pestifer/pdbrepository/<release>/solvent``),
        ready to be registered with ``PDBRepository.add_resource``.
    """
    import yaml

    from .make_solvent_box import make_solvent_box
    from ..core.resourcemanager import ResourceManager

    collection_dir = cache_release_root(release_key) / 'solvent'
    if _has_entry(collection_dir, resname):
        logger.debug(f'auto-cache hit: {resname} solvent box at {collection_dir / resname}')
        return collection_dir

    collection_dir.mkdir(parents=True, exist_ok=True)
    lock_path = collection_dir / f'.{resname}.lock'
    with open(lock_path, 'w') as lockf:
        try:
            fcntl.flock(lockf, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            logger.info(f'another process is building the {resname} solvent box; waiting for it')
            fcntl.flock(lockf, fcntl.LOCK_EX)
        # re-check under the lock: a concurrent builder may have finished while we waited
        if _has_entry(collection_dir, resname):
            logger.info(f'{resname} solvent box was built by another process; using it')
            return collection_dir

        logger.warning('=' * 70)
        logger.warning(f'pestifer: no PDB-repository box for solvent {resname!r}; building one now.')
        logger.warning(f'This is a one-time cost (pack + minimize + NPT equilibration).  The result '
                       f'is cached under {collection_dir} and reused by future builds.')
        logger.warning('=' * 70)

        work_dir = collection_dir / f'.{resname}.work'
        staged_dir = collection_dir / f'.{resname}.staged'
        for d in (work_dir, staged_dir):
            shutil.rmtree(d, ignore_errors=True)
        work_dir.mkdir()

        RM = ResourceManager(charmmff_config={'release': release_str} if release_str else {})
        CC = RM.charmmff_content
        CC.provision()
        if resname not in CC:
            shutil.rmtree(work_dir, ignore_errors=True)
            raise ValueError(f'RESI {resname} is not defined in the CHARMM force field for release '
                             f'{release_key!r}; cannot auto-generate a solvent box')

        cwd = os.getcwd()
        os.chdir(work_dir)
        try:
            info = make_solvent_box(resname, CC, RM=RM, nmol=nmol, density=density,
                                    temperature=temperature, pressure=pressure,
                                    minimize_steps=minimize_steps, npt_steps=npt_steps,
                                    seed=seed, key_atom=key_atom)
            # provenance: mark this entry as auto-generated (vs a hand-curated shipped box)
            info['quality'] = 'auto'
            with open('info.yaml', 'w') as f:
                yaml.dump(info, f)
            staged_dir.mkdir()
            for fname in ('info.yaml', info['psf'], info['pdb']):
                shutil.copyfile(fname, staged_dir / fname)
        finally:
            os.chdir(cwd)

        # atomic publish: a reader either sees no entry or a complete one
        os.rename(staged_dir, collection_dir / resname)
        shutil.rmtree(work_dir, ignore_errors=True)
        logger.warning(f'pestifer: cached {resname} solvent box '
                       f"(edge {info['box_edge']} A, density {info['density']} g/cc) at "
                       f'{collection_dir / resname}')
        return collection_dir
