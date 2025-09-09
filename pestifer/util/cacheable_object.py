
# Author: ChatGPT 5
# Assistant Author: Cameron F. Abrams <cfa22@drexel.edu>

"""
Implements the general-purpose CacheableObject class and the TarBytesFS class.
"""

from __future__ import annotations

import fsspec
import hashlib
import io
import joblib
import logging
import os
import tempfile

from filelock import FileLock
from packaging.version import Version, InvalidVersion
from pathlib import Path
from platformdirs import user_cache_dir
from typing import Iterable

logger = logging.getLogger(__name__)

from .stringthings import __pestifer_version__

class TarBytesFS:
    """
    A class for handling tarred bytes as a filesystem.
    """
    __slots__ = ("tar_bytes", "compression", "_fs")

    def __init__(self, tar_bytes: bytes, compression: str | None = None):
        self.tar_bytes = tar_bytes          # pickleable
        self.compression = compression      # pickleable
        self._fs = None                     # ephemeral, NOT pickleable

    @classmethod
    def from_file(cls, path: str | Path, compression: str | None = None):
        """
        Initializes a TarBytesFS instance from a file.
        """
        return cls(Path(path).read_bytes(), compression=compression)

    def fs(self):
        """
        Returns a FileSystem handle for the TarBytesFS object 
        """
        if self._fs is None:
            bio = io.BytesIO(self.tar_bytes)
            if self.compression is None and not hasattr(bio, "name"):
                # Let fsspec infer if you wish: bio.name = "archive.tgz"
                pass
            self._fs = fsspec.filesystem("tar", fo=bio, compression=self.compression)
        return self._fs

    def ls(self, path: str = ""):
        """
        Lists the contents of a directory in the tarred filesystem.
        """
        return self.fs().ls(path)

    def open(self, path: str, mode: str = "rb"):
        """
        Opens a file in the tarred filesystem and returns a context handle to that file.
        """
        return self.fs().open(path, mode)

    # --- Pickle hooks: drop the ephemeral filesystem
    def __getstate__(self):
        return {"tar_bytes": self.tar_bytes, "compression": self.compression}

    def __setstate__(self, state):
        self.tar_bytes = state["tar_bytes"]
        self.compression = state.get("compression")
        self._fs = None

def _latest_mtime(root: Path,
                  *,
                  ignore_names: set[str] = {"__pycache__"},
                  ignore_suffixes: set[str] = set()) -> float:
    latest = 0.0
    if isinstance(root, CacheableObject):
        root = root.resource_root
    root = Path(root)
    for p in root.rglob("*"):
        if p.name.startswith(".") or p.name in ignore_names:
            continue
        if p.is_file() and p.suffix not in ignore_suffixes:
            try:
                mt = p.stat().st_mtime
                if mt > latest:
                    latest = mt
            except FileNotFoundError:
                pass
    return latest

def _hash_resource(resource: Path) -> str:
    return hashlib.sha256(str(Path(resource).resolve()).encode()).hexdigest()[:12]

class CacheableObject:
    """
    Base/mixin providing native cache/de-cache behavior controlled by resource mtimes.

    Subclasses must implement: ``_build_from_resources(self, resource_root: Path) -> None``
    to populate `self` from files in resource_root.

    Behavior:

      - On __init__, compare cache mtime to newest file mtime under ``resource_root``.
      - If cache is fresh, hydrate ``self`` from cache and set ``self.from_cache = True``.
      - Else, call ``_build_from_resources(...)``, ensure ``self.from_cache = False``, and write cache.

    """

    APP_NAME = "pestifer"          # for per-user cache dir
    APP_VERSION = __pestifer_version__
    CACHE_PREFIX = "cacheobj"      # file prefix; class & path hash appended
    CACHE_COMPRESS = ("gzip", 3)   # joblib compression

    # Choose how versioning gates the cache:
    #   "major"        -> v{MAJOR}
    #   "major.minor"  -> v{MAJOR}.{MINOR}   (recommended default)
    #   "exact"        -> full version string
    #   None           -> no version in filename (old behavior)
    VERSION_SCOPE = "major.minor"

    @classmethod
    def _version_tag(cls) -> str:
        if not cls.VERSION_SCOPE:
            return ""
        try:
            v = Version(cls.APP_VERSION)
        except InvalidVersion:
            # fall back to raw string if non-PEP440
            return f"v{cls.APP_VERSION}"
        if cls.VERSION_SCOPE == "major":
            return f"v{v.major}"
        if cls.VERSION_SCOPE == "major.minor":
            return f"v{v.major}.{v.minor}"
        if cls.VERSION_SCOPE == "exact":
            return f"v{str(v)}"
        return ""  # unknown scope => no tag

    def __init__(self,
                 resource_root: str | Path,
                 *,
                 cache_dir: str | Path | object | None = None,
                 force_rebuild: bool = False,
                 ignore_suffixes: Iterable[str] = (),
                 **kwargs):
        if isinstance(resource_root, str):
            resource_root = Path(resource_root)
        self.resource_root = resource_root
        base_key = f"{self.__class__.__name__.lower()}-{_hash_resource(resource_root)}"
        vtag = self._version_tag()
        key = f"{base_key}-{vtag}" if vtag else base_key
        cdir = Path(cache_dir) if cache_dir else Path(user_cache_dir(self.APP_NAME))
        cdir.mkdir(parents=True, exist_ok=True)
        cpath = cdir / f"{self.CACHE_PREFIX}-{key}.joblib"
        lock = FileLock(str(cpath) + ".lock")

        resources_mtime = _latest_mtime(resource_root, ignore_suffixes=set(ignore_suffixes))
        # logger.debug("Resource mtime: %s", resources_mtime)
        # logger.debug("Cache mtime: %s", cpath.stat().st_mtime if cpath.exists() else None)
        with lock:
            # logger.debug(f'Acquired lock for {cpath}')
            # logger.debug(f'cache exists: {cpath.exists()} x force_rebuild: {force_rebuild} => load cache? {"yes" if cpath.exists() and not force_rebuild else "no"}')
            if cpath.exists() and not force_rebuild:
                # logger.debug(f"{cpath.stat().st_mtime} >=? {resources_mtime} -> {cpath.stat().st_mtime >= resources_mtime}")
                try:
                    if cpath.stat().st_mtime >= resources_mtime:
                        # logger.info(f'Loading {self.__class__.__name__} from cache: {cpath}')
                        cached = joblib.load(cpath)  # trusted cache only
                        # logger.debug("Loaded cached object")
                        self._adopt_state_from(cached)
                        # mark as loaded-from-cache (even if subclass set it earlier)
                        try:
                            object.__setattr__(self, "from_cache", True)
                        except Exception:
                            # if __slots__ disallow it, ignore
                            pass
                        # logger.debug("Loaded from cache: %s", cpath)
                        return
                except Exception:
                    # fall through to rebuild on any load problem
                    # logger.debug("Cache load failed; rebuilding...")
                    pass

            # Rebuild from resources
            logger.info(f'Rebuilding {self.__class__.__name__} from resources...')
            self._build_from_resources(resource_root, **kwargs)
            # if subclass didn't set it, default to False
            if not hasattr(self, "from_cache"):
                try:
                    object.__setattr__(self, "from_cache", False)
                except Exception:
                    pass

            # Persist atomically
            fd, tmp = tempfile.mkstemp(dir=str(cdir), suffix=".joblib")
            os.close(fd)
            try:
                logger.debug(f'Writing {self.__class__.__name__} to cache: {cpath}')
                joblib.dump(self, tmp, compress=self.CACHE_COMPRESS)
                os.replace(tmp, cpath)
                # logger.debug("Wrote cache: %s", cpath)
            finally:
                if os.path.exists(tmp):
                    os.remove(tmp)

    # ----- hooks -----
    def _build_from_resources(self, resource_root: Path, **kwargs) -> None:
        """Subclasses must populate `self` here."""
        raise NotImplementedError

    def _adopt_state_from(self, other: "CacheableObject") -> None:
        """Default hydration copies __dict__; override if using __slots__."""
        if hasattr(self, "__dict__") and hasattr(other, "__dict__"):
            self.__dict__.clear()
            self.__dict__.update(other.__dict__)
        else:
            # basic slots copy (doesn't walk MRO deeply; override for complex cases)
            for slot in getattr(type(self), "__slots__", ()):
                setattr(self, slot, getattr(other, slot))
