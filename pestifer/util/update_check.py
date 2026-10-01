# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Tells an interactive user when the pestifer they are running is older than the latest release.

The whole point of this module is that it is *optional*: a build must behave identically whether
the check succeeds, fails, or never runs.  Three rules follow from that, and each is load-bearing.

**It only runs on a terminal.**  :func:`emit_update_notice` returns immediately unless the stream
it would write to is a tty.  A redirected run -- ``pestifer build x.yaml > run.log 2>&1``, which
is how the example sweep and most cluster jobs invoke it -- therefore makes no network call at
all and writes no notice.  That keeps build logs byte-identical across runs, which matters here:
the sweep compares logs between builds, and a line that appears only when the network happens to
be up is a diff that means nothing.

**A failed check is cached like a successful one.**  The interval stamp is written whether or not
the fetch worked, so a machine with no route out pays one timeout a day rather than one per
invocation.  Without that, an unreachable PyPI would be the *expensive* case rather than the
cheap one -- precisely backwards.

**A source tree never nags.**  :data:`~pestifer.util.stringthings.__pestifer_version_from_source__`
is true when the running version came from a working tree's ``pyproject.toml``, where being ahead
of the latest release is normal.

Opt out with ``--no-update-check``, ``PESTIFER_NO_UPDATE_CHECK=1``, or by setting ``"enabled":
false`` in the cache file (which is what ``pestifer check-update --disable`` writes, and the only
one of the three that persists for a user who cannot edit the command line).
"""
import json
import logging
import os
import time

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from packaging.version import InvalidVersion, Version

from .stringthings import __pestifer_version__, __pestifer_version_from_source__

logger = logging.getLogger(__name__)

PYPI_JSON_URL = 'https://pypi.org/pypi/pestifer/json'
"""PyPI's metadata endpoint for pestifer.  ``info.version`` is the latest non-yanked release."""

CHANGELOG_URL = 'https://github.com/cameronabrams/pestifer/blob/main/CHANGELOG.md'
"""Linked from the notice: a bare version number does not tell a user whether to care."""

CACHE_PATH = Path('~/.pestifer/update-check.json').expanduser()

CHECK_INTERVAL_SECONDS = 24 * 60 * 60

FETCH_TIMEOUT_SECONDS = 2.0
"""Deliberately short.  This runs before the user's command does, so it is latency they feel."""

ENV_OPT_OUT = 'PESTIFER_NO_UPDATE_CHECK'


@dataclass
class UpdateCheckResult:
    """What :func:`check_for_update` found."""

    current: str
    """The running version."""
    latest: str | None
    """The latest release on PyPI, or ``None`` if it could not be determined."""
    origin: str
    """Where :attr:`latest` came from: ``'pypi'``, ``'cache'``, ``'unavailable'``,
    ``'disabled'``, or ``'source-tree'``."""
    update_available: bool = False

    @property
    def notice(self) -> str | None:
        """The message to show the user, or ``None`` when there is nothing to say."""
        if not self.update_available:
            return None
        return (f'A newer pestifer is available: {self.latest} (you are running {self.current})\n'
                f'  upgrade:   pip install -U pestifer\n'
                f'  changelog: {CHANGELOG_URL}')


def _fetch_latest_from_pypi(timeout: float = FETCH_TIMEOUT_SECONDS) -> str:
    """Return the latest released version string from PyPI.  Raises on any failure."""
    import requests  # deferred: keeps `import pestifer` off requests' import cost
    response = requests.get(PYPI_JSON_URL, timeout=timeout,
                            headers={'Accept': 'application/json'})
    response.raise_for_status()
    return str(response.json()['info']['version'])


def _read_cache(cache_path: Path) -> dict:
    try:
        with open(cache_path) as f:
            data = json.load(f)
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def _write_cache(cache_path: Path, data: dict) -> None:
    """Replace the cache file atomically, or give up quietly.

    Atomically because a sweep can have several pestifer processes starting at once, and a
    half-written JSON file read by the next one is a parse error on a path whose whole job is
    to be harmless.  Quietly because a read-only or full home directory must not be fatal.
    """
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = cache_path.with_suffix(f'.json.{os.getpid()}.tmp')
        with open(tmp, 'w') as f:
            json.dump(data, f, indent=2)
        os.replace(tmp, cache_path)
    except OSError as e:
        logger.debug(f'update check: could not write {cache_path}: {e}')


def _is_newer(latest: str | None, current: str) -> bool:
    """``latest > current`` under PEP 440, or ``False`` if either is unparseable.

    Strictly greater, never merely different: a working tree between a version bump and its tag
    is legitimately ahead of PyPI, and so is anything installed from a git checkout.
    """
    if not latest:
        return False
    try:
        return Version(latest) > Version(current)
    except InvalidVersion:
        return False


def is_disabled(cache_path: Path = CACHE_PATH) -> bool:
    """Whether the user has opted out, by environment variable or persistently in the cache."""
    if os.environ.get(ENV_OPT_OUT, '').strip().lower() not in ('', '0', 'false', 'no'):
        return True
    return _read_cache(cache_path).get('enabled') is False


def set_enabled(enabled: bool, cache_path: Path = CACHE_PATH) -> None:
    """Persist the opt-out decision, for a user who cannot pass a flag or set an env var."""
    data = _read_cache(cache_path)
    data['enabled'] = bool(enabled)
    _write_cache(cache_path, data)


def check_for_update(current_version: str = __pestifer_version__,
                     *,
                     from_source: bool = __pestifer_version_from_source__,
                     force: bool = False,
                     now: float | None = None,
                     cache_path: Path = CACHE_PATH,
                     fetcher: Callable[[], str] = _fetch_latest_from_pypi) -> UpdateCheckResult:
    """Compare ``current_version`` against the latest release on PyPI.

    Consults the cache first and only contacts PyPI when the cached answer is older than
    :data:`CHECK_INTERVAL_SECONDS`.  ``force`` skips both the cache interval and the
    source-tree suppression, which is what ``pestifer check-update`` wants: an explicit
    question deserves a fresh answer.

    A network failure is not an error here -- it yields ``origin='unavailable'`` and stamps the
    cache anyway, so the next invocation does not retry immediately.  Programming errors are
    *not* caught; :func:`emit_update_notice` is the guard for the path where nothing may raise.
    """
    now = time.time() if now is None else now
    if not force:
        if is_disabled(cache_path):
            return UpdateCheckResult(current=current_version, latest=None, origin='disabled')
        if from_source:
            return UpdateCheckResult(current=current_version, latest=None, origin='source-tree')

    cache = _read_cache(cache_path)
    cached_latest = cache.get('latest')
    last_check = cache.get('last_check')
    fresh = (isinstance(last_check, (int, float)) and 0 <= now - last_check < CHECK_INTERVAL_SECONDS)
    if fresh and not force:
        return UpdateCheckResult(current=current_version, latest=cached_latest,
                                 origin='cache' if cached_latest else 'unavailable',
                                 update_available=_is_newer(cached_latest, current_version))

    try:
        latest = fetcher()
        origin = 'pypi'
    except Exception as e:
        logger.debug(f'update check: could not reach PyPI: {e}')
        latest, origin = cached_latest, 'unavailable'
    cache['last_check'] = now
    if latest:
        cache['latest'] = latest
    _write_cache(cache_path, cache)
    return UpdateCheckResult(current=current_version, latest=latest, origin=origin,
                             update_available=_is_newer(latest, current_version))


def emit_update_notice(stream, enabled: bool = True, **kwargs) -> bool:
    """Write the update notice to ``stream`` if there is one, and return whether it did.

    This is the call site on the startup path, and the only thing it promises is that it cannot
    disturb the command the user actually ran: **every** exception is swallowed, including ones
    :func:`check_for_update` would let through.  A bug in the update check must not be able to
    fail a build.

    Returns ``False`` without checking anything -- without so much as a DNS lookup -- when
    ``stream`` is not a terminal.  See the module docstring for why that is a feature.
    """
    try:
        if not enabled or not callable(getattr(stream, 'isatty', None)) or not stream.isatty():
            return False
        notice = check_for_update(**kwargs).notice
        if notice is None:
            return False
        print(f'\n{notice}\n', file=stream, flush=True)
        return True
    except Exception as e:  # noqa: BLE001 -- see docstring; this one is deliberately total
        logger.debug(f'update check suppressed: {e}')
        return False
