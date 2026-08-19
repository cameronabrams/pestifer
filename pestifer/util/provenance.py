# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Environment provenance: what actually ran, recorded where the results are.

A pestifer build is a chain of external programs -- VMD, psfgen, NAMD -- driven by a Python
process, against a particular CHARMM force-field release.  Reproducing a result, or reporting it
in a paper, means knowing every version in that chain.  Until this module existed, the log
recorded only the CHARMM release, so a build's environment had to be reconstructed after the fact
from whatever happened to still be installed -- which is not evidence.

:func:`log_environment` writes the whole chain into the log at the start of a build, so any build
log is self-describing.  Nothing here may raise: a failed probe records ``unknown`` and the build
proceeds.  Provenance is worth a second of startup, never a failed run.
"""

import importlib.metadata
import logging
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile

from .stringthings import __pestifer_version__

logger = logging.getLogger(__name__)

_PROBE_TIMEOUT = 120
"""Seconds any single version probe may take before it is abandoned as ``unknown``.

Generous on purpose.  ``vmd --version`` still performs a full startup -- plugin registration and
GPU detection included -- so on a loaded machine or a cold cache it can take far longer than the
work suggests.  A tight timeout here does not save meaningful time; it just silently degrades the
provenance record, which is the one thing this module exists to prevent."""

_UNKNOWN = 'unknown'


def _run(cmd, cwd=None):
    """Run ``cmd``, returning combined output as text, or ``''`` on any failure.

    Version probes are advisory, so a nonzero exit is not an error: NAMD in particular reports its
    banner and *then* complains about the deliberately empty config it was handed.
    """
    try:
        p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           timeout=_PROBE_TIMEOUT, text=True, errors='replace')
        return p.stdout or ''
    except subprocess.TimeoutExpired:
        # Warn, not debug: the record the user keeps beside their results is now incomplete, and
        # an 'unknown' with no explanation is indistinguishable from a missing program.
        logger.warning(f'version probe {cmd[0]!r} timed out after {_PROBE_TIMEOUT}s; '
                       f'it will be recorded as unknown')
        return ''
    except Exception as e:
        logger.warning(f'version probe {cmd[0]!r} failed ({e}); it will be recorded as unknown')
        return ''


def stamp(seed=None):
    """The one-line provenance mark pestifer writes into the files it authors.

    Deliberately version-and-seed, not a wall-clock timestamp.  A timestamp says nothing about
    *what* produced a file and makes byte-identical regeneration impossible; the version and seed
    say exactly what produced it and are deterministic, so re-running one config with one pestifer
    reproduces its generated scripts byte for byte.  The build time is recorded in the log and in
    ``run-record.json``, where it does not cost that.

    Never applied to files pestifer *copies* rather than writes -- the CHARMM ``.rtf``/``.str``/
    ``.prm`` files must stay byte-identical to the release they came from, since that is how a
    reader verifies they are unmodified.
    """
    mark = f'pestifer {__pestifer_version__}'
    # _namd_seed returns 0 for 'no seed configured', so falsy means absent, not 'seed 0'
    return f'{mark}  seed {seed}' if seed else mark


def stamp_figure(fig, text=None, seed=None):
    """Write the provenance mark into the bottom-right corner of a figure.

    Figures are the artifacts most likely to travel without their build directory -- into a slide,
    a notebook, a manuscript -- so they are the ones that most need to say what made them.  The
    seed is included because a triplicate sweep otherwise produces three near-identical plots that
    nothing distinguishes.
    """
    fig.text(0.995, 0.005, text or stamp(seed), ha='right', va='bottom',
             fontsize=5.5, color='0.55')


def python_environment():
    """Interpreter identity: implementation, version, executable, and platform."""
    return {
        'implementation': platform.python_implementation(),
        'version': platform.python_version(),
        'executable': sys.executable,
        'platform': platform.platform(),
        'node': platform.node(),
    }


def _dist_names():
    """Distribution names pestifer declares as runtime requirements.

    Read from installed metadata rather than a hand-kept list, so the report cannot drift out of
    step with ``pyproject.toml``.  Entries guarded by an extra (``; extra == "test"``) are skipped:
    they are not part of a build's environment.
    """
    try:
        reqs = importlib.metadata.requires('pestifer') or []
    except Exception as e:
        logger.debug(f'could not read pestifer requirements: {e}')
        return []
    names = []
    for r in reqs:
        if 'extra ==' in r:
            continue
        m = re.match(r'^\s*([A-Za-z0-9._-]+)', r)
        if m:
            names.append(m.group(1))
    return sorted(set(names), key=str.lower)


def package_versions():
    """Resolved versions of pestifer itself and of every runtime dependency it declares.

    pestifer's own version comes from the same resolver the rest of the package uses, not from
    ``importlib.metadata``: metadata records the version at install time, so every editable
    install reports a stale number the moment a release is cut.  A provenance record that lied
    about the version of the thing producing the results would be worse than none.
    """
    out = {'pestifer': __pestifer_version__}
    for name in _dist_names():
        try:
            out[name] = importlib.metadata.version(name)
        except Exception:
            out[name] = _UNKNOWN
    return out


def vmd_version(cmd):
    """Version string reported by the ``vmd`` at ``cmd``, e.g. ``2.0.0 (March 25, 2026)``."""
    m = re.search(r'VMD for \S+, version (.+?)\s*$', _run([cmd, '--version']), re.M)
    return m.group(1).strip() if m else _UNKNOWN


def namd_version(cmd):
    """Version and build flavor of the ``namd3`` at ``cmd``.

    NAMD prints its banner only once it has a config to read, so it is handed an empty one in a
    temporary directory; it announces itself, objects to the config, and exits.  The CUDA build
    additionally reports the CUDA version it was compiled against, which is worth keeping -- a
    GPU-resident result depends on it.
    """
    with tempfile.TemporaryDirectory() as td:
        probe = os.path.join(td, 'version_probe.namd')
        try:
            open(probe, 'w').close()
        except Exception:
            return _UNKNOWN
        out = _run([cmd, '+p1', probe], cwd=td)
    m = re.search(r'^Info: (NAMD \S+ for \S+)\s*$', out, re.M)
    if not m:
        return _UNKNOWN
    version = m.group(1)
    cuda = re.search(r'^Info: Built with CUDA version (\S+)\s*$', out, re.M)
    if cuda:
        version = f'{version} (CUDA {cuda.group(1)})'
    return version


_PROBERS = {'vmd': vmd_version, 'namd3': namd_version, 'namd3gpu': namd_version}
"""Which resolved shell commands can report a version, and how.  ``charmrun`` and ``catdcd``
report none, so they are recorded by resolved path alone."""


def executable_versions(config):
    """Resolved path and reported version of each external program this build may invoke.

    ``namd3`` is always probed -- even a GPU build runs it for any task carrying
    ``cpu-override`` -- while ``namd3gpu`` is probed only in GPU mode, since probing the CUDA
    build binds a GPU device and a CPU run will never touch it.
    """
    out = {}
    gpu = getattr(config, 'namd_type', 'cpu') == 'gpu'
    for name, cmd in sorted(getattr(config, 'shell_commands', {}).items()):
        if name == 'namd3gpu' and not gpu:
            continue
        path = shutil.which(cmd) or _UNKNOWN
        prober = _PROBERS.get(name)
        version = prober(cmd) if (prober and path != _UNKNOWN) else ''
        out[name] = {'path': path, 'version': version}
    return out


def charmmff_release(config):
    """The CHARMM force-field release directory name this build resolves to."""
    try:
        return config.RM.charmmff_version_path(config['user']['charmmff'].get('release', '')).name
    except Exception as e:
        logger.debug(f'could not resolve charmmff release: {e}')
        return _UNKNOWN


def environment_report(config):
    """Assemble the full environment record as a plain dict (JSON-serializable)."""
    return {
        'python': python_environment(),
        'packages': package_versions(),
        'executables': executable_versions(config),
        'charmmff': charmmff_release(config),
        'namd_mode': getattr(config, 'namd_type', _UNKNOWN),
    }


def log_environment(config, log=None):
    """Write the environment record to the log, and return it.

    Emitted at INFO so it reaches the build log a user keeps beside their results, not only the
    debug-level diagnostics log.  Laid out one fact per line and prefixed uniformly so the block
    can be recovered from a log with a single grep.
    """
    log = log or logger.info
    try:
        report = environment_report(config)
    except Exception as e:                      # provenance must never break a build
        logger.debug(f'environment report failed: {e}')
        return {}
    py = report['python']
    log('environment: ---- build environment ----')
    log(f"environment: python      {py['implementation']} {py['version']}  ({py['executable']})")
    log(f"environment: platform    {py['platform']}  on {py['node']}")
    log(f"environment: charmmff    {report['charmmff']}")
    log(f"environment: namd mode   {report['namd_mode']}")
    for name, info in report['executables'].items():
        detail = f"{info['version']}  [{info['path']}]" if info['version'] else f"[{info['path']}]"
        log(f'environment: {name:<11} {detail}')
    pkgs = report['packages']
    log(f"environment: pestifer    {pkgs.get('pestifer', _UNKNOWN)}")
    deps = ', '.join(f'{k} {v}' for k, v in pkgs.items() if k != 'pestifer')
    log(f'environment: packages    {deps}')
    log('environment: ---------------------------')
    return report
