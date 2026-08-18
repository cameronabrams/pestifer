# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A machine-readable record of what a build actually did.

The build log has always carried these facts, but as prose: the final atom count sits in a
"System Report" table, the steps an equilibration ran sit in per-chunk log lines, and the software
versions and citations sit in their own blocks.  Anything downstream -- a Methods draft, a
comparison across replicas, a summary table for a paper -- had to scrape them back out.

This module collects them once, at the end of a build, into ``run-record.json``.

The distinction that matters is between *what was asked for* and *what happened*.  A config is a
request; ``density_equilibrate`` and ``membrane_equilibrate`` run to a convergence criterion, so
the number of steps they take is decided at run time and is not knowable in advance -- the same
config legitimately produces different step counts on different days.  A record derived from the
config would therefore state numbers the run never produced.  Everything under ``protocol`` here
is reported *by the tasks themselves* as they finish.

See ``docs/design/methods-report.md``.
"""

import json
import logging
import os

import numpy as np

from ..util.util import cell_from_xsc

logger = logging.getLogger(__name__)

RUN_RECORD_NAME = 'run-record.json'
"""Conventional filename, written into the run directory by the ``terminate`` task."""

RUN_RECORD_VERSION = 1
"""Schema version of the record itself, so a consumer can refuse a shape it does not know."""


def system_facts(state):
    """Topology and box facts for the final built system, or ``{}`` if it cannot be read.

    Reads the same PSF and XSC the ``terminate`` task reports on, so the record and the logged
    System Report cannot disagree.
    """
    from ..psfutil.psfcontents import PSFContents

    facts = {}
    if state is None:
        return facts
    psf = getattr(state, 'psf', None)
    if psf is not None and psf.exists():
        try:
            contents = PSFContents(psf.name)
            counts = {
                'atoms': contents.token_count.get('ATOM'),
                'bonds': contents.token_count.get('BOND'),
                'angles': contents.token_count.get('THETA'),
                'dihedrals': contents.token_count.get('PHI'),
                'impropers': contents.token_count.get('IMPHI'),
                'cross_terms': contents.token_count.get('CRTERM'),
            }
            facts.update({k: v for k, v in counts.items() if v is not None})
            facts['total_charge'] = round(float(sum(a.charge for a in contents.atoms)), 4)
        except Exception as e:
            logger.debug(f'could not read topology facts from {psf.name}: {e}')
    xsc = getattr(state, 'xsc', None)
    if xsc is not None and xsc.exists():
        try:
            box, _ = cell_from_xsc(xsc.name)
            if box is not None:
                facts['box'] = {
                    'a': round(float(np.linalg.norm(box[0])), 3),
                    'b': round(float(np.linalg.norm(box[1])), 3),
                    'c': round(float(np.linalg.norm(box[2])), 3),
                    'vectors': [[round(float(x), 3) for x in row] for row in box],
                }
        except Exception as e:
            logger.debug(f'could not read box facts from {xsc.name}: {e}')
    return facts


def protocol_from_tasks(tasks):
    """What each task actually ran, in order, as reported by the tasks themselves.

    A task contributes an entry only if it recorded one (see ``BaseTask.record_outcome``), so a
    task that does no simulation simply does not appear under ``protocol``.
    """
    protocol = []
    for task in tasks or []:
        outcome = getattr(task, 'outcome', None)
        if not outcome:
            continue
        entry = {'index': getattr(task, 'index', None),
                 'task': getattr(task, 'taskname', None)}
        entry.update(outcome)
        protocol.append(entry)
    return protocol


def system_facts_from_tasks(tasks):
    """The final system's facts, as captured by the task that had the state in hand.

    ``terminate`` reads them before packaging and cleanup archive the files; by the time a caller
    assembles the record those files are gone, so they cannot be re-read here.
    """
    for task in reversed(tasks or []):
        facts = getattr(task, 'system_facts', None)
        if facts:
            return facts
    return {}


def build_run_record(config, tasks, *, environment=None, citations=None):
    """Assemble the full record as a plain, JSON-serializable dict.

    ``environment`` and ``citations`` are passed in rather than recomputed, so the record carries
    exactly what the build already reported to its log, and the system facts come from the task
    that captured them rather than from files that no longer exist.
    """
    from ..util.stringthings import __pestifer_version__

    user = {}
    try:
        user = config['user']
    except Exception:
        pass
    namd = (user.get('namd') or {}) if isinstance(user, dict) else {}

    return {
        'run_record_version': RUN_RECORD_VERSION,
        'pestifer_version': __pestifer_version__,
        'title': user.get('title') if isinstance(user, dict) else None,
        'config_file': getattr(config, 'userfile', '') or None,
        'seed': namd.get('seed'),
        'namd_mode': getattr(config, 'namd_type', None),
        'environment': environment or {},
        'citations': citations or {},
        'system': system_facts_from_tasks(tasks),
        'protocol': protocol_from_tasks(tasks),
    }


def write_run_record(record, path=RUN_RECORD_NAME):
    """Write ``record`` as JSON.  Never raises: a build is not failed over its own summary."""
    try:
        with open(path, 'w') as f:
            json.dump(record, f, indent=2, sort_keys=False, default=str)
            f.write('\n')
        logger.info(f'wrote {path}')
        return path
    except Exception as e:
        logger.warning(f'could not write {path}: {e}')
        return None


def load_run_record(path):
    """Read a record written by :func:`write_run_record`.

    Raises ``ValueError`` for a shape this version does not understand, rather than silently
    misreading a record from a future pestifer.
    """
    with open(path) as f:
        record = json.load(f)
    version = record.get('run_record_version')
    if version is None:
        raise ValueError(f'{path} is not a pestifer run record (no run_record_version)')
    if int(version) > RUN_RECORD_VERSION:
        raise ValueError(f'{path} is a version {version} run record; this pestifer understands '
                         f'up to version {RUN_RECORD_VERSION}')
    return record


def load_run_records(paths):
    """Load several records -- the replica case.

    Accepts record files or directories containing one.  Returns them in the order given, so a
    caller can rely on replica order.
    """
    out = []
    for p in paths:
        candidate = os.path.join(p, RUN_RECORD_NAME) if os.path.isdir(p) else p
        out.append(load_run_record(candidate))
    return out
