# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Per-task completion manifest for a build run -- the durable substrate for in-stream restart.

A build's task pipeline hands artifacts off in memory only (see
:class:`~pestifer.core.pipeline.PipelineContext`); after an interrupt the output *files* survive on
disk but nothing records *which task produced them* or *which fileset is the current* ``state``.
:class:`RunManifest` closes that gap: the top-level controller writes ``.pestifer-manifest.json`` in
the build directory, appending one entry per **cleanly completed** task -- its index, taskname, a
hash of its resolved spec, and the ``state`` fileset current after it.  ``pestifer run --restart``
(P1b) reads this to find the last still-valid task and resume the tail from there.

The manifest is advisory bookkeeping: any failure to write it is logged and swallowed, never fatal to
the build.  Only the top-level controller (index 0) writes one; subcontroller pipelines
(e.g. ``make_membrane_system``'s inner quilt) are opaque, so their parent task appears as a single
manifest entry.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import tempfile

logger = logging.getLogger(__name__)

MANIFEST_NAME = '.pestifer-manifest.json'
_STATE_SLOTS = ('psf', 'pdb', 'coor', 'xsc', 'vel')


def spec_hash(specs) -> str:
    """Stable SHA-256 of a task's resolved spec dict (order-independent, str-coerced)."""
    return hashlib.sha256(
        json.dumps(specs, sort_keys=True, default=str).encode()).hexdigest()


def tasks_fingerprint(tasks) -> str:
    """Fingerprint of an ordered task list: hash of ``[(taskname, spec_hash), ...]``.

    A change here (task inserted/removed/reordered, or any task's spec edited) means the manifest's
    index alignment can no longer be trusted wholesale -- restart falls back to longest-matching-prefix.
    """
    pairs = [[t.taskname, spec_hash(t.specs)] for t in tasks]
    return hashlib.sha256(json.dumps(pairs).encode()).hexdigest()


def _state_files(state) -> dict:
    """``{slot: filename}`` for the non-empty members of a ``StateArtifacts`` (or ``{}``)."""
    out = {}
    if state is None:
        return out
    for slot in _STATE_SLOTS:
        art = getattr(state, slot, None)
        if art is not None:
            out[slot] = getattr(art, 'name', str(art))
    return out


class RunManifest:
    """The build's per-task completion record, persisted atomically as it grows."""

    def __init__(self, path: str, *, version: str = '', config_path: str = '',
                 fingerprint: str = ''):
        self.path = path
        self.data = {
            'pestifer_version': version,
            'config_path': config_path,
            'tasks_fingerprint': fingerprint,
            'complete': False,
            'tasks': [],
        }

    @classmethod
    def load(cls, path: str) -> 'RunManifest':
        """Load an existing manifest from ``path``."""
        with open(path) as f:
            raw = json.load(f)
        m = cls(path)
        m.data = raw
        return m

    def record(self, task, pipeline, task_spec_hash=None) -> None:
        """Append (or replace) the entry for a cleanly-completed ``task`` and persist.

        Records the current ``state`` fileset (unchanged tasks carry the prior state forward) and the
        task's provided pipeline currencies (so restart can detect a non-STATE handoff across the
        resume boundary).  Never raises -- a manifest failure must not fail the build.

        ``task_spec_hash`` must be the hash of the task's specs computed **before** it executed:
        several tasks mutate ``self.specs`` during ``do()`` (resolved defaults, artifact bookkeeping),
        so hashing here (post-execution) would not match the pre-execution hash a later ``--restart``
        recomputes.  The caller passes the pre-execution hash; we fall back to hashing now only if it
        is omitted.
        """
        try:
            state = pipeline.get_current_artifact('state')
            try:
                provides = sorted(str(c) for c in task.pipeline_contract(task.specs).provides)
            except Exception:
                provides = []
            entry = {
                'index': task.index,
                'taskname': task.taskname,
                'spec_hash': task_spec_hash if task_spec_hash is not None else spec_hash(task.specs),
                'state': _state_files(state),
                'provides': provides,
            }
            self.data['tasks'] = [e for e in self.data['tasks'] if e['index'] != task.index]
            self.data['tasks'].append(entry)
            self.data['tasks'].sort(key=lambda e: e['index'])
            self._save()
        except Exception as exc:
            logger.warning(f'run manifest: could not record task {getattr(task, "index", "?")}: {exc}')

    def resume_point(self, tasks) -> int:
        """The last task index in ``tasks`` that is *done and still current* per this manifest.

        Walks ``tasks`` from the top; task ``t`` counts iff a manifest entry with the same index
        exists, matches by ``taskname`` **and** ``spec_hash``, and every file it recorded still
        exists on disk.  Stops at the first divergence — so an edited spec, an inserted/removed task,
        or a missing/stale output invalidates from that point down.  Returns ``-1`` when nothing
        matches (fresh run) or the build is already ``complete`` (nothing to resume).  The caller
        resumes at ``resume_point + 1``.
        """
        if self.data.get('complete'):
            return -1
        by_index = {e['index']: e for e in self.data['tasks']}
        last = -1
        for t in tasks:
            e = by_index.get(t.index)
            if not e or e['taskname'] != t.taskname or e['spec_hash'] != spec_hash(t.specs):
                break
            if not all(os.path.exists(f) for f in e.get('state', {}).values()):
                break
            last = t.index
        return last

    def state_entry(self, index: int) -> dict:
        """The recorded ``{slot: filename}`` STATE fileset for a completed task index, or ``{}``."""
        for e in self.data['tasks']:
            if e['index'] == index:
                return dict(e.get('state', {}))
        return {}

    def mark_complete(self) -> None:
        """Flag the build as finished (a completed build has nothing to resume)."""
        try:
            self.data['complete'] = True
            self._save()
        except Exception as exc:
            logger.warning(f'run manifest: could not mark complete: {exc}')

    def _save(self) -> None:
        """Write the manifest atomically (temp file + ``os.replace``)."""
        d = os.path.dirname(os.path.abspath(self.path)) or '.'
        fd, tmp = tempfile.mkstemp(dir=d, prefix='.manifest-', suffix='.tmp')
        try:
            with os.fdopen(fd, 'w') as f:
                json.dump(self.data, f, indent=2)
            os.replace(tmp, self.path)
        except Exception:
            try:
                os.unlink(tmp)
            except OSError:
                pass
            raise
