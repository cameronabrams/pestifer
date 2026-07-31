# Restart / resume an interrupted build

**Status:** design, not started. Written 2026-07-31. Substrate mapped (see anchors inline).

## Why

A long multi-task build that dies partway — Ctrl-C, a crash, a killed job, a wall-clock timeout —
currently starts over from scratch. The SIGINT/SIGTERM handler (`core/command.py:82-96`, v3.11.1)
already tears children down cleanly and tells the user *"partially built files remain … the build is
incomplete and should be restarted"* — but there is no restart path. Add one that resumes from the last
cleanly-completed task instead of re-running the whole pipeline.

## Substrate (what already exists)

Favorable:
- **Deterministic basenames.** `next_basename` = `{controller:02d}-{index:02d}-{sub:02d}_{taskname}`
  (`basetask.py:398-413`); `index` is list position (`taskcollections.py:35`). Same `tasks:` list →
  identical output filenames on re-run. (Position-sensitive: inserting/removing a task shifts every
  downstream index — see staleness below.)
- **`continuation` already "starts from a fileset."** `ContinuationTask` (`tasks/continuation.py`) is an
  origin (`requires=()`, `provides=(STATE,)`) that loads a `psf` + `coor`/`pdb` (+ `xsc`/`vel`), recovers
  CHARMM stream files from the PSF `REMARKS topology` lines, and registers the `state` fileset. This is
  the natural resume primitive.
- **Seeded RNG throughout** (declash `_DECLASH_SEED`, CCD `seed`, grid packing `seed`, athermal MC
  `seed`) → a resumed tail is reproducible on a fixed machine (not bit-identical across platforms, by
  nature of wrapping NAMD).
- **Each task's resolved spec is reachable** as `task.specs` (`basetask.py:101`, after ycleptic
  `dwalk`) → hashable for staleness detection.

The gap:
- **The artifact graph is in-memory only.** `PipelineContext` keeps `head`/`history` (`core/pipeline.py`)
  but never serializes them. After an interrupt the *files* are on disk (named by basename) but nothing
  records *which task produced them* or *which fileset was the current `state`*. `do_tasks`
  (`controller.py:129-162`) is a plain loop that `break`s on failure and records nothing durable. **This
  is the substrate we must add: a persisted per-task completion manifest.**

## Decisions (the roadmap pre-committed these; settled here)

1. **Granularity: per-task.** Resume at the first task whose recorded outputs are missing or whose
   config changed. Within-task resume (e.g. mid-NAMD) is out of scope — NAMD already has its own
   restart, and not every task is idempotent.
2. **Completion detection: a persisted manifest keyed by a hash of each task's resolved spec.** A task
   counts as "done and current" iff (a) a manifest entry exists for its (index, taskname), (b) the
   recorded spec-hash matches the current resolved spec, and (c) the files it recorded still exist.
3. **Invocation: `pestifer run <cfg> --restart` (auto-detect the resume point)** as the primary path,
   with `--from <index|taskname>` as an explicit override and `--fresh` to force-ignore a manifest.
4. **Determinism: lean on the existing seeds.** The resumed tail is byte-reproducible on the same
   machine; a config edit that changes a seed correctly invalidates everything downstream of it (via the
   spec-hash), so the resume is honest.

## The manifest

A single JSON file `.pestifer-manifest.json` in the build directory, written **incrementally and
atomically** (temp + `os.rename`) by `do_tasks` after each task returns `result == 0`.

```jsonc
{
  "pestifer_version": "3.14.0",
  "config_path": "build.yaml",
  "tasks_fingerprint": "<hash of the ordered [(taskname, spec_hash), ...]>",
  "tasks": [
    {
      "index": 2,
      "taskname": "md",
      "spec_hash": "<sha256 of the canonicalized resolved spec>",
      "state": {"psf": "00-01-00_psfgen.psf", "pdb": "00-02-00_md.pdb",
                "coor": "00-02-00_md.coor", "xsc": "00-02-00_md.xsc", "vel": "00-02-00_md.vel"},
      "provides": ["STATE", "MD_OUTPUT"]
    }
    // ... one entry per cleanly-completed task
  ]
}
```

- `spec_hash` = SHA-256 of `json.dumps(resolved_spec, sort_keys=True)`. **Open substrate check:** confirm
  `task.specs` is *fully* defaulted (nested/list defaults included) before hashing, or re-resolve via
  `config.make_default_specs('tasks', taskname)` merged with the user dict (Config/ycleptic §5 of the
  map). Exclude volatile keys (none known yet; audit).
- `state` = the pipeline's current `head['state']` fileset **after** the task (unchanged tasks like
  `mdplot`/`validate` carry the prior state forward — record whatever is current).
- `provides` = the task's contract currencies (`pipeline_contract(task.specs).provides`), recorded so
  the resume algorithm can detect an unrecoverable non-STATE handoff (below).
- The interrupted (partial) trailing task has **no** entry — entries are written only on clean
  completion — so it is naturally the resume point and its half-written files get overwritten.

## Resume algorithm (`--restart`)

1. Build the task list from the (possibly edited) config as normal.
2. If no manifest in cwd → normal run (with `--restart`, warn and proceed fresh).
3. If the manifest's `tasks_fingerprint` differs from the current list's, or `pestifer_version` differs
   → **do not trust index alignment**; fall back to the longest matching *prefix* by walking
   `(taskname, spec_hash)` pairs from the top and stopping at the first divergence (a task-list edit
   invalidates from the edit point down — conservative and correct).
4. Otherwise compute **k = the last index** such that for every `i ≤ k`: `manifest[i]` matches
   `tasks[i]` by taskname + spec_hash **and** every file in `manifest[i].state` exists on disk. Resume
   point = `k + 1`.
5. **Non-STATE handoff guard (v1 scope).** Only `STATE` is restorable from the manifest. If any tail
   task (`> k`) `requires` a currency other than STATE that is `provides`-d only by a skipped task
   (`≤ k`) — e.g. `mdplot` needs `MD_OUTPUT` from a skipped `md` — back the resume point up to *before*
   that provider (re-run it), or, if that is impractical, refuse with an actionable message. Most build
   pipelines are STATE-threaded (fetch→psfgen→md→solvate→md→…), so this is the uncommon edge; handling it
   by backing up is safe if occasionally over-conservative.
6. **Restore state.** Seed the pipeline `head['state']` from `manifest[k].state` and re-establish what
   `continuation` would (copy inputs into cwd if needed, recover + register CHARMM stream files from the
   PSF REMARKS). Implementation options — pick in P1:
   - **(a) Reuse `continuation` internals directly** (call the state-restore portion of
     `ContinuationTask.do()` against `manifest[k].state`) — no index shift, keeps the tail's basenames
     identical to a fresh run. *Preferred.*
   - (b) Synthesize a `continuation` task as a new task-0 — simplest, but shifts every downstream index →
     the resumed tail gets *different* basenames than a fresh run. Rejected unless (a) proves fragile.
7. Run tasks `k+1 … end` normally (same indices → same basenames), continuing to append the manifest.
8. Proactively delete any on-disk files matching the resume task's basename glob before it re-runs, so a
   half-written partial from the interrupt can't be mistaken for a completed output.

`--from <index|taskname>` skips straight to that point (still restoring state from the manifest entry
just before it, requiring that entry to exist). `--fresh` deletes/ignores the manifest.

## Interaction with `terminate`

The default `terminate` task tars all non-`keep` file artifacts and **deletes the originals**
(`terminate.py` cleanup). So a **successfully completed** build has no loose intermediates to resume
from — and doesn't need to. Resume therefore only applies to builds that never reached `terminate`.
On a clean `terminate`, remove the manifest (or mark it `complete`) so a subsequent `--restart` on the
same directory doesn't try to resume a finished build.

## Scope / limitations (v1)

- **STATE-threaded resume only.** Non-STATE currency handoffs across the resume boundary are handled by
  backing the resume point up (§5), not by persisting every currency. Persisting the full artifact head
  is a possible v2 if a real pipeline needs it.
- **Config edits invalidate from the edit point down** (by spec-hash / fingerprint). This is intended.
- **Cross-platform determinism is not guaranteed** (NAMD). Same-machine is.
- Subcontroller tasks (`make_membrane_system`) run their own inner pipeline; v1 treats them as a single
  opaque task (resume at the outer boundary only). Inner-pipeline resume is a possible follow-on.

## Phasing

- **P1 — the manifest + `--restart` happy path.** Write the manifest incrementally in `do_tasks`;
  add `--restart` that computes the resume point for an unedited config, restores STATE via
  continuation-internals (6a), and runs the tail. Unit-test the resume-point computation (pure function
  over manifest + task list + a mocked filesystem) and the manifest round-trip. End-to-end: interrupt a
  small build (fetch→psfgen→md→solvate→md), `--restart`, confirm it skips the completed prefix, produces
  the same final artifacts as an uninterrupted run, and is byte-identical on the resumed tail.
- **P2 — staleness + edits.** spec-hash / tasks-fingerprint invalidation, `--from`, `--fresh`, the
  non-STATE handoff guard (§5), and partial-file cleanup (§8).
- **P3 — polish.** Manifest cleanup on `terminate`; `--restart` messaging (what was skipped, where it
  resumed); optional `pestifer run --status` to print the manifest/resume plan without running.

## Testing

- Pure resume-point unit tests (manifest × task-list × file-existence matrix): clean resume, edited
  task, inserted/removed task, missing state file, non-STATE handoff back-up.
- Manifest atomic-write / round-trip.
- End-to-end interrupt+resume on a fast soluble example, asserting final-artifact equality and
  resumed-tail reproducibility.
