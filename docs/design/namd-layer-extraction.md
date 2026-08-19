# Checkpoint: extracting a NAMD layer from pestifer

Status: **not started, not decided.** This is a handoff note for a fresh session, written
2026-08-19. It records the evidence gathered so far, the shape the extraction would take, and the
questions that have to be answered before any code moves. Nothing here has been implemented.

## Where the idea came from

The question was "is pestifer getting too big?", prompted by noticing that it ships subcommands
mixing analysis of *builds* with analysis of *production runs it did not manage* —
`make-namd-restart` being the clearest case.

The first measurement mostly disproved the framing. Of ~48,600 lines, the code serving runs
pestifer did not manage is tiny, and the mass sits in `tasks` (8.5k), `core` (7.0k), `charmmff`
(6.9k) and `molecule` (6.1k) — all core mission. **Pestifer is not big because of scope creep.**

What *is* real is a coherent NAMD-shaped layer inside it that other NAMD users would want and that
has nothing to do with system preparation. That is what this note is about. The precedent is
strong and local: **`pidibble` and `ycleptic` are both prior extractions from this codebase**, both
independently useful, and both worked.

## What would move, and what cannot

This is the finding that matters, and it is sharper than the original "extract namdtools" idea.
There are two populations, and only one is extractable.

### Extractable: NAMD file formats and analysis (~1,000 statements)

A layer that reads what NAMD writes and computes things from it. Its coupling back into pestifer
is light or nil:

| Module | Stmts | Imports from pestifer |
|---|---|---|
| `util/densityprofile.py` | 128 | **nothing** |
| `util/density_convergence.py` | 362 | `densityprofile` only (internal to this group) |
| `logparsers/namdlogparser.py` | 396 | `logparser` base, `util/progress`, `util/stringthings` |
| `util/namdrestart.py` | 126 | `namdlogparser`, `scripters/generic`, `stringthings` |

Between them: parse `.log` and `.xst`, compute density, lateral area and area-per-lipid,
autocorrelation-corrected convergence assessment, density profiles along an axis, and generate a
restart configuration from a checkpoint. Every one of those is meaningful to someone who has never
heard of pestifer.

Three small shared pieces would have to come along or be inlined: `logparsers/logparser.py` (289
lines, the generic base), `util/progress.py` (88 lines, of which `NAMDProgress` is the relevant
subclass), and parts of `util/stringthings.py`.

### Not extractable: composing a NAMD configuration

`scripters/namd.py` (193 stmts) is the piece one would most want to move and the one that cannot,
because of what it imports:

```
from ..charmmff.charmmffprm import CharmmParamFile
from ..psfutil.psfcontents import PSFContents
from ..core.command import Command
```

Writing a NAMD config for a **CHARMM system built by pestifer** is inherently pestifer's job — it
needs the force field and the PSF. Likewise `scripters/namd_colvar_input.py` pulls in
`molecule/atom`.

So the honest split is:

> **reading and analysing NAMD output** is a general capability worth extracting;
> **writing NAMD input for a CHARMM system** is pestifer's own domain and stays.

## What consumes the extractable layer today

```
pestifer/subcommands/density_profile.py
pestifer/tasks/{mdtask, mdplot, density_equilibrate, equilibrate_base, membrane_equilibrate}.py
pestifer/tasks/make_membrane_system.py (via the convergence criterion)
```

Note this means pestifer would **depend on** the new package, exactly as it depends on pidibble —
it does not shed the code, it consumes it from outside. That is the correct model and it should be
stated plainly, because "pestifer is too big" instincts sometimes expect the line count to drop.
It would drop by roughly 1,000 statements out of ~23,800, i.e. about 4%.

## The case for, and the case against

**For.** The layer is coherent and independently useful. NAMD log parsing, `.xst` handling, an
honest autocorrelation-corrected convergence criterion, and restart generation are things NAMD
users write badly and repeatedly. The convergence criterion in particular is *good* — it is
autocorrelation-corrected, documented in `docs/design/density-equilibrate.md`, and validated
against an every-step reference. It is arguably the most reusable thing in this repository. The
extraction also forces the interface to be honest, which is what happened with pidibble.

**Against.** Three repositories instead of two means a third release cycle, a third version pin,
and a compatibility matrix. Two of the four modules are also *build* dependencies, so pestifer
consumes its own extraction rather than being freed of it. And the timing is bad while a
manuscript pins a version and a sweep is pending.

## The decision test

Not "would this make pestifer smaller" — it barely would. The test is:

> **Would you want this package if pestifer did not exist?**

If the answer is yes, extract it; the reuse argument carries it, exactly as it did for pidibble. If
the answer is only "it would tidy pestifer up", don't — 4% is not worth a third release cycle.

That question is about future work and is the user's to answer, not something to reason into from
the code.

## If it goes ahead: sequencing notes

1. **Take the leaves first.** `densityprofile` has zero pestifer imports and `density_convergence`
   depends only on it. Those two alone (490 statements, 94% and 74% covered) are a viable v0.1.0
   and prove the interface without touching the log parser.
2. **`namdlogparser` needs its base class and progress reporting resolved** — decide whether the
   generic `LogParser` goes too (it also serves psfgen and pdb2pqr parsers, which stay) or whether
   the NAMD parser is rewritten standalone.
3. **Coverage is uneven and should be levelled before the move, not after.**
   `density_convergence` is at 94% but `namdrestart` is at **13%** — the least-tested module in the
   candidate set, and the one whose behaviour is most obviously public API for an external user.
4. **Naming.** `namdtools` is taken on PyPI territory-wise; check before committing to a name.
5. **The extracted package should keep the design docs** — `density-equilibrate.md` and
   `membrane-equilibrate.md` document the criterion, not the pipeline, and would belong with it.

## What has changed since the first analysis (2026-08-18)

Worth knowing, because it adds NAMD-adjacent surface that did not exist when the question was
first asked:

- `util/provenance.py` probes the NAMD and VMD binaries for their versions — a small NAMD-aware
  piece that would plausibly belong in the extracted package.
- `MDTask._namd_seed` derives a per-invocation NAMD RNG seed from a build-wide base. NAMD-specific
  but pipeline-coupled; probably stays.
- `core/run_record.py` records what each MD task actually ran. Pipeline, stays.
- `equilibrate_base` and the density/membrane hooks are now unit-tested without NAMD (23% → 80%,
  and 42% → 73% for `membrane_equilibrate`), which materially lowers the risk of moving the
  convergence layer: its behaviour is now pinned by tests that do not need the toolchain.

## Related

- `docs/design/density-equilibrate.md` — the convergence criterion, the most reusable asset here.
- `docs/design/membrane-equilibrate.md` — the two-stage protocol and area convergence.
- `docs/design/roadmap.md` — "Split the stable CHARMM resources into a separate data package" is a
  *different* extraction idea, about the bundled force field rather than code. They are
  independent; do not conflate them.
