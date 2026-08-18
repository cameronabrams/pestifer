# Design: `report-methods` — a drafted Methods section from a completed build

Status: **P0-P3 implemented** (`run-record.json`, `methods.bib`, `methods.tex`/`methods.md`,
`pestifer report-methods`). Replica merging is in from the start, per the decision below. The
optional networked `--enrich` path is **not** implemented and remains as designed.

Analogue: GROMACS ships `gmx report-methods`, which emits a short LaTeX/plain-text summary of a
`.tpr` — system composition and run parameters. The pestifer version can say considerably more,
because a pestifer build knows not just the final system but *how it was produced*: which entry it
started from, what was modelled in, what force field, what equilibration actually ran, under which
software versions, from which seed.

## The contract, first

This feature writes prose that a scientist will put in a paper under their own name. That makes
the failure mode asymmetric: a *missing* sentence costs the author five minutes, while a
**plausible-sounding but wrong** sentence is a defect in the literature. Everything below follows
from that.

Three rules:

1. **Every generated sentence traces to a recorded fact from the run that produced it.** No
   sentence is generated from the config alone where the run could have diverged from it — and it
   routinely does: `density_equilibrate` and `membrane_equilibrate` run to a convergence criterion,
   so the number of steps they took is not in the config and is not knowable in advance.
2. **Pestifer never describes what it did not do.** It prepared a system; it did not run the
   production simulation, choose the ensemble for the science, or analyse anything. Those sections
   are emitted as explicit, visible `\TODO{...}` placeholders, never as prose.
3. **The output announces that it is a draft.** A header comment in the `.tex`, and a marker
   that has to be deleted deliberately.

Anti-goal: a button that produces a paragraph nobody reads before pasting.

## The blocker: there is no structured record of a run

The facts a Methods paragraph needs are, today, scattered and mostly unstructured:

| Fact | Where it lives now | Machine-readable? |
|---|---|---|
| Final atom count, box vectors, total charge | "System Report" block in the build log | no — prose |
| Steps actually run per MD task | the archived per-chunk `.namd` files (`run N`) | only by summing |
| Converged density / area, and the criterion met | build-log lines from the convergence gate | no — prose |
| Temperature, pressure, ensemble, timestep | archived `.namd` configs | partially |
| Solvent, ion species, salt concentration | the config | yes |
| Software versions, force-field release | the `environment:` block (new in 3.17.x) | no — prose |
| Citations owed | the `citations:` block (new in 3.17.x) | no — prose |
| Task list, task state files | `.pestifer-manifest.json` | yes |
| Seed | the config, plus NAMD logs | yes |

So the first phase is not the report. It is **recording the facts**.

## Phase 0 — a structured run record (`run-record.json`)

Written by `terminate` alongside the existing outputs. One JSON object per build:

```jsonc
{
  "pestifer_version": "3.18.0",
  "config": { "title": "...", "path": "bpti1.yaml", "sha256": "..." },
  "environment": { /* exactly what util.provenance.environment_report() already returns */ },
  "citations":  { /* exactly what util.citations already computes, software + coordinates */ },
  "seed": 27021972,
  "system": {
    "atoms": 14130, "bonds": 9799, "total_charge": 0.0,
    "box": {"a": 55.97, "b": 49.40, "c": 50.02, "shape": "orthorhombic"},
    "composition": {"protein": 892, "water": 4358, "ions": {"SOD": 6, "CLA": 5}}
  },
  "provenance": {
    "sources": [{"id": "6pti", "origin": "rcsb", "role": "starting structure"}],
    "modifications": {"mutations": [...], "grafts": [...], "ssbonds": [...]}
  },
  "protocol": [
    {"index": 3, "task": "md", "ensemble": "minimize", "steps": 1000},
    {"index": 7, "task": "density_equilibrate", "ensemble": "NPT",
     "steps": 82210, "temperature": 300.0, "pressure": 1.0,
     "stopped_because": "converged", "final_density": 1.0299,
     "criterion": "drift 1.77e-03 < tol 2.00e-03, three consecutive"}
  ]
}
```

Two things to note. `environment` and `citations` are *already computed* — this phase mostly
gives them a machine-readable home rather than only a log line. And `protocol[].steps` /
`stopped_because` must be **recorded by the tasks as they finish**, not reconstructed afterward;
that is the actual code change, and it is small (each MD task already knows what it ran).

**Phase 0 is independently worth doing** even if the Methods report is never built. It is what
makes a run queryable — for the manuscript's own sweep table, for regression comparison across
replicas, for `mdplot`. I would not bundle it into the report feature's justification.

## Phase 1 — `methods.bib`

Two populations of citation, with different problems.

**Software citations** are a fixed, curated set. `util.citations` already holds them as rendered
strings; this phase adds a hand-written BibTeX entry per catalog entry with a stable key
(`pestifer:namd2005`, `pestifer:charmm36m`, …). Fully under our control, correct by construction,
no network.

**Coordinate citations** are the hard half, and the reason the current report prints DOI and PMID
rather than reference strings: a PDB-format file stores authors and titles in upper case and the
original capitalization is not recoverable (`A.B.MCDERMOTT` → `McDermott, A.B.` is a guess, and
guessing sentence case from an ALL-CAPS title mangles `HIV-1`). Options, in order of preference:

1. **`@misc` with the DOI only** (default). Cannot be wrong; the author's reference manager or
   journal style resolves the rest from the DOI. Honest and offline.
2. **Full entry when the build fetched mmCIF**, which carries real capitalization. Free when
   `fetch.source_format: cif` — and worth a documentation nudge toward `cif` for this reason.
3. **`--enrich` (opt-in, networked)**, resolving the DOI via the RCSB Data API or Crossref for a
   properly-cased title, journal and page range. Never the default: pestifer's stated position is
   that a build touches no web service unless asked.

An `@software` entry for pestifer itself, keyed to the version and Zenodo DOI, comes from the same
place the Data Availability statement does.

## Phase 2 — `methods.tex`

A **fragment**, not a document: no preamble, `\input`-able, so it drops into an existing
manuscript. Sections, each emitted only if the run supports it:

- *System preparation* — starting coordinates and their provenance (entry, assembly, chains kept
  or excluded), then what was modelled: loops closed and by what method, mutations reverted or
  introduced, disulfides, grafts, fusions, cleavages. Every clause backed by a `modifications`
  entry, so a build with no mutations emits no mutation sentence.
- *Force field and solvation* — CHARMM36 release, water model, ion species and concentration, box
  geometry and final atom count. The non-aqueous case says so, and says whether the box shipped or
  was generated on demand.
- *Membrane*, when present — composition per leaflet, symmetric or asymmetric, area per lipid,
  how the protein was positioned.
- *Equilibration* — the ladder that actually ran, with real step counts, and for the adaptive
  tasks the criterion that was met, not the ceiling that was configured.
- *Software and reproducibility* — versions from the environment record, the seed, the input file
  name and hash, and a pointer to the archived release.
- `\TODO{Describe the production simulation: ensemble, length, timestep, ...}` — deliberately
  ugly, deliberately not prose.

Also emit **`methods.md`** from the same fact set (many users are not writing LaTeX) and echo the
underlying `run-record.json`. The generator must have exactly one source of truth; three renderers
over one record, never three independent template sets.

## Phase 3 — integration

- **`pestifer report-methods`**, run in a completed build directory (or pointed at a
  `run-record.json`). Naming deliberately echoes `gmx report-methods`.
- Not emitted by default at `terminate`. `terminate` writes `run-record.json`; the prose is
  produced on request. Generating a Methods paragraph unbidden invites exactly the
  paste-it-unread failure this design is trying to avoid.
- Flags: `--format tex|md|json`, `--output`, `--enrich`.

## Open questions

1. ~~**Replicas.**~~ **Decided: merge.** `report-methods` takes a list of records and describes
   them as a set. Facts every replica agrees on are stated once; a fact they disagree on is
   reported as a range (final box dimensions) rather than dropped, because silently omitting a
   fact the reader needs is the wrong resolution of a disagreement.
2. ~~**How much prose is too much?**~~ **Decided: visibly machine-generated.** Short declarative
   sentences, no connective flourish, no rationale. Fluency reads as authority and invites being
   pasted unread; text that obviously came from a tool invites the editing it needs.
3. **Verification.** Implemented as proposed: the tests assert on statements, not style — that a
   claim appears only when the record supports it, and that a claim which would be false does not
   appear. The test that earns its keep is the one asserting an unconverged run is never described
   as converged; the first implementation got that wrong.

   Still open: nothing compile-checks the generated LaTeX, because this machine has no TeX
   installation. The tests check brace balance and that task-name underscores are escaped, which
   caught the real defect, but that is not the same as compiling.
4. **Scope of the citation set.** Does the Methods draft cite everything the run touched, or only
   what a reader needs? Probably the former, since the author can delete — but it means a small
   BPTI build arrives with six software citations, which may read as padding.

## Why not just a template in the docs

Considered and rejected as the whole answer, but worth keeping alongside: a documented worked
example of a Methods paragraph costs almost nothing and helps immediately. The generated version
earns its complexity only because the adaptive equilibration means the true numbers are not
knowable from the config — which is exactly the part an author would otherwise get wrong by hand.
