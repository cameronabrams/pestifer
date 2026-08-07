# Design: bring an incoming (foreign) PSF into the pipeline

Status: **P1 implemented** (Unreleased). Roadmap item: *Ideas / unsorted → "Accept an incoming PSF
… (build onto an existing topology)."*

## Goal

Let a user bring a **pre-built topology** — a PSF plus its coordinates, from CHARMM-GUI, a prior
pestifer build, or another tool — into pestifer's downstream machinery (membrane equilibration,
production, packaging) **without rebuilding**, and later **extend** it (add a ligand, graft a glycan,
patch, mutate).

Two workflows:

- **Interop / round-trip** (P1) — take a foreign system straight into pestifer's
  `md` / `membrane_equilibrate` / `terminate` pipeline.
- **Modify a finished build** (P2/P3) — reopen a foreign or finished system and edit its topology.

## Architecture: state enters through the pipeline, never a psfgen source key

**Decision (superseding an earlier draft):** an incoming PSF must enter through a task that *provides*
the pipeline's `STATE` currency — not through a `source`-level key on `psfgen`. `psfgen`'s contract is
"consume `SOURCE`, build a system"; it is not a state door. Keeping a single data-flow model for where
systems come from (`fetch`→`psfgen`, `continuation`, `merge`) is the point.

The contract-correct entry is **`continuation`** (`pestifer/tasks/continuation.py`, a subclass of
`PsfgenTask`). It already:

- reads an external `psf` + `pdb`/`coor` + optional `xsc`/`vel` from disk and copies them into the run
  directory;
- converts `.coor`→`.pdb` via `coor_to_pdb`;
- resolves the topology stream files the PSF records — **explicitly handling PSFs built outside
  pestifer** (warns and drops an unshipped stream file rather than aborting);
- registers the `STATE` fileset (`provides=(STATE,)`, `discards_state=True` — an origin task).

So the pass-through half of the interop workflow already worked; the entry point is not `fetch`
(which provides `SOURCE` = raw coordinates for `psfgen` to build), because a PSF is *already-built* =
`STATE`.

## P1 — what was added to `continuation` (Unreleased)

A `verify_parameters` flag (schema, **default `True`**) gates a one-time parse of the incoming PSF
that does two things a bare state hand-off never did:

1. **Composition report + unknown-resname audit.** Classify each residue's segtype from its resname
   (`Labels.segtype_of_resname`), log a per-segid inventory (segid → segtype, residue count, sample
   resnames), and **warn** about any resname unknown to pestifer's segtype table (it classifies as
   blank, which weakens downstream protein/lipid/water/ion selections — the topology is still carried
   through).

2. **Force-field-consistency preflight (hard-error).** Verify the build's CHARMM release resolves
   **every atom type and every bond/angle/dihedral/improper type-tuple** in the PSF, and hard-error
   up front — naming the residue and term — on a miss. This turns the most likely real failure (a
   CHARMM-GUI PSF built against a *different* CHARMM version) from a cryptic mid-NAMD abort
   (`DIDN'T FIND vdW PARAMETER FOR ATOM TYPE …`) into an actionable message. The parameter set is the
   release standard params (via the NAMD scripter's `fetch_standard_charmm_parameters`) merged with
   any parameter/stream files the PSF records (already staged by `continuation`).

pestifer's **own internal continuations** — `make_membrane_system` (the quilt/patch continuation),
`make_pdb_collection`, and `make_solvent_box` — pass `verify_parameters=False`, preserving the
zero-parse fast path on systems pestifer itself built (which are trivially consistent, and can be
large or numerous).

### Components

- **`pestifer/charmmff/psf_param_check.py`** (new, NAMD-free, unit-tested) —
  `check_psf_parameters(psf, param) -> MissingParameters` and `format_missing(...)`. Matching mirrors
  CHARMM's own lookup: exact atom-type membership (vdW); unordered bonds; reversible angles;
  reversible dihedral quartets with the `X A B X` central-pair wildcard indexed and rarer
  partial-wildcards matched directly; permissive `X`-wildcard improper matching in either direction.
  Terms are **deduplicated by atom-type tuple**, so the scan is cheap even on a multi-million-atom
  system. Validated: **zero false positives** on a real pestifer-built BPTI PSF against the full feb26
  parameter set (all atom types, bonds, angles, dihedrals, impropers resolve).

- **`continuation._verify_incoming_topology` / `._report_incoming_composition`** — orchestrate the
  parse, report, param merge, and hard-error.

- **`PSFAtom._from_psf_atomline`** — the resname→segtype lookup changed from a hard `[...]` (which
  `KeyError`-crashed the whole parse on any unknown resname) to `.get(resname, '')`. A foreign PSF may
  legitimately carry a custom-ligand resname; it now classifies as blank rather than crashing. This is
  behavior-preserving for every known resname.

### Usage

```yaml
tasks:
  - continuation:
      psf: system.psf
      pdb: system.pdb            # or coor + xsc
  - membrane_equilibrate: { temperature: 310 }
  - terminate: { basename: my_system, package: { basename: prod } }
```

### Testing

- **Unit** (`tests/unit/test_charmmff/test_psf_param_check.py`): synthetic param set + crafted minimal
  PSF — all terms resolve (incl. wildcard dihedral/improper); missing atom type is named with its
  residue; missing bond flagged; missing wildcard dihedral flagged; unknown resname parses (blank
  segtype) without crashing.
- **Integration** (`tests/unit/test_tasks/test_continuation.py`): the existing BPTI continuation now
  runs with the default-on check and passes; a PSF with a bogus atom type hard-errors (naming the
  term); `verify_parameters=False` skips the check.

## P2 / P3 — editing an incoming topology (future)

Editing stays consistent with the state-through-pipeline principle: `continuation` provides the
`STATE`, and a downstream `psfgen` consumes that `STATE` and edits it. That requires `psfgen` to gain
a **readpsf-preserve** consume-a-`STATE` mode (`load_project` → `readpsf`, `guesscoord=False`,
`regenerate=False`, as `merge`/`ligate`/`ring_check` already use) — *not* the segment-rebuild path,
which would re-derive topology from the PDB and discard the foreign PSF's patches/custom residues.
Note the existing internal `prebuilt` reconstruction in `psfgen` *does* re-segment, so it is not the
basis for this.

- **P2 — additive edits.** `patches`/`ssbonds`/`links`/`grafts` layered via readpsf + patch +
  `regenerate`; coordinate mods (orient/rotate) in numpy.
- **P3 — mutating edits.** `mutations`/`deletions`/`insertions` via per-chain re-segmentation surgery
  on the preserved topology; its own design pass.

Fetch-metadata mods (`biological_assembly`, `SEQADV` mutations, `REMARK 465` healing,
`terminal_tails`) are meaningless for a foreign PSF (no source metadata) and should hard-error in the
edit path.

## Open items

- **Unknown-resname policy** stays warn-and-proceed (readpsf carries the topology; a downstream task
  that needs to *select* those residues is where a missing segtype would bite). Revisit if a consumer
  needs a hard requirement.
- **Large-system cost.** The preflight parses the whole PSF once; fine for a one-time ingest guard and
  off for internal callers. If a user ingests a very large foreign system repeatedly, `verify_parameters:
  false` is the escape hatch.
