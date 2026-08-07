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

## P2 / P3 — editing an incoming topology

Editing stays consistent with the state-through-pipeline principle: `continuation` provides the
`STATE`, and a downstream `psfgen` consumes that `STATE` and edits it:

```yaml
tasks:
  - continuation: { psf: sys.psf, pdb: sys.pdb }
  - psfgen:                       # sees STATE, no SOURCE -> readpsf-preserve
      mods:
        patches: [ ... ]
  - md: { ... }
```

`psfgen` gains a **readpsf-preserve** mode (`load_project` → `readpsf`, `guesscoord=False`,
`regenerate=` only when a mod changed topology) — *not* the segment-rebuild path, which would
re-derive topology from the PDB and discard the foreign PSF's patches/custom residues. (The existing
internal `prebuilt` reconstruction *does* re-segment, so it is not the basis for this.)

### Mode selection — infer from the pipeline (decided)

`psfgen` preserves when it follows a `STATE`-provider with **no** `SOURCE`; it rebuilds when it
follows `fetch`. No new user syntax. This makes the **pipeline contract pipeline-aware**, which today
it is not: `validate_pipeline` tracks an `available` currency set but calls
`task.pipeline_contract(specs)` with specs only. The change:

- Thread the `available` set into the contract call. All 17 `pipeline_contract(cls, specs)` methods
  share one signature; the low-churn route is for `validate_pipeline` to pass `available` **only** to
  contracts that declare the parameter (`inspect.signature` check), so just `psfgen`'s override gains
  `available=frozenset()`. (Uniformly updating all 17 is the alternative; more churn, no `inspect`.)
- `psfgen.pipeline_contract(specs, available)` returns:
  - **preserve** — `SOURCE ∉ available and STATE ∈ available` →
    `requires=(), provides=(STATE, MOLECULE), discards_state=False`.
  - **build** (unchanged) — else → `requires=(SOURCE,), provides=(STATE, MOLECULE),
    discards_state=True`. Neither present → the existing "no `fetch` precedes" error still fires.

At runtime `psfgen.do()` branches on the same condition (STATE present, no `base_coordinates`
SOURCE): the preserve branch runs `_psfgen_preserve()` instead of `ingest_molecules()` + segment
`psfgen()`. `MOLECULE` is built lazily (via `ensure_base_molecule`), as `continuation` already does —
patches operate on segids/resids directly and don't need the full parse; grafts do and trigger it.

### Staging

- **P2.0 — plumbing (no mods). [DONE, Unreleased]** Pipeline-aware contract (`validate_pipeline`
  passes `available` to contracts that declare it, via `inspect.signature`; `psfgen.pipeline_contract`
  branches preserve-vs-build) + `PsfgenTask._psfgen_preserve()` that `load_project`s the incoming
  PSF/PDB and writes a fresh state (`guesscoord=False, regenerate=False`). `psfgen.do()` branches on
  `_incoming_state_without_source()`. Mods present in preserve mode hard-error ("not yet supported").
  Verified: `continuation → psfgen` validates and **reproduces the system** atom-for-atom (14127→14127
  on BPTI); build path (fetch→psfgen) unregressed; the old "continuation→psfgen is malformed" contract
  test updated to assert it is now valid. Establishes the contract/branch every later mod rides on.
- **P2.1 — patches. [DONE, Unreleased]** `_psfgen_preserve` builds an `ObjManager` from `mods`,
  emits each `patches` entry via `Patch.return_TcL()` (`patch NAME SEG:RESID`; `use_after_regenerate`
  patches through `addpostregenerateline`), then `writescript(guesscoord=True, regenerate=True)` so a
  patch that adds atoms (e.g. a protonation) has them placed and the connectivity rebuilt. Any mod
  outside `_PRESERVE_SUPPORTED_MODS = {('seq','patches')}` hard-errors, naming it. Verified: `ASPP:A:3`
  on an incoming BPTI PSF protonates ASP A:3 (12→13 atoms, total +1); an unsupported `ssbonds` mod
  hard-errors.
- **P2.2 — ssbonds. [DONE, Unreleased]** `ssbonds` route to a `patch DISU c1:r1 c2:r2` on the loaded
  project (the ssbond's chainIDs are the PSF segids in preserve mode — no assembly transform to remap),
  before `regenerate`. Allowlist now `{('seq','patches'), ('topol','ssbonds')}`. Verified on BPTI
  (routing/acceptance/clean run; its cysteines are already bonded, so a free-thiol formation case would
  need a different fixture).
- **P2.3 — coord torsion rotations. [DONE, Unreleased]** `irotations`/`crotations` apply on the
  just-written state via `coormods()`. They run per biological-assembly image, so the preserve path
  builds the base molecule lazily first (`ensure_base_molecule()` → identity single-image assembly for
  a pre-built system) and sets `self.base_molecule`, then `coormods()` routes them through
  `_apply_python_crotations`. Allowlist gains `('coord','irotations')`, `('coord','crotations')`.
  Verified: `CHI1,A,3,60.0` on an incoming BPTI PSF pivots the ASP A:3 sidechain about CA–CB (CB
  fixed <0.05 Å; OD2 swings >1 Å), topology unchanged. This also establishes the lazy-molecule pattern
  P2.4 reuses. (Note: `orient` is a `manipulate`-task movetype, not a psfgen `mods` key, so it is not
  applicable here; the psfgen coord mods are the torsion rotations plus VMD `rottrans`.)
- **P2.4 — links. [DONE, Unreleased]** A `links` mod (shortcode `C1_R1_A1-C2_R2_A2`, no user patch)
  has its patch **resolved from residue geometry**: the lazy `base_molecule`'s residues are assigned
  (`LinkList.assign_residues` → `set_patchname` → `ic_reference_closest` over the ASN/SER/THR–glycan IC
  patterns), then emitted via the scripter's `write_links` using the identity assembly transform
  (single-image, so `chainIDmap` is identity). Unresolvable links (residues not found) are warned and
  skipped. Allowlist gains `('topol','links')`. Verified on a committed glycoprotein fixture: the
  ASN A:61 – glycan V:1304 attachment resolves to **NGLA** and emits `patch NGLA A:61 V:1304`.
- **P2.5 — grafts (adds topology).** Unlike the patch-on-existing-topology edits above, a graft **adds**
  a glycan (donor `segment{}` + `coordpdb` + linkage patch + regenerate/guesscoord), deep in the
  build-path segment machinery — its own milestone.

### Mod routing (allowlist)

Fetch-metadata mods — `biological_assembly`, `SEQADV` mutations, `REMARK 465` healing,
`terminal_tails` — are meaningless for a foreign PSF (no source metadata) and **hard-error** in the
preserve path. Re-segmenting mods (`mutations`/`deletions`/`insertions`) are **P3** (per-chain
re-segmentation surgery on the preserved topology) and until then error as "not yet supported" rather
than silently no-op.

- **P3 — mutating edits.** `mutations`/`deletions`/`insertions` via per-chain re-segmentation surgery;
  its own design pass.

## Open items

- **Unknown-resname policy** stays warn-and-proceed (readpsf carries the topology; a downstream task
  that needs to *select* those residues is where a missing segtype would bite). Revisit if a consumer
  needs a hard requirement.
- **Large-system cost.** The preflight parses the whole PSF once; fine for a one-time ingest guard and
  off for internal callers. If a user ingests a very large foreign system repeatedly, `verify_parameters:
  false` is the escape hatch.
