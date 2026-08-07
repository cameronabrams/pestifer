# Design: accept an incoming PSF in the `psfgen` task

Status: **plan** (nothing implemented yet). Roadmap item: *Tooling / packaging → "Accept an
incoming PSF in the `psfgen` task (build onto an existing topology)."*

## Goal

Let a user bring a **pre-built topology** — a PSF plus its coordinates, from CHARMM-GUI, a prior
pestifer build, or another tool — into pestifer's downstream machinery (mods, membrane embedding,
equilibration, packaging) **without rebuilding the topology from scratch**, and optionally **extend**
it (add a ligand, graft a glycan, patch, mutate) later.

Two workflows this unlocks:

- **Interop / round-trip.** Take a CHARMM-GUI membrane-protein system straight into pestifer's
  `membrane_equilibrate` / production / packaging pipeline.
- **Modify a finished build.** Reopen a completed pestifer (or foreign) system and edit its topology
  without a from-scratch rebuild.

This is distinct from the existing `continuation` task, which runs MD on a *fixed* topology from a
state fileset. Here `psfgen` may *edit* the incoming topology.

## What already exists (substrate — we are not starting from zero)

Reconnaissance of the current tree (2026-08, `ccb43a40`):

1. **Two psfgen modes already coexist.**
   - **Build mode** — `PsfgenTask.psfgen()` → `PsfgenScripter.set_molecule` → `write_segments`:
     parses a structure into `Segment`s and *rebuilds* topology via `segment{ }` stanzas. Requires
     every residue to be a known CHARMM `RESI`. This is the normal fetch→build path.
   - **readpsf mode** — `PsfgenScripter.load_project(psf, pdb)` emits `readpsf <psf> pdb <pdb>`,
     loading an existing topology *verbatim*. Already used by `ligate`, `ring_check`, and — the
     closest cousin — **`merge`**, which reads external PSFs and merges them with
     `writescript(guesscoord=False, regenerate=False)` (merge.py:199-205). **This is the mechanism
     the incoming-PSF feature is built on.**

2. **A `prebuilt` source mode already exists** (`Molecule(source={'prebuilt': {psf, pdb, xsc}})`,
   molecule.py:121-128) — but today it is an **internal pipeline mechanism**: a psfgen task
   reconstructs its base molecule from a *prior pestifer STATE* (basetask.py:335, psfgen.py:838-843,
   919-951). It parses the *PDB* for structure and carries the PSF only as a file reference; when it
   then runs `psfgen()` it goes through **build mode** (re-segment + rebuild), which is fine only
   because pestifer itself produced that PSF (all RESIs/patches are known and re-derivable). It does
   **not** preserve a foreign topology. So the name is taken, but the capability we want is not there.

3. **Segtype classification is resname-driven** — `ResidueList.apply_segtypes()` maps
   `resname → Labels.segtype_of_resname` (residue.py:939). This works on residues parsed from a PDB
   regardless of origin, so classifying an incoming system works **as long as its resnames are known**
   to pestifer's label tables. Unknown resnames classify as `''` (blank) today — a gap to handle.

4. **Missing-parameter → actionable error** machinery exists for solvent boxes
   (`make_solvent_box`; see `on-demand-solvent-parameters.md`): it turns a cryptic mid-run NAMD
   parameter abort into an up-front error naming the exact residue and atom-type term. The
   force-field-consistency check below should reuse/generalize that, not reinvent it.

## Core design decision: **preserve, don't re-segment**

The incoming-PSF path must use **readpsf mode** (`load_project`, `guesscoord=False`,
`regenerate=False`), *not* build mode. Rationale:

- The whole point is "without rebuilding." Re-segmenting would re-derive topology from the PDB and
  **discard** the incoming PSF's patches, non-standard residues, custom bonds, and any topology
  pestifer can't reproduce — defeating the purpose.
- A foreign PSF may contain residues or patches pestifer has no `RESI`/knowledge of; build mode would
  hard-fail on them, readpsf mode carries them through.
- `merge` already proves readpsf-preserve works end-to-end.

Consequence: the incoming-PSF path is a **new, separate `psfgen()` code path** (call it
`psfgen_from_incoming()`), parallel to the segment-rebuild `psfgen()`, selected by the source mode.
Mods that *edit* topology are then expressed as psfgen operations layered on the loaded project
(patches, `regenerate` after patching), or as pure-numpy coordinate ops — never by re-segmenting the
whole system.

## Proposed design

### 1. Source specification (schema) — DECIDED

**Reuse the existing `prebuilt` source key** (rather than mint a new one) and **branch on
provenance**: `prebuilt` declared in the user's config → the new readpsf-preserve path; `prebuilt`
*injected at runtime* by the pipeline (basetask.py:335, psfgen.py:838) → the existing
state-continuation behavior, unchanged. Working shape:

```yaml
tasks:
  - psfgen:
      source:
        prebuilt:
          psf:  system.psf        # required
          pdb:  system.pdb        # coordinates as a .pdb ...
          # coor: system.coor     # ... or a NAMD binary restart .coor ...
          # xsc:  system.xsc      # ... paired with its .xsc cell
      mods:
        ...                        # the subset that is meaningful (see §4)
```

Provenance discriminator: the readpsf-preserve path fires when `prebuilt` is present in the **parsed
config specs** for this task; the continuation path is the runtime-injected case where
`ingest_molecules` adds `prebuilt` because there is no `base_coordinates` artifact but there is a
prior STATE. Concretely, gate on "was `prebuilt` in `self.specs['source']` before `ingest_molecules`
touched it" (e.g. a `self._user_prebuilt` flag captured at task provision), so the two never collide.

**Coordinate form — DECIDED: accept both.** `prebuilt.pdb`, *or* `prebuilt.coor` + `prebuilt.xsc`.
Since the `Molecule` constructor derives structure from a *PDB*, normalize the binary-restart case
up-front: when only `.coor`+`.xsc` are given, run pestifer's existing `coor_to_pdb(psf, coor)` to
materialize a `.pdb` (restoring chain IDs from nothing available yet — segid-derived is fine here),
then the ingest is uniform. The `.xsc` also seeds the initial cell for downstream NPT/NPgT.

The task's `pipeline_contract` for this mode: `requires=()` (no fetched SOURCE),
`provides=(STATE, MOLECULE)`, `discards_state=True` — it originates a system, like the fetch-build
psfgen, but from files on disk rather than a fetched ID. A `continuation`/`fetch` predecessor is
**not** required (today psfgen asserts one; that assertion is relaxed for this mode).

### 2. Ingest + classify

- Build the base `Molecule` from the incoming **PDB** for structure (residues/segments/coords) and
  attach the **PSF** as the authoritative topology reference. (Reuses the molecule.py `prebuilt`
  constructor logic; may share code but a distinct source key.)
- **Classify segtypes from resnames** via the existing `apply_segtypes()`. Add an **unknown-resname
  audit**: any residue whose resname is not in `segtype_of_resname` is reported up-front (grouped by
  resname + segid) rather than silently classified blank. Unknown *is allowed* (the PSF carries its
  topology), but the user is told, because segtype drives downstream selections (membrane embedding's
  protein selection, water/ion handling, etc.).
- **Segid inventory.** Read segids from the PSF (`PSFContents`), report the segid → (segtype, size,
  resname sample) table — the same courtesy `new-system --inspect` gives for a fetched structure.

### 3. Force-field consistency (fail up-front, not mid-NAMD) — DECIDED: hard-error

Before any downstream MD, verify the build's CHARMM release resolves **every atom type and every
bond/angle/dihedral/improper type-tuple** present in the incoming PSF against the loaded parameter
set. On a miss, **hard-error** with an actionable message naming the exact residue and atom-type term
(generalize the `make_solvent_box` missing-parameter reporter). No silent proceed, no auto-fill in
P1 — a CGenFF auto-fill path is deferred until the CGenFF-automagic roadmap item exists, and would
compose here later. This is the single most valuable guard: a CHARMM-GUI PSF built against a
*different* CHARMM version is the most likely real-world failure, and today it would surface as a
cryptic NAMD abort deep in equilibration. Target message:

```
ERROR: incoming PSF atom type 'CGX' (residue LIG 401) has no parameter in charmmff feb26.
Bond term 'CGX-CG2O1' also unresolved.
```

### 4. Which mods are meaningful (and error cleanly on the rest)

The incoming PSF has **no source metadata** (no `SEQADV`, `REMARK 465`, `BIOMT`). Partition the mod
pipeline explicitly:

- **Applicable (layer on via readpsf/patch/coords):**
  - `patches`, `ssbonds`, `links` → psfgen `patch` on existing segids, then `regenerate`.
  - `grafts` (add a glycan), ligand attachment → `coordpdb` the donor + patch + `regenerate`.
  - `crotations`/`irotations`, `orient` → pure-numpy coordinate ops (already numpy, post-load).
  - membrane embedding / solvation / equilibration → downstream tasks, not psfgen mods; they consume
    the STATE this task provides.
- **Meaningless for a foreign PSF → hard-error with a clear message:**
  - `biological_assembly` expansion (no BIOMT), `SEQADV`-derived mutations, `REMARK 465` missing-loop
    healing, `terminal_tails` building — all require fetch-source metadata that a PSF lacks.
- **Delicate → defer (out of scope for P1; error as "not yet supported"):**
  - `mutations`, `deletions`, `insertions` that require *re-segmenting* an affected chain. On a
    preserved topology these mean delete-residue + rebuild-residue psfgen surgery; correct handling is
    its own milestone. Fail cleanly rather than silently no-op.

The partition is enforced in code (a mode-specific allowlist), so the failure mode is a precise
up-front error, never a wrong build.

### 5. Segid / chain reconciliation

- **Preserve incoming segids verbatim** (readpsf keeps them). CHARMM-GUI conventions (`PROA`, `MEMB`,
  `WATA`, `IONS`, ...) flow through; classification is by resname, not segid, so this is fine.
- Reconcile with pestifer's **chain-ID persistence** machinery (`restore_chain_ids_from_reference`,
  keyed by `(segid, resid)`): the incoming PDB is the reference, so chain IDs survive downstream
  PSF→PDB regenerations. Report (don't silently remap) any segid that collides with pestifer's
  reserved conventions.

## Milestones

- **P1 — ingest + verify + pass-through (no topology edits).** New source key + readpsf ingest path;
  segtype/segid classification with unknown-resname audit; force-field-consistency preflight; provide
  a clean STATE. Acceptance: take a CHARMM-GUI (or pestifer-produced foreign) PSF/PDB straight into a
  `md`/`membrane_equilibrate`/`terminate` pipeline and run it end-to-end. This alone delivers the
  interop workflow and is the bulk of the value.
- **P2 — additive edits.** `patches`/`ssbonds`/`links`/`grafts` layered via readpsf + patch +
  `regenerate`; coordinate mods (orient/rotate). Acceptance: graft a glycan onto an incoming PSF.
- **P3 — subtractive/mutating edits.** `mutations`/`deletions`/`insertions` via per-chain
  re-segmentation surgery on the preserved topology. Its own design pass; likely reuses build-mode
  segment machinery for just the touched chain.

## Testing

- **Unit:** the mod-allowlist partition (each disallowed mod raises the precise error); the
  force-field-consistency checker (missing atom type / missing bond-type tuple → named error) on a
  crafted PSF; segtype/unknown-resname audit on a PSF with a non-standard resname.
- **Integration:** a small fixture PSF/PDB (e.g. BPTI built by pestifer, then re-ingested via the new
  path — must reproduce an equivalent STATE) run through `md` → `terminate`. A CHARMM-GUI-style
  fixture (foreign segids `PROA`/`WATA`/`IONS`) for the classification + FF-consistency path.
- **Determinism:** lean on the seeded-RNG substrate so a re-ingested build is reproducible.

## Decisions

1. **Source key — DECIDED:** reuse `prebuilt`, branch on provenance (user-config `prebuilt` →
   readpsf-preserve; runtime-injected `prebuilt` → existing state-continuation).
2. **Coordinate form — DECIDED:** accept `.pdb` *or* `.coor`+`.xsc`; normalize the binary-restart
   case to a `.pdb` via `coor_to_pdb` before ingest.
3. **FF-consistency strictness — DECIDED:** hard-error naming the residue/term; no auto-fill in P1
   (CGenFF compose deferred).
4. **Unknown resname policy — open, leaning warn-and-proceed:** readpsf carries the topology so the
   RESI is not needed, but segtype `''` weakens downstream selections — report each unknown resname
   (grouped by segid), proceed, and let the user register a segtype if a downstream task needs it.
