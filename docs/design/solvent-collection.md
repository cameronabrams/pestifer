# Design: a `solvent` PDB-repository collection & non-water solvation

Status: **draft for review** (no code yet)

## Problem

Pestifer can only solvate with **TIP3P water**. The `solvate` task and
`bilayer_embed.tcl` both call VMD's `solvate` plugin with just `-minmax`, so they
tile the plugin's built-in pre-equilibrated water box. To bulk-solvate with any
other solvent, VMD `solvate` requires you to hand it a **pre-equilibrated periodic
box** of that solvent as `-spsf`/`-spdb` plus its exact edge length (`-ws`). We
have nowhere to store such boxes and no way to build them.

## Two consumers, two artifact kinds (the crux)

"Solvent coordinates" are consumed two different ways, needing two different
artifacts:

| Consumer | Needs | Today |
|---|---|---|
| **Grid membrane packer** (`bilayer.py`) | a **single molecule** to place on its 3-D chamber lattice; already species-agnostic via `solvent_specstring`/mole-fractions | `TIP3` + ions, as single-molecule conformer PDBs in `water_ions` |
| **VMD `solvate`** (`solvate` task, embed slabs) | a **pre-equilibrated periodic box** (`-spsf/-spdb/-ws`) | only the plugin's built-in TIP3 box |

So a `solvent` collection must hold entries of two **kinds**, made explicit in
`info.yaml` so a consumer never grabs the wrong one.

Leverage: because the grid packer is already species-agnostic, a *single-molecule*
solvent entry unlocks that solvent for **grid membrane builds immediately**; the
**box** is specifically the unlock for the **VMD-solvate** path. A fully-supported
solvent generally wants both (and the single molecule can be extracted from the
equilibrated box).

## Spike findings (step 3 — DONE, the biggest unknown is resolved)

Ran VMD 2.0.0's `solvate` plugin directly (read `solvate1.7/solvate.tcl` + empirical
round-trips) to characterize the custom-box path. Results:

- **The `-spsf/-spdb/-ws/-ks` custom-box path works and is solvent-agnostic.** Given
  the built-in water box explicitly, it reproduces the default solvation exactly
  (4,454 waters, identical). It simply tiles whatever periodic box you hand it.
- **`solvate` does NOT clean solvent–solvent seam overlaps.** Its only trimming is
  "delete solvent within `-b`/rwat of the *solute*, or outside the target box." So the
  box must be periodic-clean by construction, and **`-ws` must be the box's *exact*
  edge**: tiling the clean 65.42 Å box at the correct `-ws` over a 3×3×3 region gave
  **0** O–O clashes < 2.0 Å; tiling the same box at `-ws 60.0` (a 5 Å seam overlap)
  gave **29,011** clashes with no cleanup.
- **`-ks` (key atom) is required for a non-water box and must name a one-per-residue
  atom.** With a valid non-`OH2` key it works identically; with a key that matches
  nothing, the trim selection matches nothing and *no* solvent is deleted (neither
  solute-overlapping nor out-of-box) → a broken, clashing system. The default key is
  `name OH2`, so non-water solvents MUST pass `-ks`.

**Implications for the builder (step 4):** a `box` entry must ship (a) a genuinely
PBC-equilibrated periodic cube (clean at its faces), (b) its **exact** equilibrated
edge (→ `-ws`; small errors are catastrophic, no cleanup), and (c) a **key atom name**
that occurs once per molecule (→ `-ks`). These become required `info.yaml` fields.

## Locked decisions

- **The collection holds BOTH kinds, permanently — it is *not* boxes-only.** The
  single-molecule (`kind: molecule`) entries stay: the **grid membrane packer's
  chamber fill** checks out the water molecule (`checkout('TIP3')`) and ion molecules
  and lattice-places them (`bilayer.py`/`write_grid_pdb`; `salt_concentration` just
  sets how many ions to drop). They were *not* only for the removed packmol packer. So
  `solvent` = molecule entries (grid chamber) **+** box entries (VMD `solvate`).
- **Grid chamber-fill keeps single-molecule lattice placement** (it's validated and
  fast). We are *not* migrating it to `solvate`/`autoionize`; therefore the molecule
  entries are load-bearing and are not deleted.
- **Single-species boxes first.** A box is one solvent RESI. Cosolvent *mixtures*
  (e.g. 80/20 water/glycerol) are deferred.
- **Ions stay under the `solvent` umbrella.** `solvent` = "chamber-fill species:
  water + ions + cosolvents" — matches how the grid packer already treats chamber
  species. No separate `ions` collection.
- **Water keeps using VMD's built-in box.** Boxes are used only for *non-water*
  solvents; the default water path is untouched (zero regression risk).

## Design

### 1. Collection rename: `water_ions` → `solvent`

- `resources/charmmff/<ver>/pdbrepository/water_ions.tgz` → `solvent.tgz`, and its
  single top-level directory `water_ions/` → `solvent/` (the collection streamID is
  the tarball basename, so both must change together).
- Update consumers that name it: the segtype→collection default in
  `ResourceManager.add_pdb_entry` (`ion`/`water` → `solvent`), any doc/text
  mentions, and the `show-resources`/`make-pdb-collection` help. `water_ions` is a
  pestifer-coined name (there is no CHARMM `water_ions` *stream*), so nothing in the
  force field depends on it.
- Internal-only change; no user config references the collection by name.

### 2. `info.yaml` schema: a `kind` discriminator

Add `kind: molecule | box` (a loader that sees no `kind` treats it as `molecule`,
so every existing entry is backward-compatible).

- `kind: molecule` (existing shape): `charge`, `conformers: [{pdb, head-tail-length,
  max-internal-length}]`, `reference-atoms`, `defined-in`, `parameters`, `synonym`.
- `kind: box` (new): a single pre-equilibrated periodic cube.
  ```yaml
  kind: box
  resname: GROL
  nmol: 216
  density: 1.26          # g/cc, measured post-equilibration
  box_edge: 24.13        # Å, the EXACT equilibrated cubic edge -> VMD solvate -ws
  key_atom: C1           # an atom name occurring once per molecule -> VMD solvate -ks
  psf: GROL-box.psf      # -> VMD solvate -spsf
  pdb: GROL-box.pdb      # -> VMD solvate -spdb
  parameters: [...]      # param files required to read/relax the box
  ```
  (`box_edge` and `key_atom` are load-bearing per the spike: a wrong edge leaves
  uncleaned seam clashes, and a wrong/absent key atom disables all trimming.)
  `PDBCollection.build_from_resources` learns to load a `box` entry (psf+pdb + the
  scalar box metadata) instead of a conformer list.

### 3. `make-pdb-collection`: a `solvent` box-build mode

A distinct pipeline from the current per-RESI conformer sampling (which pulls all
RESIs from a CHARMM *stream*). "Solvent" is a **category**, not a stream: a curated
solvent named explicitly with box parameters.

Sketch: `pestifer make-pdb-collection solvent --resname GROL --density 1.26 --nmol 216`

Pipeline:
1. Build one solvent molecule (psfgen from the RESI).
2. Pack `nmol` copies into a cube sized to the target density (grid/random seed).
3. **Minimize + NPT-equilibrate under PBC** at target T/P so the box relaxes to its
   equilibrium density and edge.
4. Wrap into the primary cell; **verify no cross-boundary (seam) clashes** — the box
   must tile seamlessly, the same quality as VMD's built-in watbox.
5. Record the equilibrated **edge length** and density; write the box psf/pdb and a
   `kind: box` `info.yaml`.

Output is an entry directory `<RESN>/` installable with the machinery we already
have: `modify-package pdb-repo add-entry <dir> --collection solvent`.

Uses NAMD (like the existing conformer sampling), so it's heavier than a molecule
entry.

### 4. `solvate` task integration

- Add a `solvent` spec (resname; default `TIP3`).
- `TIP3`/default → current behavior (VMD `solvate` with its built-in box).
- Non-water → look up the `kind: box` entry in the `solvent` collection, check out
  its psf/pdb, and call `solvate … -spsf <box.psf> -spdb <box.pdb> -ws <box_edge>
  -ks "name <key_atom>" -minmax …`. Autoionize (ions) still runs afterward.
- `bilayer_embed.tcl` non-water slabs are a later extension; initial scope is the
  standalone `solvate` task.

### 5. How the pieces connect (little new plumbing)

- `modify-package pdb-repo add-entry` already installs entries into a collection
  tarball; once the schema supports `kind: box` it installs boxes with no new code,
  and the segtype→collection default routes solvent species to `solvent`.
- Grid packer: single-molecule `solvent` entries are picked up by the existing
  species-agnostic chamber logic.

## Open questions / risks

- ~~**VMD `solvate` `-ws`/`-spdb`/`-spsf` exact semantics.**~~ **RESOLVED by the spike**
  (see "Spike findings" above): the path works, is solvent-agnostic, does *not* clean
  seam overlaps, and needs the exact edge (`-ws`) plus a one-per-residue key atom
  (`-ks`).
- **Box quality is the remaining real work.** Because `solvate` won't fix seams, the
  step-4 builder's PBC equilibration + wrap must yield a box that is clean at its faces
  and whose recorded edge is exactly the tiling period. This is where the effort goes.
- **Edge vs density trade-off.** An NPT box fluctuates; we tile at a single recorded
  edge, so density has a small error. Record the mean edge; quantify the error. (A
  final short NVT at the target edge, or averaging the edge over the NPT tail, keeps
  the tiled density on-target.)
- **Do we ever make water a `box` entry?** No — keep VMD's default water path;
  boxes are non-water only. (Revisit only if we want reproducible water boxes.)

## Staged implementation plan

1. ~~**Rename** `water_ions` → `solvent` (collection + inner dir + consumers + docs).~~
   **DONE.** Repacked `water_ions.tgz` → `solvent.tgz` (inner `water_ions/` → `solvent/`)
   for both releases; updated the `ion`/`water` → `solvent` default in `add_pdb_entry`,
   docstrings, docs, and tests. Left the CHARMM *stream* label `water_ions`
   (`CHARMMFFStreamID`, `segtype_classifier`, `toppar_water_ions.str`) untouched — the
   pestifer collection and the CHARMM stream intentionally diverge now.
2. ~~**Schema**: add `kind:` to `info.yaml` + loader (default `molecule`).~~ **DONE.**
   `PDBCollection.build_from_resources` branches on `kind` (box → load psf+pdb+metadata;
   molecule → conformers as before; missing `kind` = molecule), and `PDBInput` exposes
   `is_box`/`get_box_edge`/`get_key_atom`/`get_box_psf`/`get_box_pdb`. Backward-compatible.
3. ~~**Spike**: prove the VMD-`solvate` `-spsf/-spdb/-ws` round-trip.~~ **DONE** — see
   "Spike findings"; the mechanism works, and the box requirements (periodic-clean,
   exact `-ws` edge, one-per-residue `-ks` key atom) are pinned down.
4. **`make-pdb-collection solvent`**: the pack + NPT-equilibrate + seam-check + edge
   pipeline, for one test solvent.
5. **`solvate` task**: the non-water box path.
6. **Docs** + a worked example; extend to `bilayer_embed.tcl` slabs if desired.
