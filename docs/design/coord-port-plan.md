# Bringing coordinate transformations back into Python

pestifer offloads much of its coordinate geometry — rigid-body moves, superpositions,
symmetry (BIOMT) image generation, torsion rotations — to VMD/Tcl. The Cfusion work
(v3.10.0) showed that doing this geometry in pure numpy is clean, removes a
Python→VMD→disk→Python round trip, and sidesteps a whole class of column-corruption bugs.
This doc inventories what is still offloaded and lays out a phased plan to move as much as
possible into Python.

**North star:** make psfgen *topology-only*. VMD/psfgen does segments, patches, and
`regenerate`; pestifer authors **every coordinate** in Python and hands psfgen
already-correct PDBs to `coordpdb`.

---

## 1. What is offloaded today

### Group A — build-time transforms, coordinates already in a Python `Molecule` (highest value)
These run during a psfgen build while the `Molecule`/`asymmetric_unit` still exists in
Python, so the numbers we need are already on hand.

- **BIOMT / biological-assembly image generation** — `$sel move {tmat}` at
  `scripters/psfgen.py:299` (`write_polymer_stanza`) and `:471` (`write_generic_stanza`),
  bracketed by `backup_selection`/`restore_selection`. The 4×4 (`Transform.tmat`) and the
  atoms it acts on are both already in Python; VMD is doing a matrix multiply we could do in
  numpy. **Biggest footprint, lowest risk.**
- **Cfusion orientation** — the Tcl in `scripters/psfgen.py:631-694`
  (`write_cfusion_presegment`), added in v3.10.0. Base coords are in Python; the donor is a
  fresh VMD load.
- **Graft alignment** — `scripters/psfgen.py:600-627` (`write_graft_presegment`), a
  `measure fit` of the graft onto a target stub (a Cfusion sibling).
- **Loop-N seeding** — `caco_str` (`molecule/residue.py:850`), emitted at `psfgen.py:355`.
  *Already* computes the N position in numpy (`util/coord.positionN`); only the final `coord`
  assignment crosses into Tcl. This is the template for the whole effort.

### Group B — post-psfgen rigid-body / alignment ops (`scripters/vmd.py`)
Operate on the built psf/pdb, where coordinates live **only in VMD**. Driven by user
directive objects (`objs/rottrans.py`, `orient.py`, `align.py`, `crot.py` — all just emit Tcl).

- `RotTrans`: TRANS (`moveby`), ROT (rotate about COM), AXISANGLE (three-atom normal),
  ALIGN (carry one direction onto another) — `vmd.py:396-469`.
- `Orient`: principal-axis reorientation — `vmd.py:356-378` (uses the third-party `orient`
  Tcl package + `La` SVD).
- `Align`: `measure fit` onto a loaded reference — `vmd.py:490-525`.
- `TransferCoords`: optional `measure fit` pre-align + a direct `$mob set {x y z} [$donor get {x y z}]`
  coordinate copy — `vmd.py:527-568`.

### Group C — torsion / internal-coordinate rotations (`resources/tcl/pkg/pestifer/`)
A genuinely different problem: needs a bond graph + a clash metric, not just matrices.

- `crot.tcl`: `brot`, `Crot_*`, `SCrot_*` (dihedral rotations), `fold_alpha`/`new_alpha`,
  `glycan_pendant_rotate`→`rotate_pendant` (bond-graph BFS to find the moving set).
  Called from `scripters/vmd.py`, `scripters/psfgen.py`, `tasks/ringcheck.py:283`.
- `declash.tcl`: `declash_loop`, `declash_pendant` — Metropolis MC using `measure contacts`.

### Group D — analysis-only Tcl, NOT called from Python (low priority)
`autools.tcl` (`make_au_from_aligning`, `overlay_protomers` — `measure fit`), `multimer.tcl`
(opening angles, rotation matrices), `axes.tcl` (symmetry-axis vecmath), and principal-axis
via the third-party `orient.tcl` + `la.tcl` (a whole Golub–Reinsch SVD library used only to
diagonalize a 3×3 inertia tensor — replaceable by one `numpy.linalg.eigh`). Not in the build
path.

### Dead weight to delete
- The entire `resources/tcl/archive/` tree — unreferenced legacy (predecessors of the live
  code and of already-ported Python).
- A dangling `PestiferPierce → checkpierce.tcl` entry in `pkg/pestifer/pkgIndex.tcl`
  (ring-piercing was already ported to `psfutil/ring_resolve.py`; the file doesn't exist in
  the live tree).

---

## 2. What already exists in Python (build on this)

- **`util/coord.py`**: `build_tmat` (4×4 from R+t), `rotate_points_about_axis` (Rodrigues,
  vectorized, non-mutating), `positionN` (NeRF amide-N placement), `measure_dihedral`,
  `mic_shift`, `lawofcos`, `coorddf_from_pdb` (PDB→DataFrame read), `pdb_replace_coords`
  (coordinate-only byte-surgery write), `standardize_pdb_columns` (wide-resname column repair).
- **`psfutil/loop_ccd.py`** (explicitly VMD-free): `place_atom_nerf`, `optimal_ccd_angle`
  (closed-form single-axis least-squares fit), `end_rmsd` (RMSD), `set_dihedral`,
  `dihedral_deg`, `ccd_close`, Ramachandran sampling.
- **`molecule/transform.py`**: the `Transform` class (4×4 `tmat`), built from BIOMT records —
  but it only has `write_TcL()`; **no method applies it to an `AtomList`.**
- **`molecule/atom.py`**: `Atom` (scalar `x/y/z`), `overwrite_position(s)`, `pdb_line()`
  (per-record formatting; currently unused in live code), `AtomList.from_pdb/from_cif`.
- **Two complete VMD-free pipelines already ship**: loop CCD (`tasks/ligate.py`) and glycan
  ring de-piercing (`psfutil/ring_resolve.py`).
- numpy / scipy / pandas are all hard dependencies. The everyday `Atom` model is scalar
  (no `(N,3)` store); the numpy arrays appear only inside the offline kernels above.

---

## 3. The enabling primitives (Phase 0)

Small, pure-Python, tool-free (good for the CI core). They unblock Groups A and B.

- [x] **`Transform.apply(coords)`** — apply the 4×4 to an `(N,3)` array (non-mutating) or an
      `AtomList` (mutates positions in place). Landed in `molecule/transform.py`, backed by the
      pure `util/coord.apply_tmat`.
- [x] **Kabsch/SVD superposition** → `util/coord.kabsch` returns `(R, t, rmsd)`;
      `Transform.superpose(mobile, ref)` wraps it to return `(Transform, rmsd)`. Reflection
      correction via `sign(det(V·Uᵀ))`; verified atom-for-atom against VMD `measure fit`/`move`
      (agreement at the single-precision/PDB-3-decimal rounding floor, ~7e-4).
- [x] **Vectorized `coords` accessor** on `AtomList` — `coords` property (getter gathers scalar
      `Atom.x/y/z` into `(N,3)`; setter scatters back as plain floats). Model stays scalar.
- [x] **PDB writer** — `molecule/pdbwrite.py` (`write_atoms_pdb`) + `AtomList.write_pdb()`.
      Reconstructs pidibble record shims from pestifer `Atom`s and calls
      `pidibble.pdbwrite.PDBWriter` under the **CHARMM dialect** (pidibble ≥1.8.0). Verified:
      wide resnames (`BGLCNA`/`ANE5AC`) keep x/y/z pinned at cols 31-54 and segID at 73-76, and
      `parse(write(x,'charmm'),'charmm')` round-trips 294 glycan atoms with 0 mismatches. *Note:*
      pestifer's `Atom._adapt` still derives `segname` from `chainID` on read and does not yet
      ingest the segID column — reading segIDs back (to fully retire `standardize_pdb_columns` on
      the read side) is a small follow-up on the parse side.
      **Decision (below): offload the line formatting to pidibble's new CHARMM write dialect**
      rather than growing pestifer's own; keep `pdb_replace_coords` as a pestifer-local
      coordinate-only operation.

---

## 4. Phased plan

- [x] **Phase 0 — primitives** (above). All four landed with unit tests: the three numpy
      primitives (`apply_tmat`/`Transform.apply`, `kabsch`/`Transform.superpose`,
      `AtomList.coords`) with VMD parity checks, plus the pidibble-backed CHARMM PDB writer
      (`AtomList.write_pdb`). Bumped the pidibble pin to ≥1.8.0.
- [~] **Phase 1 — Group A** (coords already in Python). *Cfusion + BIOMT done; graft deferred to
      Phase 3 (entangled with `glycan_pendant_rotate` torsion pre-conditioning).*
      - [x] **Cfusion orientation** — `write_cfusion_presegment` is now pure Python: parse the
            donor via `AtomList.from_pdb`, select the fusion domain in serial/sequence order
            (so an in-chain CRO chromophore stays put — `from_pdb` otherwise tails all HETATM),
            renumber, orient with `util/coord.orient_peptide_fusion` (numpy port of the
            translate+Rodrigues-rotate), and write the segfile with `write_pdb(dialect='standard')`
            (standard columns match psfgen's reader; the donor's ≤4-char resnames don't need the
            CHARMM dialect). Validated atom-for-atom vs VMD on example 27: segfile 1771 atoms, 0
            identity mismatches, 0.001 Å max coord diff; full build PSF 4982 atoms/5033 bonds
            identical. Unit + `needs_tools` integration tests added.
      - [x] **BIOMT / assembly image generation** — `_author_subsegment_pdb` authors every
            resolved polymer subsegment and every generic (glycan/ion/water) segment directly
            from *copies* of the AU atoms: relabel segname/chain, apply the image `tmat` via
            `Transform.apply`, write standard-dialect PDB in serial order. Deleted the
            `atomselect ... [set chain] set resid [list ...] move {tmat} writepdb` path and both
            backup/restore brackets (copying the AU makes restore unnecessary). The atoms already
            carry the resid the segment expects, so no positional `set resid` is needed. Validated
            atom-for-atom vs VMD on 4zmj assembly 1 (3-fold: identity + 2 images, protein +
            glycan) from **both PDB and mmCIF**: 69 intermediate PDBs / 14568 atoms, 0 identity
            mismatches, 0.001 Å max diff; final PSF (30747 atoms/31134 bonds/56313 angles) and PDB
            identical. Added a `needs_tools` image-invariant test; the glycan-bond integration
            test still passes.
      - [ ] **Graft alignment** — deferred to Phase 3: the presegment's `measure fit` is a clean
            `Transform.superpose` swap, but it is preceded by `glycan_pendant_rotate` dihedral
            pre-conditioning (Group C), so it can't shed VMD until the torsion port lands.
- [~] **Phase 2 — Group B** (`vmd.py` post-psfgen ops). *transrot/align/transfer_coords done;
      orient deferred (see below).*
      - [x] **`transrot` / `align` / `transfer_coords`** — `molecule/coordmanip.py`
            (`CoordManipulator`) loads the built psf+pdb into an `AtomList` plus a per-atom mass
            array and bond graph from `PSFContents`, translates the restricted VMD atomselect
            subset actually used (`all`/`protein`/`backbone`/`name`/`chain`/`segid`/`resid N`/
            `resid N to M`/`fragment N`, joined by `and`) into a boolean mask, applies the
            transform in numpy (rigid move about mass-COM, three-atom AXISANGLE, minimal
            vector-ALIGN, Kabsch `measure fit`, coord copy), and writes the pdb back. The
            `ManipulateTask` routes these three objtypes to the Python engine; the bond-crossing
            disconnection guard and atom-count congruency guards are reproduced. Verified
            atom-for-atom vs VMD on the 6PTI fixture: ROT and transfer exact (0.0 Å), vector-ALIGN
            / AXISANGLE / Kabsch-align at the 0.001 Å rounding floor; all 12 selection strings
            match VMD's atom sets exactly. Added pure-Python unit tests + the existing
            `needs_tools` integration tests pass on the new path.
      - [ ] **`orient`** (principal-axis) — left on VMD: it is unused in every shipping config, and
            VMD's `mevsvd_br` returns eigenvectors in an unsorted, sign-quirked order that resists
            faithful `numpy.linalg.eigh` reproduction; a best-effort port would risk a silent
            behavior change. The membrane-embed "orient" already goes through `transrot` ALIGN
            (now ported). Revisit alongside Phase 4's `la.tcl`/`orient.tcl` retirement.
      - [ ] **`crotations` / `irotations`** (torsion) — Phase 3 (Group C), stay on VMD for now.
- [~] **Phase 3 — Group C** (torsion + declash). *irotations/crotations (`brot`/`fold_alpha`)
      done; pendant/SCrot/declash pending.*
      - [x] **`brot` (PHI/PSI/OMEGA/CHI1/CHI2) + `fold_alpha` (ALPHA)** — added to
            `CoordManipulator`: a VMD-matching per-atom `residue` index, the exact per-angle moving
            sets (residue-index ranges + main-chain-atom exclusions) and rotation axes, plus
            `fold_alpha` (delta-rotate each residue's φ/ψ/ω toward the α targets, reading current
            angles via `loop_ccd.dihedral_deg`). Both the psfgen task's post-build crotation
            coormods (iterating BA-image `chainIDmap`s) and the manipulate task now route
            irotations/crotations to the numpy engine. Verified atom-for-atom vs VMD: `brot`
            phi/psi exact (0.0 Å), `fold_alpha` at the 0.001 Å floor on 6PTI, and the full psfgen
            crotation path on 4tvp (env trimer, 3 images, crotations on B→B/D/F, a −180° loop
            flip moving atoms up to 147 Å) reproduced to 0.003 Å across all 33 693 atoms with 0
            identity mismatches. Pure-Python unit tests added.
      - [ ] **`rotate_pendant` / `glycan_pendant_rotate`** — bond-graph BFS pendant rotation;
            needed to unblock the deferred Phase 1 **graft** pre-conditioning. `ring_resolve.py`
            already has the numpy pendant machinery to reuse.
      - [ ] **`SCrot_chi1/chi2`** — ringcheck side-chain trial rotations (distinct from the `brot`
            CHI path).
      - [ ] **`declash_loop` / `declash_pendant`** — greedy contact-count declash on the psfgen
            build path (protein loops / NA loops / glycans). `loop_ccd` already supersedes
            `declash_loop` for the ligate path; `psfring.clash_count` (cKDTree) is the reusable
            contact counter for a numpy `declash_pendant`.
- [ ] **Phase 4 — cleanup**. Delete the `archive/` tree + dangling package entry; retire
      `la.tcl`/`orient.tcl` behind `numpy.linalg`; decide whether Group D utilities get thin
      Python ports or stay.

---

## 5. Key decisions & risks

1. **Numerical parity with VMD.** Kabsch handedness (reflection correction) and transform
   order must reproduce `measure fit`/`move` within tolerance. Mitigation: a regression
   harness comparing Python vs VMD on a few representative builds (assembly image, graft,
   align) — cheap since VMD is on hand.
2. **Where the "pipeline coordinate state" lives for Phase 2.** Options: keep a live
   in-memory `Molecule` across tasks, or load/transform/save per task (matches the current
   file-passing style). Load/save-per-task is lower-risk and incremental — start there.
3. **Bond-graph source for Phase 3.** Post-psfgen we have PSF connectivity; the crot
   "which atoms move" logic must come from that, not VMD `getbonds`.

## 6. PDB-writer decision — offload to pidibble

pestifer's `Atom` is already built *from* pidibble records, so a pestifer-own writer means
two independent encodings of PDB column layout that must agree forever — the exact drift that
caused the wide-resname corruption. Offloading the *structure-authoring* writer to pidibble
(the writer is in active development there) gives a single source of truth and lets us
**retire `standardize_pdb_columns`** entirely, because one library reading and writing with a
consistent column model can't desync with itself.

The condition that makes this safe: pidibble must offer a **CHARMM write dialect** that
round-trips faithfully — `parse(write(x, "charmm"), "charmm") == x` — including the segID
column (73–76), un-truncated ≥6-char resnames (`BGLCNA`, `ANE5AC`), altLoc, and iCode, with
the numeric block (x/y/z at 31–54, occ, beta, segID, element) pinned regardless of resName
width. (Full spec handed to the pidibble effort separately.) pestifer will reconstruct
records from its own `Atom`s and call the writer — no need to retain the parse.

**Keep in pestifer regardless:** `pdb_replace_coords` (copy a PDB verbatim, swap only cols
30–54). It's a different operation with a provable non-destructiveness guarantee, not a
"writer."
