# Design: offload all mmCIF parsing to pidibble

Status: **implemented.** Roadmap item "Offload all mmCIF parsing to pidibble" (Tooling / packaging).
mmCIF is now parsed by pidibble (`PDBParser(input_format='mmCIF')`); `util/cifutil.py` and the direct
`mmcif` / `mmcif-pdbx` dependencies are gone. Label-primary chain identity preserved; the golden
oracle diffs to only the intentional symmetry-operator normalization. See the phase log at the end.

## Problem

pestifer carries **two structure parsers**. PDB input goes through pidibble
(`PDBParser(...).parse().parsed` → a `PDBRecordDict`); mmCIF input goes through the `mmcif` /
`mmcif-pdbx` libraries (`from mmcif.api.PdbxContainers import DataContainer`,
`from mmcif.io.IoAdapterCore import IoAdapterCore`), wrapped by `util/cifutil.py`'s `CIFload` +
`CIFdict`. Every structure object then carries **divergent** `from_pdb` and `from_cif` paths
(`molecule/{atom,residue,asymmetricunit,bioassemb}.py`, `objs/{ssbond,link,seqadv}.py`).

Two parsers means two dependencies (`mmcif`, `mmcif-pdbx>=1.1.2`), two field-mapping code paths per
object to keep in sync, and a type-based dispatch (`isinstance(parsed, DataContainer)`) threaded
through `Molecule`/`AsymmetricUnit`/`BioAssembList`. The goal: parse mmCIF with pidibble too, unify on
one parser, and drop the `mmcif` / `mmcif-pdbx` dependencies.

## What changed: pidibble 1.7.0 (the enabler)

The roadmap's original framing assumed we'd write *new* pidibble-based `from_cif` methods that read raw
mmCIF categories (`atom_site`, `struct_conn`, …). **pidibble 1.7.0 makes that obsolete.** Its mmCIF
parser (`PDBParser(input_format='mmCIF')`) normalizes mmCIF into the **same PDB-record namespace** as
the PDB path. Parsing 4zmj as mmCIF yields these `PDBRecordDict` keys:

```
ATOM, HETATM, LINK, SSBOND, SEQADV, SEQRES, HELIX, SHEET, COMPND, SOURCE,
REMARK.465, REMARK.350.BIOMOLECULE1.TRANSFORM{1,2,3}, HEADER, TITLE, EXPDTA, KEYWDS, CRYST1
```

So the migration is not "map raw CIF categories onto pestifer objects" — it is **"delete pestifer's
raw-CIF layer entirely and feed pidibble's normalized records into the record-adapter machinery pestifer
already has."** pidibble deliberately populates these records to mirror PDB, and — crucially — retains
**both** identities on every residue-bearing record:

- `.residue` / `.residue1` / `.initRes` … = **label** numbering (`label_asym_id`, `label_seq_id`)
- `.residue_auth` / `.residue1_auth` … = **author** numbering (`auth_asym_id`, `auth_seq_id`,
  `pdbx_PDB_ins_code`)

On the **PDB** path there is only `.residue`, and it already holds author numbering. This dual-identity
availability is what makes the migration behavior-preserving (next section).

pidibble 1.7.0 is on PyPI, so pestifer can depend on it and still cut releases.

## The pivotal decision: keep label-primary chain identity

pestifer's current CIF path is **label-primary**: it sets the primary `chainID` from `label_asym_id`
and stashes `auth_*` as metadata (`molecule/atom.py` `_adapt` CIFdict branch;
`residue.py:map_chainIDs_label_to_auth`, `residue.py:cif_residue_map`). Its PDB path is **auth-primary**
(`chainID` = author chain). This asymmetry is **deliberate and load-bearing**, not an accident:

- CIF label chains cleanly separate entities. In 4zmj every N-linked glycan branch is its own
  `label_asym_id` (N, O, P, …), but they all share `auth_asym_id` **G** with the protein they decorate.
  Label-primary gives pestifer one clean segment per entity directly from the file; auth-primary would
  **collapse glycans and ligands onto their parent protein chain**, breaking segtype-based
  segmentation, glycan grafting, and ion/water separation for CIF input.
- pestifer preserves the author numbering (the numbering users know from the paper/PDB) as metadata and
  emits a user-facing translation artifact, `cif_residue_map` (`residue.py:1029`, written in
  `tasks/psfgen.py:697` when `cif_residue_map_file` is requested):
  `label_chain:label_resid → auth_chain:auth_resid+icode`.

**Decision: preserve label-primary semantics.** pidibble hands us both identities per record, so we set
the primary `chainID`/`resid` from `.residue` (label) and stash `auth_*` from `.residue_auth`, exactly
as today. The entire downstream label→auth machinery (`map_chainIDs_label_to_auth`, `cif_residue_map`,
seqadv's label/auth `assign_residues` stash) stays **untouched**.

**Subtlety to preserve (found in Phase-0 baseline):** `cif_residue_map` is emitted for *both* formats,
not CIF-only — `Residue.auth_asym_id` is always a field, so the `hasattr` guard is always true. The
difference is the *value*: PDB input leaves `auth_asym_id=None`, so the map reads `A:1 → None:NoneNone`;
CIF input fills it, so the map reads `A:1 → A:1` (6pti, label==auth) or `A:1 → G:31` (4zmj,
label≠auth). The migration must keep exactly this split — PDB → `None` auth values, CIF → real author
numbering from `.residue_auth`. (So dispatch must be by explicit `source_format`, not by sniffing
`auth_asym_id` presence — which is always present.)

Rejected: **auth-primary unification** (route CIF through `from_pdb`, making CIF and PDB of the same
entry produce identical AsymmetricUnits, deleting the label→auth machinery). It is the more elegant end
state on paper, but it destroys the per-entity chain separation that pestifer relies on for glycans and
ligands. Not worth it.

## Field-parity findings (what reuses cleanly vs. needs a rewrite)

Empirically compared, per record, CIF-parsed vs PDB-parsed (4zmj):

- **Loop records — ATOM/HETATM, SSBOND, LINK, SEQADV** and **missing residues (REMARK.465)**: the CIF
  records carry the same core field names as their PDB counterparts and simply **add**
  `residue_auth`/`auth_*`. Each `from_cif` becomes a straightforward walk of pidibble's
  `PDBRecordList` reading `.residue*` (primary, label) + `.residue*_auth` (auth stash) — a direct
  swap for the old `CIFdict[...]` raw-key access.
- **Bioassembly (REMARK.350)** is the exception: CIF and PDB record **shapes differ**. CIF transform
  records expose `row1/row2/row3` (each `{m1,m2,m3,t}`) plus scalar `BIOMOLECULE`/`TRANSFORM`; PDB
  records expose a `row` list plus `tokens` metadata. pidibble ships `pdbparse.get_symm_ops(rec)` that
  returns the `(M, T)` 3×3/3-vector from **either** shape. So bioassembly gets a small rewrite onto
  that helper, not a free reuse of the PDB parser.

Confirmed non-issues (pestifer reads none of these from CIF today, so nothing regresses): missing
*atoms* (`pdbx_unobs_or_zero_occ_atoms`), DBREF, MODRES, SITE, ANISOU, CISPEP, CONECT — all unmapped on
pidibble's CIF path, and all already unused by pestifer's CIF path.

## Design

Two-parser removal, PDB path left untouched (it is the dominant path — no reason to risk it). Each CIF
adapter is repointed from `CIFdict`/`DataContainer` onto pidibble's normalized records; the raw-CIF
layer and the two dependencies are deleted.

### 1. Entry point & dispatch

- `molecule/molecule.py` (~lines 100–114): replace the `CIFload(...)` branch with
  `PDBParser(input_format='mmCIF', ...).parse().parsed`. Both branches now return a `PDBRecordDict`.
- Dispatch can no longer key on `isinstance(parsed, DataContainer)`. Thread the already-known
  `source['file_format']` from `Molecule` down into `AsymmetricUnit.__init__`
  (`asymmetricunit.py:131`) and `BioAssembList.__init__` (`bioassemb.py:148`) as an explicit
  `source_format` argument. (Preferred over sniffing `residue_auth` presence — explicit and readable.)

### 2. Loop-record adapters (5)

Rewrite each `from_cif` to iterate the pidibble `PDBRecordList`, reading label fields as primary and
`*_auth` as the auth stash, dropping `CIFdict` entirely:

| Object | File | pidibble record(s) | notes |
|---|---|---|---|
| Atom | `molecule/atom.py:343` | `ATOM` + `HETATM` | CIF now splits HETATM into its own key; iterate both |
| ResiduePlaceholder | `molecule/residue.py:298` | `REMARK.465.tables['MISSING']` | rows already carry `auth_chainID/auth_resName/auth_seqNum` |
| SSBond | `objs/ssbond.py:225` | `SSBOND` | normalize `serNum` `'disulf1'` (str) → int |
| Link | `objs/link.py:476` | `LINK` | pidibble already folds covale+metalc into `LINK` (matches old filter) |
| Seqadv | `objs/seqadv.py:280` | `SEQADV` | `struct_ref_seq_dif` → `SEQADV`; author strand id + `pdbx_auth_seq_num` semantics unchanged |

### 3. Bioassembly

Replace `BioAssembList.from_data_container` (`bioassemb.py:228`, the one place that used
`getObj`/`selectIndices`/matrix-element access directly) with a reader over the pidibble
`REMARK.350.BIOMOLECULE*.TRANSFORM*` records, using `pidibble.pdbparse.get_symm_ops` for `(M, T)`.

**Chain-list caveat (label vs auth — found running Example 8).** A transform's
`applies_chainIDs` must be in the *label* namespace to match the label-primary AU. pidibble's
`REMARK.350.header` is **author**-mapped (correct for PDB; lossy for label-primary mmCIF — it
collapses many label chains onto fewer author chains and drops non-polymer ones). pidibble 1.7.1
adds `header_label` (the raw label `asym_id_list`); `Transform._adapt` (`molecule/transform.py`)
prefers `header_label` when present, else `header`. Without this, assembly copies leave
waters/ions/glycan label-chains unremapped → psfgen `duplicate segment key` fatals. This is why
the migration requires `pidibble>=1.7.1`.

### 4. Delete the dead layer

- Delete `util/cifutil.py` (`CIFdict` + `CIFload`) and every `from mmcif…` / `DataContainer` import.
- Drop `mmcif` and `mmcif-pdbx` from `pyproject.toml` dependencies. (`mmcif-pdbx` provides the `mmcif`
  import namespace; both go.)

## Locked decisions

- **Label-primary is preserved** (see above). The `cif_residue_map` artifact and behavior are unchanged.
- **PDB adapters are not touched** in this migration. Merging each `from_pdb`/`from_cif` pair into a
  single adapter (now possible, since both consume pidibble records: always read `.residue` as primary
  and stash `.residue_auth` when present) is an **optional follow-up refactor**, deliberately deferred to
  keep the dominant PDB path out of the blast radius.
- **pidibble requirement** bumps to `>=1.7.0`.

## Open questions / risks

1. **VMD `.cif` strip. RESOLVED (found running Example 8): pestifer DOES hand the raw `.cif` to VMD.**
   The psfgen task loads the source structure with `mol new <id>.cif`, and VMD's pdbx plugin
   mis-parses a `pdbx_audit_revision_item` loop ("coordinate fields not found" — it reads only the
   audit rows as atoms). `CIFload` used to strip that block from the on-disk file as a side effect;
   the migration removed it, so VMD failed. Restored as `util.util.strip_cif_category` — a
   dependency-free text strip (no `mmcif` reimport) called on the on-disk `.cif` right after pidibble
   parses it (`molecule.py`). Verified: raw → VMD reads 41 atoms and errors; stripped → 12188 atoms.
2. **Fetch/caching location.** `PDBParser(input_format='mmCIF')` downloads `<id>.cif` into CWD. Confirm
   it lands where the artifact machinery expects (the `FetchTask` mmCIF path,
   `test_tasks/test_fetch.py:39`, produces a `CIFFileArtifact`).
3. **`HEADER.depDate` format.** ISO (`YYYY-MM-DD`) on pidibble's CIF path vs `DD-MON-YY` on PDB. Only
   matters if anything consumes depDate from CIF (believed nothing — verify with a grep).
4. **Regression surface for tests** that assert `DataContainer` types or call `CIFload`/`CIFdict`
   directly (`test_util/test_cif.py`, the `*_from_cif` unit tests). These rewrite to the pidibble
   record shape; the PDB↔CIF equivalence test (`test_asymmetricunit.py:30`, on 6pti) is the key gate.

## Staged plan

- **Phase 0 — pin + baseline. DONE.** Bumped `pidibble>=1.7.0`. Captured golden snapshots on the
  *current* (mmcif-backed) code to `~/devtests/mmcif-offload-golden/` (reproducer: `capture.py` there):
  AsymmetricUnit invariants — segments (chainID/segtype/counts), residue+atom totals, ssbonds, links,
  missing residues, and the full `cif_residue_map`. This is the equivalence oracle for every later
  phase. Baseline values:
  - **4zmj (mmCIF, glycosylated):** 20 segments (2 protein A/B + 18 glycan, each glycan tree its own
    segment — the label-primary entity separation), 659 residues, 4856 atoms, 11 ssbonds, 25 links,
    61 missing residues, `cif_residue_map` of 659 (label→auth, e.g. `A:1 → G:31`).
  - **6pti (mmCIF vs PDB equivalence pair):** both give 3 segments (protein A / ion B=PO4 / water C),
    132 residues, 536 atoms, 3 ssbonds (`A:5-A:55`, `A:14-A:38`, `A:30-A:51`), 0 links, 1 missing.
  - Full CIF-touching unit suite green (70 tests).
- **Phase 1 — dispatch seam. DONE.** Threaded an explicit `source_format` into
  `AsymmetricUnit`/`BioAssembList`, dispatch prefers it (isinstance fallback retained). Kept
  `CIFload`/`DataContainer` in place — zero behavior change, golden byte-identical. (Split from the
  doc's original Phase 1: the entry-point swap moved to Phase 2 since it's coupled to the adapters.)
- **Phase 2+3 — entry point, adapters, bioassembly. DONE (landed together — coupled).** `molecule.py`
  now parses mmCIF with `PDBParser(input_format='mmCIF')`; all five loop-record `from_cif` methods
  (Atom, ResiduePlaceholder, SSBond, Link, Seqadv) consume pidibble records (SSBond/Link/Seqadv/
  ResiduePlaceholder just delegate to `from_pdb`, since pidibble emits the same record keys;
  `Atom.from_cif` mirrors `from_pdb` with an ATOM+HETATM split sorted by serial). `_adapt`/`_from_pdbrecord`
  detect an mmCIF record by the presence of `residue*_auth` and read label numbering as primary +
  author as the stash. Bioassembly: CIF now routes through `from_pdb_record_dict` (pidibble's
  `get_symm_ops` handles the CIF `row1/row2/row3` shape); `from_data_container` deleted. The
  `DataContainer` isinstance fallback is gone. CIF unit tests rewired to the pidibble API. Result:
  golden diff is **only** the two intentional deltas below; full unit suite green.

  Two intentional deltas (behavior changes, both benign and documented):
  1. **Symmetry operators normalize `1_555` → `1555`.** pidibble rectifies mmCIF symmetry operators to
     numbers; `symmetry_str` (in `util/util.py`) coerces them to the PDB string form. This *aligns* CIF
     with the PDB path (which already yielded `1555`) — the two formats previously disagreed. `sym1`/
     `sym2` are only ever emitted in `pdb_line`; no functional use.
  2. **Type coercions for pidibble's numeric rectification.** `charge` (mmCIF `pdbx_formal_charge`, e.g.
     `1`) and purely-numeric chain ids (e.g. auth chain `'0'` on an 8fae glycan, rectified to `int 0`)
     are coerced back to `str` in the adapters — the old CIFdict path returned everything as strings.
- **Phase 4 — delete. DONE.** Deleted `util/cifutil.py`; removed the `mmcif` and `mmcif-pdbx`
  declarations from `pyproject.toml`. Both packages remain available transitively (`mmcif` via
  pidibble, which uses it for raw tokenizing; `mmcif-pdbx` via pdb2pqr), so nothing at runtime breaks —
  pestifer simply no longer declares or imports either. Re-resolved the env and re-ran the CIF suite +
  golden: green, only-sym deltas.
- **Phase 5 — final verify. DONE.** Full unit suite: 729 passed / 21 skipped; the only failures are 6
  pre-existing `test_pdbrepository.py` env failures that fail identically on the original code. Golden
  diff = only the symmetry-operator normalization.

## Bonus now available (out of scope; roadmap follow-up)

pidibble 1.7.0 also maps CIF **SEQRES** (`pdbx_poly_seq_scheme`), **HELIX/SHEET**, and
**COMPND/SOURCE** — none of which pestifer reads from CIF today. A later item could let pestifer use a
true reference sequence instead of reconstructing it from resolved + missing residues. Note COMPND/SOURCE
arrive as **flat per-entity records**, not PDB `tokengroups`, so any future consumer must branch on
format.
