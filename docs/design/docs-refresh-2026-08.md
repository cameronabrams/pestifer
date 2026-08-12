# Documentation refresh, August 2026

Checkpoint for the post-membrane-builder documentation update. Started 2026-08-08.

**Motivation:** `docs/source` was last touched 2026-07-23, while the membrane work and the
incoming-PSF feature ran through 2026-08-07 (releases v3.13 → v3.15). Several pages described
behavior that no longer existed.

## Status

- **Tier 1 (generated trees): DONE, committed** (`27f70da3`).
- **Tier 2 (prose corrections): DONE, committed** (`6e0e00a0`).
- **Tier 3 (new pages): DONE for `membrane_equilibrate`;** the optional "importing a pre-built
  system" narrative page is still open.
- **Tier 4 (examples + README): DONE.** ex17 and ex16 refreshed with figures from the clean
  2026-08-11 sweep run (including the `membrane_equilibrate` convergence plots, which did not
  exist when those pages were written) and prose brought onto the current membrane path; README
  "Key capabilities" gained the phase/self-terminating-equilibration detail plus new
  incoming-PSF and restartable-build bullets.

Verified after Tier 3: `make config-ref` is idempotent against the committed tree (only the footer
date stamp changes), and `sphinx-build` emits **zero** non-autodoc warnings — the `domainswap`
dangling reference that was the lone warning is gone with the page.

Note: the many untracked `tests/unit/test_tasks/**` files in `git status` are test artifacts that
predate this work; they are unrelated.

## What was found

### `make config-ref` was failing outright (fixed)

Ten schema nodes added during the membrane work had no `text:` key, and ycleptic's doc generator
raises `KeyError: 'text'` on those. The offending nodes were `constraints.atoms/k/consfile` and
`other_parameters.useflexiblecell/useconstantratio`, duplicated across the `patch` and `quilt`
entries of `make_membrane_system.bilayer.relaxation_protocols`.

Because `.readthedocs.yaml` regenerates the config reference in `pre_build`, **the RTD build had
been failing at that step since those keys landed.** Fixed in `pestifer/schema/base.yaml`.

### Corrections to the initial survey

Two things assumed in the first pass turned out to be wrong, and are recorded here so they are not
re-derived:

- `docs/source/api/` is **gitignored**. There are no committed API stubs to refresh; RTD regenerates
  them in `pre_build`. The stubs visible locally are untracked leftovers.
- RTD regenerates `config_ref` too, so the stale *committed* tree was never what got published. The
  staleness was repo hygiene, not a broken published reference — but the generator *crash* above was
  real and did break the published build.

### `domainswap` is retired but still documented

Regenerating deleted `docs/source/config_ref/tasks/domainswap.rst`, because commit `1302b0bc`
retired the task (code kept importable, removed from schema and registry). Left behind:

- `docs/source/subs/buildtasks/domainswap.rst` still exists, is **orphaned from every toctree**, and
  cross-references `config_ref tasks domainswap`, which no longer resolves — the one non-autodoc
  Sphinx warning in the build.

**Resolved 2026-08-09: deleted.** The page was orphaned from every toctree and its cross-reference to
`config_ref tasks domainswap` no longer resolved; git history keeps the prose, and the code stays
importable. This cleared the build's only non-autodoc warning.

## Work completed (uncommitted)

### Schema (`pestifer/schema/base.yaml`)

- Added the 10 missing `text:` entries described above.
- Rewrote the `packer` description, which was itself stale doc-bearing text: it still described
  solvent placed on a 3-D lattice, and did not mention the orthohexagonal lipid lattice.

### Regenerated `docs/source/config_ref/`

Picks up `membrane_equilibrate`, `density_equilibrate`, `upper/lower_leaflet_phase`,
`verify_parameters`, `quilt_grid_slack`, `lipids_per_leaflet`, `two_stage`, and the chunking and
area-convergence keys. Deletes the retired `domainswap` page.

Every file shows as modified because the generator stamps a date into each page footer; the
substantive diff is much smaller than the file count suggests.

### Prose

| File | Change |
| --- | --- |
| `subs/build.rst` | Removed the false claim that the controller rejects `continuation` → `psfgen` (that is now the incoming-PSF path). Added a "Resuming an interrupted build" section: `--restart`, `--from`, `--fresh`, run manifest. |
| `subs/buildtasks/make_membrane_system.rst` | Orthohexagonal lattice and why; MC-sampled conformers with per-lipid ensemble draw as the default; new leaflet-phase (Lo/Ld) subsection; `relaxation_protocols` now documents `membrane_equilibrate` stages and the `constant_ratio` calibration trap; new "Solvent chambers" section; `quilt_grid_slack`; replaced the stale example with the current ex17 build; updated Task Flow prose and mermaid diagram. Also fixed a pre-existing indentation bug rendering the `lipids`/`mole_fractions` bullet as a nested blockquote, and the malformed duplicate-key YAML in the old example. |
| `subs/buildtasks/continuation.rst` | New "Importing a foreign PSF" section documenting `verify_parameters`. |
| `subs/buildtasks/psfgen.rst` | New "Build mode vs. editing an incoming topology" section: additive vs. re-segmenting edits, the rebuild's RESI requirement. Corrected stale 2-digit task basenames (`00-00-00_` → `00-01-000_`). |
| `subs/make-pdb-collection.rst` | New "Choosing a conformer sampler" section: `--sampler`, cylinder sizing, sterols, pointer to the Lo ensembles. |

One claim was written and then removed as unverifiable: psfgen does **not** log which mode it took at
startup. It only errors, listing unsupported mods, when edit mode cannot apply a requested mod.

## Open decisions

Both resolved on 2026-08-09:

1. **Commit granularity** — split as suggested: the schema fix plus the regenerated `config_ref`
   landed as `27f70da3`, the prose as `6e0e00a0`.
2. **`domainswap` page** — deleted (see above).

## Remaining work

### Tier 3 — new pages

- ~~`subs/buildtasks/membrane_equilibrate.rst`~~ **Done.** Modeled on `density_equilibrate.rst`:
  why density *and* area is the minimal complete observable set for an anisotropic cell, the
  two-stage (const-area → tensionless) protocol and the co-compaction failure it avoids,
  `constant_ratio` on calibration patches, why the run is chunked, protein-corrected per-leaflet
  APL, the `area_*` tolerance family and the `area_plateau_tol` gate, and the report/plot outputs.
  Listed in `build.rst`'s toctree after `density_equilibrate`. Also fixed the now-stale closing
  sentence of `density_equilibrate.rst`, which told membrane users to keep a hand-written NPgT
  protocol; it now points at the sibling task.
- Still open (optional): a short "importing a pre-built system" narrative page tying `continuation`
  and `psfgen` together. The incoming-PSF feature is the one most likely to be missed by someone
  scanning task pages.

### Tier 4 — examples and front matter

- `examples/17/hiv-mpertm3-membrane2.rst` (last touched 2026-07-07) predates the entire membrane
  rework: hex packing, phase-aware calibration, solvate chambers, `membrane_equilibrate`.
- `examples/16/…rst` predates the solvate-chamber speedup (post-embed ~3h → ~1h49m).
- `README.md` "Key capabilities" (included verbatim into `index.md`): the membrane bullet has no
  phase or self-terminating-equilibration story, and there is no incoming-PSF or restart bullet.

**Figure regeneration: resolved.** The 2026-08-11 sweep ran all 27 examples clean, so ex16 and ex17
figures were harvested from it (ex17 8.8 h, ex16 5.2 h). Density profiles were regenerated from each
run's production system with `pestifer density-profile`. The one figure deliberately *not*
regenerated is `my_6e8w_viral.png` (and ex16's `my_6e8w_pc.png`): hand-composed VMD renders whose
captions document specific per-species colors, which the automatic snapshot renderer does not
reproduce.

The `membrane_equilibrate` plot titles were clipped at both ends (the task name plus description
overran the figure width); fixed at the source to a two-line, smaller title, and the already-rendered
figures were cropped for the docs.

## Reproducing the verification

Sphinx is not installed in `.venv`. A throwaway venv works:

```
uv venv docsenv --python 3.12
uv pip install --python docsenv/bin/python sphinx myst-parser sphinx-copybutton \
    sphinxcontrib-mermaid furo sphinx-rtd-theme alabaster pillow
uv pip install --python docsenv/bin/python -e . --no-deps   # conf.py needs pestifer metadata
docsenv/bin/sphinx-build -b html docs/source /tmp/docsbuild
```

`--no-deps` makes every autodoc import fail, so filter those out; what remains is the real warning
set:

```
docsenv/bin/sphinx-build -b html docs/source /tmp/docsbuild 2>&1 \
    | grep -E "WARNING|ERROR" | grep -viE "autodoc|failed to import|AttributeError"
```

As of 2026-08-09 that filter is **empty** — the `domainswap` warning was the last one.
