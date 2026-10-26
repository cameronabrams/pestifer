# Design: build the membrane by tiling relaxed patches (patch quilting, revived)

Status: **proposed, 2026-10-26**. Not implemented. Approved in outline: tile each leaflet from its
own relaxed patch and **crop** to the target box (option (a) below), 100 lipids per leaflet per
patch by default.

## The decision in one paragraph

`make_membrane_system` relaxes small calibration patches, reads a number off them (the area per
lipid), throws them away, and then places the whole membrane from scratch on a grid at that number.
The thing it throws away is the only part of the build that is physically packed. The thing it keeps
is a scalar that, for ordered or cholesterol-rich leaflets, cannot be measured in a build's time
budget. This design keeps the patches instead: each leaflet of the full membrane is tiled from the
relaxed leaflet of its own patch and cropped to the box. Nothing is placed, so nothing is threaded,
and the full membrane starts at the patches' packing instead of condensing ~17% from a grid.

## Why: what the grid path costs

Every failure below is one consequence of gridding a full membrane at a measured APL.

**The APL is measured too early, and replicas disagree only about when they stopped.** Example 17
in the 3.29.8 sweep. All three replicas start patchA (upper leaflet: PSM/POPC/CHL1, 47% CHL1, Lo) at
APL 50.0 Å². Their paths through the tensionless stage are the same; their stopping points are not:

| step | rep-01 | rep-02 | rep-03 |
| ---: | ---: | ---: | ---: |
| ~176-181k | 47.74 | -- | 47.78 |
| 196k | 47.67 | -- | 47.61 |
| 166k | -- | **47.67, stopped CONVERGED** | -- |
| 216k | **47.48, stopped CONVERGED** | -- | 46.93 |
| 1.0M | | | 43.59 |
| 2.5M | | | **41.86, stopped at the step ceiling** |

rep-01 and rep-02 were declared converged after 35k and 40k steps of the tensionless stage. rep-03
kept going, contracted another 6 Å², and was refused by the calibration gate, which made it the
only replica correctly judged. The two that shipped were gridded with an upper leaflet ~14% too
loose.

**The plateau cannot be extrapolated either.** Fitting a single exponential to rep-03's trace using
only its first part:

| stage-2 steps fitted | fitted asymptote (Å²) |
| ---: | ---: |
| 200k | 44.38 ± 0.16 |
| 500k | 44.26 ± 0.02 |
| 1.0M | 43.29 ± 0.03 |
| 2.33M (all) | 41.61 ± 0.03 |

The estimate falls as the window grows while its error bar stays small. The relaxation has more
than one timescale, so each fit is confidently wrong. 41.6 Å² is an upper bound, not the APL.
No stopping rule tuned on a build budget can fix this. Data: the `[stage 2]` chunk lines of
`make_membrane_system-membrane_equilibrate-patchA` in each replica's `build.log` from the 3.29.8
sweep.

**Gridding at a condensed APL threads rings; gridding loose costs a long condensation.** An
asymmetric PIP2 build (another group, 3.29.7) pinned the upper leaflet at 41.58 Å², gridded at
43.7 Å² (the default 5% slack), and its first check found 108 piercings, 254 after rotation. Its
quilt minimize died in rattle2. Placing loose avoids that but hands the quilt the ~17% lateral
condensation that the patch already went through, on a system 3-7x larger (ex17's quilt is ~580
lipids, the PIP2 quilt 1,460, a patch 200). That stage (`md-quilt`) is where the v3.27.0 ex17
RATTLE failures and the PIP2 build's rattle2 crash both happened. The piercing machinery of
3.29.x (rotation, slide, spin, the deferral ceiling, the keep-best pass) exists to manage what grid
placement creates.

## What existed before, and why it was retired

Until June 2026 the quilt was tiled from patches (`bilayer_quilt.tcl`, last at `3eaad886^`). It was
retired in `3eaad886` (symmetric), `ef2032ec` (asymmetric) and `c0141b98` (asymmetric + embedding)
for the grid packer, on speed (a POPC test went from 2:00 to 1:16), tiling seams, and fitting the
box to a protein. None of those commits weighed the loss of the patches' packing.

The old asymmetric quilt also had a real flaw, which this design does not repeat. It tiled **both**
leaflets on one lattice spaced by the **larger** patch's box, so the leaflet from the smaller patch
had a strip of empty space between every pair of tiles, and it then equalized areas by deleting
randomly chosen lipids from the larger leaflet, punching holes through its interior.

## The design

### Leaflets come from their own patches, on their own lattices

- **Symmetric build:** one patch; both leaflets from it.
- **Asymmetric build:** patchA (symmetric, upper composition) supplies the upper leaflet; patchB
  (symmetric, lower composition) supplies the lower. Unchanged from today: both are relaxed with
  the `patch` protocol at 100 lipids per leaflet by default (`patch_nlipids`).

Each leaflet is tiled with **its own** patch's lateral box vectors, so within a leaflet the tiling
is seamless: a periodic patch's images already match at its edges.

### Crop to the target box

1. **Target box** `Lx × Ly`: the protein footprint plus `embed.xydist` when embedding (as today,
   `_protein_box_dims`); otherwise `npatch` tiles of the **upper** patch (asymmetric) or of the
   patch (symmetric), or `dims` if given.
2. **Tile** each leaflet `ceil(Lx/ax) × ceil(Ly/ay)` times, where `ax, ay` are that patch's box
   lengths.
3. **Keep** a lipid if its reference atom (the headgroup atom used for leaflet assignment) lies in
   `[g/2, Lx - g/2) × [g/2, Ly - g/2)`, where `g` is a small seam gap. Every other lipid is discarded.
   This is where the "excess lipids" go: lipids from partial tiles are built and thrown away. None
   reach the system.

Leaflet counts are then **whatever the patch density puts in the box**. `n_upper ≈ Lx·Ly / APL_A`
and `n_lower ≈ Lx·Ly / APL_B` fall out of the cropping; no APL number is computed or used to size
anything. The count ratio is the stress-free ratio of the two patches' states, exactly as now.

**Seams.** Where the target box does not hold a whole number of tiles, the crop edge meets the
periodic image of the opposite edge, so each leaflet has one seam in x and one in y. The gap `g`
makes that seam a slight gap instead of an overlap (a lipid's tails lean laterally by a few Å, so a
center-only crop with `g = 0` overlaps tails across the seam). A gap closes in NPT; an overlap is a
clash, and possibly a threading. Seams are the only place the tiled quilt is not at the patch's
packing. `g` is a parameter, `bilayer.quilt_seam_gap`; its default is to be measured (start from
2 Å). When the box is a whole number of tiles, as for a symmetric build without embedding, there is
no seam and `g` is not applied.

### Mechanics

- **Wrap first.** NAMD writes unwrapped coordinates. Before tiling, wrap each patch **by residue**
  into its cell, so a lipid straddling the patch boundary is moved whole, never split.
- **Leaflet assignment** by a headgroup reference atom's z relative to the patch midplane: the P
  atom for phospholipids, O3 for CHL1, with a per-species table where needed. Do not use the center
  of mass: a sterol that has sunk toward the midplane is still in the leaflet its hydroxyl points
  out of. A lipid within ~2 Å of the midplane (a flip-flop in progress) is reported, not silently
  assigned.
- **z registration.** Shift each patch so its midplane (mean of the two leaflets' reference-atom z)
  is at z = 0 before taking a leaflet from it. Leaflets of different thickness then meet at the
  midplane, which is where the grid path also puts them.
- **Lipids only.** Water and ions are not tiled. The assembled PDB contains lipids only, which is
  exactly what `grid_patch` hands to `do_psfgen` today, and `equilibrate_bilayer` solvates it as it
  does now. This is the integration point, and it means `bilayer_patch.tcl`, solvation and the
  relaxation path are unchanged. The tiler replaces `grid_patch` for the full membrane, nothing else.
- **Python, not Tcl.** The tiler is numpy over a parsed PDB (`pestifer/molecule/quilt.py`),
  consistent with the coordinate-transform port. Watch the wide-resname column issue for six-character
  CHARMM names (`SAPI24`): write through the same PDB writer the grid path uses.
- **Embedding** is unchanged: `embed_protein` takes `quilt_state` and deletes what overlaps the
  protein.

### Relaxation of the tiled quilt

The quilt no longer melts a grid, so its protocol should be short: minimize, a short NVT, then NPT
(or NPgT) long enough to close seams. The spliced `ring_check` guard stays as it is. The expected
piercing count after tiling is near zero (only at seams), so if the first check finds many, that is
a defect report, not routine. The `quilt` relaxation protocols in examples 16 and 17 were written to
melt a grid, and should be cut back once measured.

### What this does not fix

**The count ratio is still only as good as the patches' state.** A patchA stopped at 47.5 Å², like
rep-01's, gives an upper leaflet about 14% sparser than rep-03's would, now physically packed but
sparse. The tiled quilt keeps relaxing as one system, and production NPT continues that. What it
cannot do is change the number of lipids in each leaflet, and that number is fixed at the crop.
Measuring the result (`diagnose_differential_stress`) is the check on it. The docs must say
plainly that a build delivers a physically packed **starting** membrane with the patch-derived leaflet
ratio, not an equilibrated asymmetric bilayer.

**Tiling replicates the patch's lateral structure** with period equal to the patch length
(~65-70 Å at 100 lipids per leaflet), including any incipient domains or cholesterol clusters. This
is the standard price of tiling. It decays during production, and it is a reason not to make the
patch smaller.

## Costs

- **Symmetric builds gain a patch stage.** They currently grid the full membrane directly. Under
  this design they relax one 100-per-leaflet patch first (tens of minutes for a fluid mixture), then
  tile. The quilt's own relaxation should shrink by more than that, but this is to be measured, not
  assumed.
- **Lipids discarded at the crop** are computed and placed but never simulated, so they cost nothing
  that matters.
- **Grid-path code that only manages grid-placement consequences** (the first check's
  `on_unresolved: continue`, `quilt_grid_slack`, the deferral ceiling) becomes dead for the full
  membrane. The patches are still gridded, so it stays for them.

## Decisions (2026-10-26)

1. **The calibration gate warns under tiling; it no longer stops the build.** An unconverged patch
   is still the best-packed thing available. Its drift is logged and recorded in the run record,
   and `diagnose_differential_stress` is where a wrong leaflet ratio is measured.
2. **`bilayer.leaflet_apl` is retired.** Nothing is sized from an APL, so it has nothing to inject
   into. Setting it raises a clear error saying why.
3. **The full membrane may still be gridded directly, as an opt-out** (`bilayer.quilt_method:
   grid`; the default is `tile`), while tiling is validated. The patches are gridded either way.
   On the grid route the calibration gate stays fatal and `leaflet_apl` keeps working, because that
   route still sizes the membrane from an APL; decisions 1 and 2 apply to the tiled route. The
   opt-out is to be removed once tiling has carried the examples through a sweep.
4. **`quilt_seam_gap` starts at 2 Å**, to be confirmed by measurement.

## Validation plan

1. **Unit tests (pure numpy, run in CI):** wrapping by residue; leaflet assignment including a
   sunk sterol and a midplane lipid; crop counts equal to `round(area · density)` within one
   lipid; whole-tile boxes produce no seam; a seam gap of `g` leaves no pair of reference atoms
   across the seam closer than the in-tile nearest-neighbor distance minus a tolerance.
2. **Negative control:** tile with `g = 0` and confirm the seam check flags overlaps, so the gap
   test can fail.
3. **Builds, local:** symmetric POPC (seamless case); an asymmetric POPC / POPC:CHL1 build; the
   piercing count after tiling, which must be near zero and only at seams. Check with the
   independent segment-plane test, not only `ring_check`.
4. **Examples 16 and 17, then a sweep:** quilt relaxation survives with a shortened protocol;
   leaflet counts match the patch densities; wall time compared with the 3.29.8 sweep.
