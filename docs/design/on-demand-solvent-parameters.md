# Design: parameter provisioning for on-demand solvent boxes

Status: **proposed** — promoted from the roadmap; investigation done, Sub-problem-B policy
decided (B2 — detect and refuse; B1 rejected).

## Problem

On-demand solvent-box generation (v3.5.0/v3.6.0) packs a fresh box of a requested solvent and
equilibrates it under NPT (`make_solvent_box` in `pestifer/charmmff/make_solvent_box.py`). For some
solvents the NPT step dies with a `UNABLE TO FIND ... PARAMETERS` fatal. Observed building Example 24:
**acetone (`ACO`) builds; acetonitrile (`ACN`) fails** with

```
UNABLE TO FIND DIHEDRAL PARAMETERS FOR HGA3 CG331 CG1N1 NG1T1
```

## Investigation — what is actually happening

- `make_solvent_box` appends the residue's **defining topology file** (`meta['charmmfftopfile']`) to
  both the psfgen topology list and the NAMD parameter list (lines ~287 and ~349). It does **not**
  separately discover the residue's *parameter* file(s).
- But the NAMD equilibration already loads the **full default standard set** — all five base
  parameter files, including `top_all36_cgenff.rtf` / `par_all36_cgenff.prm` (see `schema/base.yaml`
  and `namd.py:fetch_standard_charmm_parameters`). So for a base CGenFF residue the CGenFF parameters
  *are* present.
- For **`ACN`** specifically: the residue is defined in `top_all36_cgenff.rtf`, and
  `par_all36_cgenff.prm` (which is loaded) carries its **bond** (`CG1N1 CG331`) and both **angles**
  (`CG331 CG1N1 NG1T1`, `CG1N1 CG331 HGA3`) — but **no dihedral** for `HGA3-CG331-CG1N1-NG1T1`. That
  torsion is absent from the **entire bundled force field** — no explicit line and no `X CG331 CG1N1 X`
  wildcard. This is not an accidental gap: acetonitrile's C–C≡N heavy-atom skeleton is **linear**
  (collinear, C₃ᵥ molecule), so the reference nitrogen `NG1T1` lies *on* the `CG331–CG1N1` bond axis and
  the H–C–C≡N dihedral angle is **geometrically degenerate / undefined**. CGenFF correctly omits a
  parameter for a meaningless torsion — but psfgen's `autogenerate dihedrals` still emits the dihedral
  from connectivity, and NAMD refuses to run a term with no parameter. The term should not need a
  parameter at all; there is nothing physical to fill in.

So the roadmap's original framing ("a CGenFF solvent whose parameters live outside the default set
fails; discover and add the RESI's parameter file") is only half the story. The work splits into two
distinct sub-problems.

### Sub-problem A — residue parameters in a non-default file (tractable)

A solvent defined in a collection/stream file that is **not** part of the default standard set: appending
only `meta['charmmfftopfile']` can leave its bonded parameters unloaded when they live in a **separate
`.prm`** (the classic `rtf`+`prm` split) or in a stream not in the default set. (A self-contained `.str`
already carries both topology and parameters, so that case happens to work today; `rtf`+`prm` and
non-default streams do not.)

**Fix:** when generating a box, discover the residue's **full** topology+parameter file set — its
defining stream/rtf plus any companion parameter file(s) — and add all of them to the psfgen topology
list *and* the NAMD parameter list before equilibrating. This mirrors, and generalizes, the existing
single-topfile append.

### Sub-problem B — genuine force-field gap (needs a policy decision)

A solvent whose connectivity generates a bonded term that has **no parameter anywhere in the bundle**
(`ACN`'s methyl–nitrile dihedral). No amount of file discovery fixes this.

**Decided: B2 — detect and refuse.** After the box psfgen, surface any missing term up front with an
actionable message (the offending residue and atom-type quartet, and that no shipped parameter defines
it) instead of a cryptic NAMD fatal mid-equilibration. Do **not** fabricate parameters.

Rejected alternatives:

- **B1 — inject a zero-force-constant term** (e.g. `X CG331 CG1N1 X 0.0 1 0.0`). Rejected: for the
  `ACN` case the torsion is geometrically degenerate (see above), so a zero-`k` line just papers over a
  meaningless term; more generally, silently editing the force field to make a build proceed risks
  masking a real error and is exactly what B2 refuses to do.
- **B3 — declare such solvents flatly not auto-boxable.** Subsumed by B2: the clear up-front message is
  the "not auto-boxable as shipped" signal, and it tells the user precisely what is missing.

Note for a possible future follow-up (out of scope here): the deeper cause for `ACN` is that
`autogenerate dihedrals` emits a dihedral across a linear (degenerate) skeleton. A cleaner long-term fix
would suppress generation of such degenerate torsions rather than requiring a parameter for them — but
that reaches into the psfgen dihedral-generation path and is well beyond on-demand solvent boxing.

## Proposed plan

1. **Parameter-set discovery (Sub-problem A).** Given the solvent RESI, resolve its complete defining
   file set (topology + parameter/stream files) from the force-field metadata and add every file to
   both scripters' standard sets in `make_solvent_box`, not just `meta['charmmfftopfile']`. Covers the
   `rtf`+`prm`-split and non-default-stream cases.
2. **Dry-run parameter check.** After the box psfgen (before NPT), enumerate the bonded types the PSF
   actually requires and diff them against the loaded parameters. This both (a) confirms A is fixed and
   (b) surfaces B **precisely** (exact missing term + residue) instead of as a mid-run NAMD fatal.
3. **Sub-problem B (B2).** Turn the dry-run diff into a clear pre-equilibration error naming the
   residue and the unparameterized atom-type term(s); do not fabricate parameters.
4. **Regression.** Build end-to-end: `ACN` (the FF-gap case) and at least one `rtf`+`prm`-split /
   non-default-stream solvent (the Sub-problem-A case). Tie into the on-demand cache so a first build
   pays the cost once.

## Notes

- Discovered building Example 24 (subtilisin in acetone vs. acetonitrile).
- Related roadmap item "Solvent blends / mixtures" will union multiple components' parameter files —
  it should reuse whatever discovery mechanism Sub-problem A introduces here.
