# Working on pestifer

Notes for an agent (or a person) changing this repository. For how to *use* pestifer, see the
docs; for contributor setup and how to run the test suites, see `docs/source/contributing.rst`,
which is the maintained account -- this file only records things that are easy to get wrong and
that nothing else will tell you.

## A local test pass proves less than it looks

CI runs `tests/unit` on GitHub runners that have **no VMD and no NAMD**. Tests needing them are
marked `needs_tools` and skip there, so a developer machine with the toolchain runs roughly 200
tests that CI never sees -- and, more dangerously, a fresh checkout lacks untracked state that a
working machine has accumulated.

Two separate failures have reached `main` this way: fixture symlinks committed with absolute
paths, which become dangling links in any other checkout, and a per-module test working directory
that was never tracked, so `conftest.py`'s `chdir` fell back a level and every path in that module
resolved wrongly. Both passed locally, for the same reason each broke everywhere else.

Before pushing anything that could affect CI, run it the way CI will see it -- a clean export,
with the tools hidden:

```bash
export CISIM=~/devtests/pestifer/cisim && rm -rf $CISIM && mkdir -p $CISIM
git archive HEAD | tar -x -C $CISIM && cd $CISIM
PATH="$(echo "$PATH" | tr : '\n' | grep -v /usr/local/bin | paste -sd:)" \
    uv run --extra test pytest tests/unit -q
```

The export goes under `~/devtests/pestifer/` rather than `/tmp`: `/tmp` is shared with every other
session's scratch, which makes anything there unsafe for anyone else to clean up, and a `/tmp` that
fills has taken builds down with it.

`uv run` is what makes this a clean *environment* and not just a clean source tree: run under a
bare `python`, the export borrows the repo's own `.venv`, so any package hand-installed there over
the months is invisible to the check and absent on the runner -- the same accumulated-state
failure the export is meant to catch. uv ignores the active venv (it says so: `VIRTUAL_ENV ... will
be ignored`) and resolves the export's `pyproject.toml` into a fresh one. The cost is that this
check now needs the network, and a resolution failure is a different red than a test failure.

## A check that passes is not evidence the check ran

The section above is one instance of a wider failure, and the wider one cost more time here in a
week than any bug did. **A command can complete successfully without checking the thing you
meant**, and the result then looks exactly like a pass. Seven of these in early September 2026,
most of them self-inflicted (two more added since):

| the check | what it actually verified |
| :--- | :--- |
| `uv run --with /path/to/ycleptic` + a green suite | PyPI's ycleptic. The `>=2.3.0` floor was already satisfied so uv never used the path, and the branch reported the same version number. |
| a docs build reporting `0 warnings` | nothing. sphinx exited 127 -- its venv had been cleaned out of `/tmp` -- and zero warnings meant zero of anything. |
| a CI waiter reporting "finished" | a failed `gh` query. The loop treated an empty result as "nothing incomplete". |
| conformer-cache tests passing | the helper function, not that the cache guard *calls* it. Reverting the call site left every test green. |
| the SLURM cpu-detection test passing | nothing: this box has 24 cores, the expected value, so `os.cpu_count()` returned the right answer and the test passed against the bug it existed to catch. |
| `git push`, exit 0, "Everything up-to-date" | that a detached HEAD had nothing to push. The commit was not on `main`. |
| a peer's `grep -E` positive control | the pattern's *syntax*. The alternation was written BRE-style with an escaped pipe, which in ERE matches a literal pipe character, so every search could only return zero -- and the control still passed. |
| "0 unset coordinates" on a sugar, grepping the psfgen log for `AGALNA` | nothing. psfgen logs coordinate warnings under the residue's *PDB* name (`A2G`); the grep could not match, so the control without the fix also read 0. (2026-09-14) |
| VMD selecting the right number of lipid atoms after regenerating `macros.tcl` | VMD's *built-in* `lipid` keyword. One hyphenated name made VMD reject pestifer's whole macro and fall back silently; the count was identical. Read the macro back with `atomselect macro lipid`. (2026-09-14) |

Two habits catch all of them.

**Gate on positive evidence of the specific claim, never on the absence of a contrary one.** "No
failures appeared" is not "it ran". Before trusting a dependency test, diff the installed file
against the source you meant to test -- a version string cannot tell you (ycleptic's branch and
its release both said 2.3.0). Before trusting a wait, require the thing you are waiting for to
say it is done. Before trusting a push, read its output, not its exit code.

**Prove the check can fail, at the level the fix actually lives.** Revert the fix and confirm the
test goes red. Rows 4 and 5 are the sharp ones: in both, a negative control existed and was
aimed one level off -- at the helper rather than its call site, and at affinity rather than
`cpu_count` as well. A negative control that cannot fail is the same bug as the one it is
guarding against.

This is why the clean-export gate above is written as it is, and why
`tests/unit/test_core/test_processor_info.py` is deliberately **not** marked `needs_tools`: it
is pure logic over environment variables, and it is exactly the code that misbehaves on a
cluster CI cannot reach, so it has to run everywhere. Marking a test `needs_tools` that does not
need them is how a guard against silent failure becomes silent itself.

## The test suite dirties the working tree

Several tracked files under `tests/unit/test_tasks/` are generated test output, rewritten on every
run -- and, since generated scripts carry a version watermark, they change on every release too.
`git status` is therefore dirty after running the suite.

**Stage by explicit path. Never `git add -A` after running tests.**

## Do not build inside the repository

Builds write their outputs into the current directory. Run them somewhere else; the sweep scripts
default `OUTROOT` to the current directory, so set it or `cd` first. A build started in the repo
root scatters system files through the working tree and rewrites the tracked fixtures above.

## Releases

Cut releases with `./scripts/release.sh <version>`. It gates on the integration tests, rotates the
CHANGELOG, bumps `pyproject.toml` and `CITATION.cff`, tags, and pushes. Do not hand-roll any of
that: **pushing the tag is what publishes to PyPI**, so a hand-made tag publishes without the
gate.

Do not release on a red badge. The tag is permanent; the failure gets baked into it.

## Fixed: packaged NAMD config pointed at a parameter file the tarball did not contain

Found 2026-08-21 against v3.19.1 while benchmarking a packaged system on Picotte; fixed the same
day. Kept here because the shape of the bug -- two names for one file, agreeing everywhere except
inside the tarball -- can recur anywhere packaging renames things.

A `terminate` task with both a `basename:` and a `package: basename:` wrote a tarball whose NAMD
config referenced the *package* basename for its parameter file, while the file actually tarred
carried the *terminate* basename:

```
$ tar tzf prod_vansc_native_r1.tar.gz
...
prod_vansc_native_r1/prod_vansc_native_r1.namd
prod_vansc_native_r1/vansc_native_r1_minimal.prm     <-- shipped

$ grep ^parameters prod_vansc_native_r1.namd
parameters prod_vansc_native_r1_minimal.prm          <-- referenced
```

`structure`, `coordinates`, and `extendedSystem` were all correct; only `parameters` was wrong,
and NAMD aborted on the missing file. It stayed invisible in the build directory, where a
consolidated `.prm` under the package basename also exists and everything resolves.

`NAMDScripter.consolidate_params()` (`pestifer/scripters/namd.py`) named its output from the
scripter's *current* basename, which inside `TerminateTask.make_package()`
(`pestifer/tasks/terminate.py`) has deliberately been set to the package basename, while
`TarballContents.append(min_artifact)` tars the artifact named by `copy_state_to_basename()` --
the terminate basename. `make_package`'s docstring says which one is intended ("State files are
included in the tarball under their existing names (the terminate basename); only the tarball
itself and the NAMD config script use the package basename"), so the fix is on the
`consolidate_params` side: when the parameter set is already a single `*_minimal.prm`, keep that
file under its own name instead of re-deriving one from the basename in force.

Two regression tests guard it:
`tests/unit/test_scripters/test_namdscripter.py::TestConsolidateParams` pins the naming rule, and
`tests/unit/test_tasks/test_terminate.py::TestPackagedConfigIsSelfContained` runs a real
`make_package` and asserts every file named in the packaged `.namd` is present in the tarball --
the invariant that would catch the next one of these too.

## Fixed: `new-system --inspect` chain summary counted waters and ions as protein

Found 2026-08-21 against v3.19.1, alongside the packaging bug above; fixed 2026-08-27. Kept
because the shape of it -- a label and a number computed over different sets -- is easy to
reintroduce anywhere a chain id is treated as one molecule.

The scaffold config emitted by `pestifer new-system --inspect` summarizes each chain like this,
for PDB 8DX0:

```
#   A: protein (263 residues) — HISTIDINE KINASE
#   B: protein (246 residues) — HISTIDINE KINASE
```

Neither number is a residue count of the protein. Both chains contain **139** amino-acid
residues. The printed figures are the count of *every* distinct residue in the chain, waters
and ions included:

```
chain A:  139 protein + 123 HOH + 1 MG = 263
chain B:  139 protein + 105 HOH + 2 MG = 246
```

So the label says "protein" and the number counts solvent. Reproduce with any structure
carrying chain-tagged waters:

```bash
pestifer new-system 8dx0 --inspect
```

Cross-checks for 8DX0 chain A, none of which agree with 263: 139 resolved (ATOM records),
15 missing per REMARK 465 (208-209, 313-322, 359-361), 154 in SEQRES/DBREF (208-361), and 149
in the system pestifer actually builds (210-358, 139 from the crystal plus the 10-residue loop
it rebuilds).

Why it was worth fixing rather than documenting: `--inspect` output is advisory text a user reads
once and copies into a methods section, so a wrong count propagates silently into writing and is
never checked again. It reached a published page that way before being caught.

`_chain_identities()` (`pestifer/core/system_inspector.py`) now buckets residues **by segtype**
before counting anything, and `ChainIdentity` carries a `composition` dict so `describe()` can
name what else shares the id rather than folding it into the headline number:

```
#   A: protein (139 residues; also 123 water, 1 ion) — HISTIDINE KINASE
#   B: protein (139 residues; also 105 water, 2 ion) — HISTIDINE KINASE
```

Two further faces of the same bug were fixed with it, neither in the original report:

- **`resnames` was sampled across the whole chain too**, so a glycan chain carrying waters
  advertised `glycan (NAG, BMA, HOH...)`. The sample is now segtype-restricted.
- **The segtype itself was chosen by counting *distinct resnames***, not residues. That is why a
  protein chain still read "protein" despite the solvent (20 amino-acid names beat one `HOH`) --
  but a DNA chain carrying five kinds of ion would have been classified an *ion* chain, since DNA
  has only four resnames. Classification is now: any polymer residues win; only among non-polymer
  chains does residue count decide.

`tests/unit/test_core/test_system_inspector.py::TestChainComposition` pins all three against
synthetic atom tables (no network), including the 8DX0 shape.

## A wrong-but-valid `type:` in the schema silently disables a whole subtree

Found 2026-08-27, chasing a report that a `validate` task with a mistyped `measure:` logged an
error and emitted no check. The class-level fallthrough was real, but it was the second problem.
The first: `validate` was declared

```yaml
      - name: validate
        type: list          # <-- its payload is a mapping
```

Tasks are elements of the `tasks` list, and ycleptic's `lwalk` dispatches a list element on its
declared type: scalars are ignored, `dict` descends via `dwalk`, and **anything else falls to a
debug-level "ignored"**. So the entire `validate` subtree was never walked. Every `choices:`
under it was inert -- `measure`, `connection_type`, all of it -- and a typo reached the task
untouched. `validate` was the only task in the schema not declared `type: dict`.

What makes this worth recording is that nothing flags it:

- **`yclept check-spec` passes.** It verifies that keys and type *names* are recognized, not that
  a declared type matches the shape it describes. It reported "no unrecognized keys or types"
  before and after the fix.
- **Every example config still parsed**, before and after -- a skipped subtree raises nothing.
- The one visible symptom was in a code path far away: an unsupported `measure` reaching
  `ResidueTest`, which could only have happened if the schema had not rejected it first.

Diagnosing it requires *walking* the schema, not reading it. The regression test does exactly
that: `tests/unit/test_tasks/test_validate.py::TestSchemaEnforcesTheSameSpecsAsTheCode` calls
ycleptic's `dwalk` on the real `base.yaml` with a deliberately mistyped spec and requires it to
raise -- and separately asserts each `choices:` list equals the matching class-level `*_supported`
set, so the two gates cannot drift apart.

If you add a task to the schema, declare it `type: dict`, then prove a bad value in it is
rejected. A passing parse is not evidence that anything was checked.

## Validation tests must not be written against chain or segment letters

Letters in a built system are assigned by pestifer, not inherited from the input. Each input
chain is split into one segment per segtype, and each segment past the first takes the next
*unused* letter -- so excluding a chain frees its letter for the next segment that needs one.
Verified 2026-08-27 by building 8DX0 both ways:

```
no exclusion                      exclude: [chainID == 'B']
A protein / B protein             A protein
C MG      / D MG                  B MG      <-- chain A's magnesium
E water   / F water               C water
```

Both `segname B` and `chain B` name that magnesium, so the obvious test -- "I excluded chain B,
so chain B should be empty" -- reports FAIL on a correct build. Whether you get away with it
depends on where the excluded letter sits in the pool: example 5 excludes chain `P` and tests
`chain P`, which is safe by accident, not design. Guidance is in
`docs/source/subs/buildtasks/validate.rst` ("Test the molecule, not the letter").

Note the reported cause was "the chainIDmanager reassigns segids" -- i.e. that retained chains get
renamed. They do not. Chain A stays A; it is the *derived* water and ion segments that draw from
the freed pool. Same symptom, different mechanism, and the difference decides the fix.


## Partner order in a Link is an invariant, not a convention two places happen to share

Found 2026-09-09 from a Rosetta-written PDB. Both the PDB `LINK` convention and the CHARMM
carbohydrate `PRES` definitions put the **anomeric carbon second** -- `O4 -> C1`, `ND2 -> C1`,
`O6 -> C2` for sialic acids. Other producers sometimes write the reverse. Nothing about the bond
changes, and the atom identifiers are perfectly good, but two independent things downstream read
the order and *both* fail silently:

- `Link.set_patchname()` keys every glycan branch on `name2` being the anomeric carbon. A
  reversed link matches nothing, falls to the terminal `else`, and ends `UNFOUND`. At
  `scripters/psfgen.py:848` that writes a comment and **no `patch` line**, so the build
  *succeeds* with the glycosidic bond simply absent from the PSF.
- `Residue.link_to()` is directional -- `self.down.append(other)`, so partner 1 is the parent and
  partner 2 the child. A reversed link builds that branch of the glycan tree **upside down**, and
  anything walking it with `get_down_group()` works from the wrong end. The mutation-driven
  pruning in `segment.py` then deletes the wrong subtree, and that one does not even warn.

Both consumers sit downstream of `LinkList.assign_residues`, so one canonicalization there fixes
both: `Link.canonicalize_glycan_orientation()` is called after segtypes and resnames are set and
**before** `link_to` and `set_patchname`. Do not move it after either.

Two things not to redo the hard way:

- **Do not fix this in `set_patchname` alone.** It is the visible half. The inverted tree is the
  worse half precisely because it warns about nothing.
- **Do not hand-list the fields to swap.** `Link` carries 16 partner-indexed pairs in two naming
  shapes (`chainID1`/`chainID2` and `ptnr1_label_asym_id`/`ptnr2_label_asym_id`), and a partial
  swap is how "the atom identifiers are fine" stops being true. `Link._paired_fields()` derives
  them from `model_fields`, so a pair added later cannot be forgotten.

`segment.py`'s glycan-adjacency builder was already order-agnostic (it handles either partner
being the protein); the pruning code was not. That inconsistency is why normalizing at the
boundary beats teaching each consumer both orders -- every future consumer gets it for free.

`tests/unit/test_objs/test_link.py::TestLinkOrientationAgainstRealStructure` pins the invariant
against 4zmj's 25 real glycan links across seven patch types: the same bonds written either way
round must give the same patches *and* the same parent/child directions. It asserts the canonical
run resolved real patches, so the comparison cannot pass as two piles of `UNFOUND` agreeing with
each other. With the canonicalization disabled, all 25 come back `UNFOUND` with inverted parents.

## An xsc that exists is not a promise of a periodic cell — and sweep by call site, not module

`cell_from_xsc` returns `(None, None)` for any xsc it cannot get a cell out of. The common case is
not a corrupt file: a run with **no periodic boundaries** — a vacuum minimize — writes a perfectly
valid origin-only xsc, `step o_x o_y o_z`, 4 columns, where a cell needs at least 13. So

```python
box = cell_from_xsc(xsc)[0] if xsc is not None else None   # None for a cell-less xsc
if xsc is not None:                                        # still True: the path exists
    sidelengths = np.diagonal(box)                         # np.diagonal(None) -> ValueError
```

is wrong in a way that reads as right. **Guard on the parsed cell, never on the path.** Every other
`cell_from_xsc` call site in the tree already does (`if box is not None`); the two that did not were
`make_membrane_system` (fixed 2026-08-27, `c9ab48fd`) and `RingChecker.check` (fixed 2026-08-28).

The second one is the part worth remembering. The 08-27 fix came with a sweep of the other call
sites, and that sweep **cleared `tasks/ringcheck.py:218`, correctly** — it does guard. But
`ringcheck.py` is the task wrapper; it only passes the xsc path down. The site that *consumes* the
cell is `RingChecker.check()` in `psfutil/psfring.py`, one layer below and in a different module,
and it was on nobody's list. Clearing the wrapper closed the question a layer too high, and the
same build hit the same defect the next day at the next call site down.

Sweep this pattern with `grep -rn cell_from_xsc` and check every hit, including the ones in modules
you have already decided are fine. A module is not a unit of correctness here; a dereference is.

`tests/unit/test_tasks/test_ringcheck.py::TestRingCheck::test_ring_check_cell_less_xsc` pins it,
with an xsc that *exists and parses* but yields no cell — the pre-existing non-periodic test passes
`xsc=None`, takes the other branch, and would never have caught this.

One unguarded deref remains, deliberately: `make_solvent_box.py:432` indexes `final_box[0][0]`
straight from `cell_from_xsc`. It reads the last xsc of a solvent-box NPT equilibration that
pestifer itself just ran, so a missing cell there means the pipeline is already broken and there is
no user input that reaches it. Left alone rather than papered over.

## CHARMM atom types are case-insensitive, and the shipped release relies on it

Found 2026-09-10 chasing "phosphotyrosine has no CHARMM parameters", which was wrong. A single
shipped file, `toppar_all36_prot_na_combined.str`, spells one type three ways:

```
MASS  -1  ON2B     15.99940 O ! ...        <- the type declaration
ATOM  OH  ON2B     -0.36                   <- RESI PTR (and PRES TP2) -- what reaches a PTR PSF
ATOM  OH  ON2b     -0.36                   <- PRES TP1 only
CA    ON2b  340.0   1.38                   <- every bonded parameter
ON2B  0.0  -0.1521  1.77                   <- the vdW record
```

**psfgen does not normalize the case.** An earlier version of this note said psfgen writes the
`MASS` spelling into the PSF. It does not: pestifer runs psfgen under `psfcontext mixedcase`,
where a residue's `ATOM` type is taken verbatim and must match a `MASS` record exactly. A PTR PSF
carries `ON2B` because `RESI PTR` spells it that way. `PRES TP1`, which spells it `ON2b`, cannot be
applied at all -- psfgen stops with `unknown atom type ON2b` / `MOLECULE DESTROYED BY FATAL ERROR`
(found 2026-09-12 building every PTM route; `TP2`, spelled `ON2B`, builds). CHARMM itself would
accept `TP1`. Across the release that mismatch otherwise appears only in model compounds and in
lowercase-named protonation patches in a stream pestifer does not load.

CHARMM does not care. Anything in pestifer that compares an atom type as a Python string does,
and gets it **half** right, which is worse than getting it wrong. A real built phosphotyrosine
PSF carries `ON2B` (checked against one, rather than reasoned about -- the first version of this
note had the direction backwards, and the second had the reason wrong; see above). So its *vdW* lookup succeeds and every one of its *bonded*
terms fails: `('CA','ON2B')` is not `('CA','ON2b')`. That is precisely why
`extract_for_atomtypes` used to drop phosphotyrosine's parameters from the minimal file while
leaving a result whose own counts looked self-consistent.

What is actually in the shipped release, measured rather than assumed: **56** `MASS`-declared
types contain a lowercase letter, and every one is a metal ion (`Ag1p`, `Fe2p`, `Ni1p`, ...
all ending in `p`); `ON2B`/`ON2b` above is the case *split*, where one type is declared, used
and parameterised in three different spellings. An earlier version of this note, and the code
comment in `extract_for_atomtypes`, claimed "13 types ... including `Br` and `Cl`" and a
lowercase `x` wildcard. Both are wrong and were corrected 2026-09-11: `Br`/`Cl` on those `MASS`
lines are the *element-symbol* column, not the atom type (the types are `BRGA1`, `CL`, both
upper-case), and every lowercase standalone `x` in the release sits after a `!`, inside a
comment. Check the column before counting a token as a type.

**Upper-case both sides of every atom-type comparison.** The fixes so far:
`CharmmParamFile.extract_for_atomtypes` and its five `_*_key` dedup functions
(`charmmff/charmmffprm.py`), and every comparison in `charmmff/psf_param_check.py`.

The instructive part is the second one. `extract_for_atomtypes` was fixed in `5b2b6fd8`; the
`_key` functions were fixed in the *same file* at the same time -- but only `_bond_key` and
`_angle_key`, leaving `_dihedral_key`, `_improper_key` and `_nbfix_key` case-sensitive. And
`psf_param_check.py` builds its own lookup structures from the same `CharmmParamFile` and was
not touched at all, so it shipped the identical false positive to `ContinuationTask` for a
month. This is the `cell_from_xsc` lesson again in a different costume: a module is not a unit
of correctness, and neither is a commit. Sweep with
`grep -rn 'type1\|atomtype' pestifer/charmmff/` and check every *comparison*, not every file.

## The consolidated .prm is verified against the PSF before NAMD sees either

`NAMDScripter.consolidate_params` writes `{basename}_minimal.prm` and then calls
`_verify_params_cover_psf`, which raises `PestiferBuildError` if the file does not resolve every
atom type and every bond/angle/dihedral/improper type-tuple in the PSF. The engine is
`charmmff/psf_param_check.py`, already written for incoming foreign PSFs and now used on the
normal build path too.

Why raise rather than warn: a missing parameter is fatal to NAMD regardless, so the only
question is whether it fails somewhere it can be diagnosed. NAMD reports one term, by atom
serial, after minimization setup:

```
FATAL ERROR: UNABLE TO FIND ANGLE PARAMETERS FOR C CTL2 CTL2 (ATOMS 8 10 13)
```

The check reports all of them, by residue, before NAMD launches -- six terms for that same
1UPH build, because NAMD dies on the first one it meets and hides the rest.

The two causes are separated because the remedies are opposite. Missing from the *merged* set
means no file in `self.parameters` defines it, and the build needs another stream: the case
that motivated this is `GLYM`, defined in `toppar_all36_lipid_prot.str` but taking
`C CTL2 CTL2` from `toppar_all36_lipid_sphingo.str`, a file with no other reason to be loaded.
Missing from the consolidated file *only* means `extract_for_atomtypes` dropped a record it
should have kept -- a pestifer bug, and the message says so, because otherwise it gets debugged
as a user configuration problem.

Validated before shipping the raise, since a false positive here stops every build: 63 real
builds under `~/devtests`, 62 clean, and the one flagged was a build NAMD had already killed on
the same term -- a THR whose CB a patch retyped `CT2` while leaving HB as `HA1`, and `CT2 HA1`
exists nowhere in the release.

CMAP cross-terms are checked too, added 2026-09-11 after the rest. They are matched on the full
8-type tuple, **exactly**: the shipped release's six CMAP records carry no wildcards, and the
tuple is a phi quartet followed by a psi quartet, so it is directional in a way a bond or angle
is not -- reversing it names a different thing. Every cross-term in every build swept matched a
record exactly, with no reversal needed, which is the evidence for matching that strictly.
`CharmmParamFile.merge` also gained a `_cmap_key`; it was the one section merged by
`list.extend` with no dedup key, so an overlapping stream would have duplicated records. No
duplicate arises from the default file set, so that one is latent, not a fix for an observed
failure.

## A residue's terminal patch is decided by topology load order unless something says otherwise

Found 2026-09-12 building example 29. psfgen gives each residue the `DEFA FIRS ... LAST ...`
default **in force when that residue's definition is read**, not when the segment is built. Most
CHARMM files set their own default, so this never shows. `toppar_all36_lipid_prot.str` does not
-- its `DEFA` line is commented out -- so `GLYM`, `LYSM`, `CYSP` and the rest inherit whatever
the previously loaded file left behind. For `GLYM` that matters: it acylates its own backbone N
and must start a chain with no `NTER`.

pestifer auto-loads that stream last, after `top_all35_ethers.rtf` ends on `DEFAULT FIRST NONE`,
so every build was correct and nothing hinted it was an accident. Built directly with
`prot_modify_res.str` (`DEFA FIRS NTER`) read just before it, psfgen produced:

```
default order:   N:NH1  HN:H                          <- amide, correct
NTER in force:   N:NH3  HT1:HC HT2:HC HT3:HC  C1:C    <- NH3 still bonded to the acyl carbon
```

with no error and no warning. The fix does not touch load order: `PsfgenScripter` writes
`first none` when a segment's first built residue is in `Labels.backbone_acylated_resnames`,
derived from `_fusible_ligands` (entries whose bond lands on `N`), after applying any mutation at
that position. Guarded in `tests/unit/test_scripters/test_psfgenscripter.py::TestBackboneAcylatedFirstResidue`.

If you reorder topology loading, or add a residue defined in a stream with no `DEFA`, check its
terminal patches in the built PSF -- not in the script, which will look fine either way.

**Known gap: `CYSL`.** It also acylates its backbone N (`BOND N CA1` -- the N-palmitoyl,
S-diacylglyceryl cysteine of bacterial lipoproteins) and lives in the same stream, but the set is
derived from the fusion table, which only has `GLYM`/`LYSM`. So an N-terminal `CYSL` is still
right by load order alone. Deriving the set from the topology instead -- a residue whose `N` bonds
to a carbon other than `CA` -- would cover it and anything added later; not done, because no build
has needed an N-terminal `CYSL`. Found 2026-09-12 by the PTM audit, where `CYSL` mid-chain
(correctly) failed the parameter check on the junction it cannot have.

## Glycan link patches are chosen by what the sugar IS, not by how the link is twisted

Found 2026-09-13 while asking whether the example set covered the week's glycan fixes. CHARMM names
each glycosidic patch by **configuration** -- axial or equatorial on the ring, and its PRES comments
say so ("(i)1->4(i-1) axial at C1 and equatorial at C4" is `14ab`). pestifer used to pick whichever
patch's reference **dihedral** was closest to the deposit's. A torsion about the glycosidic bond is
conformation; it cannot see which face of C1 the bond is on. On 4zmj it gave 45 of 54 beta-GlcNAc-Asn
links the alpha patch. The periodicity fix in `fc2b881f` made that lookup mathematically correct and
changed which wrong answers it gave -- a correct fix to the wrong method. Nothing broke, because the
alpha/beta patches in a family differ only in their IC tables, so a wrong one changes no force-field
term for a fully resolved glycan. That is exactly why it survived: nothing downstream could see it.

`Link.set_patchname` now takes the anomer from the CHARMM residue name and the acceptor position
from `_AXIAL_POSITIONS`, a table **generated** by building every CHARMM pyranose from its internal
coordinates and measuring which substituents are axial. Geometry is the fallback, not the primary,
because deposited rings are not always chairs: 4zmj's Man A5 is boat-like and reads equatorial despite
alpha handedness. A disagreement between identity and geometry is logged.

Two things to keep:

- **The oracle is CHARMM, not a textbook.** `tests/inputs/glycan_patch_reference_geometry.json` holds
  all 23 patches built by psfgen from ICs alone; every one must come back under its own name, through
  `set_patchname` itself. That test caught a real bug in the fix: CHARMM names 1<->1 links from
  residue 1 first, the reverse of every other family.
- **Do not re-pin a set of "patches the geometry selects".** The test that did so enshrined the wrong
  answer for a week, with a comment explaining why it was right. Pin per-link expectations derived
  from what the residues are.

## A residue alias renames every residue of that name -- including CHARMM's own

`_residue_aliases` becomes `pdbalias residue` in every psfgen script, and psfgen applies it to any
residue it reads with that name, including one pestifer wrote itself. So an alias from a PDB code
that is *also* a CHARMM residue name turns one molecule into another. Found 2026-09-14 before it
shipped: `GLA` is alpha-galactose in the PDB and gamma-linolenic acid in CHARMM
(`toppar_all36_lipid_detergent.str`), and "GLA AGAL" would have made every such fatty acid a sugar.

Check a new alias against `CHARMMFFContent(...).resi_to_topfile_map` before adding it.
`tests/unit/test_core/test_labels.py::TestResidueAliasesDoNotReclassifyCharmmResidues` fails if an
alias moves a derived CHARMM residue into a different segtype, and lists the four deliberate
collisions it tolerates.

`BGLC BGLCNA` was one of them, and it was a bug: it existed because pestifer's own segment PDBs keep
four columns of a residue name, so GlcNAc is written as `BGLC` -- and it turned every real
beta-glucose named `BGLC` into GlcNAc. Removed 2026-09-14. Its replacement, `molecule/resname_repair.py`,
decides by the residue's heavy atoms which CHARMM residue a four-character name was cut from, and
gives a long result the PDB code aliased to it (`BGLCNA` -> `NAG`) so it survives the writer.
**Do not re-add a name-based alias for a truncated stem**: the stub is ambiguous for ~70 residues.
