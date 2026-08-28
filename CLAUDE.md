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
git archive HEAD | tar -x -C /tmp/cisim && cd /tmp/cisim
PATH="$(echo "$PATH" | tr : '\n' | grep -v /usr/local/bin | paste -sd:)" \
    python -m pytest tests/unit -q
```

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
