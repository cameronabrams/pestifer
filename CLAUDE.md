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

Four separate failures have reached `main` this way. Fixture symlinks committed with absolute
paths, which become dangling links in any other checkout. A per-module test working directory
that was never tracked, so `conftest.py`'s `chdir` fell back a level and every path in that module
resolved wrongly -- **twice more** on 2026-10-14, `test_tasks/test_make_membrane_system/` and
`test_core/test_tcl/`, each of which existed locally with zero tracked files in it; a `.gitkeep`
in each fixes it, and a directory the suite needs now has to carry one. And `test_merge.py` and
`test_desolvate.py`, which read their inputs out of `scratch/builds/1/` -- a **gitignored**
leftover of a build someone had once run on that machine. All passed locally, for the same reason
each broke everywhere else.

A test directory that holds only generated output, and a fixture that is a file nothing in the
repo produces, are the same defect wearing different clothes: the suite depends on state the
repository does not carry. Both are invisible to `git status` -- one because the directory is
untracked, the other because it is ignored.

Before pushing anything that could affect CI, run it the way CI will see it -- a clean export,
with the tools hidden:

```bash
export CISIM=~/devtests/pestifer/cisim && rm -rf $CISIM && mkdir -p $CISIM
git archive HEAD | tar -x -C $CISIM && cd $CISIM
PATH="$(echo "$PATH" | tr : '\n' | grep -v /usr/local/bin | paste -sd:)" \
    uv run --extra test pytest tests/unit -q
```

**This gate alone is not enough, and the reason is the subject of the next section.** It hides the
toolchain on purpose, so every `needs_tools` test skips in it -- correctly, since that is what CI
sees. If your ordinary repo run also omits `--runslow`, the 11 `slow` tests in `tests/unit` run in
**neither** place. Five of them are the membrane builds, the only tests anywhere that build a
bilayer from a composition. Two filters, each right on its own, whose intersection is empty: in
October 2026 that pair hid four real defects in `test_make_membrane_system.py` through about
twenty full-suite runs, during a week spent almost entirely on the membrane builder. Run the repo
suite as `uv run --extra test pytest tests/unit -q --runslow` (about 9 minutes), and report a run
without the flag as "N passed, with the build tests skipped" -- which is what it is.

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
| `sha256sum` on the artifacts a release should have written, diffed against a reference | nothing. The files did not exist, so `sha256sum` wrote its complaint to *stderr* and nothing to stdout -- the comparison was `"" = ""`, and it printed MATCH for three missing files in a row. (2026-10-01) |
| a same-seed rerun of example 11 coming back with no ring piercing, "confirming" a fix | nothing. `ligate`'s CCD closure sits downstream of a minimize, and parallel NAMD is not bitwise reproducible: same `.namd`, byte-identical input, `+p24` twice gives two different `.coor` (measured 2026-10-21; `+p1` twice is identical). The rerun drew a different loop conformation that happened not to pierce. A controlled probe on the one structure was needed instead. (2026-10-21) |
| `uv run --extra test pytest tests/unit -q` reporting "2353 passed, 0 failed", about twenty times in a week spent on the membrane builder | the 2367 tests that are **not** the membrane builds. Those five are `slow`, so they need `--runslow`; the clean-export gate that would otherwise run them hides VMD and NAMD on purpose. Two filters whose intersection is empty, and four real defects sat inside it. (2026-10-14) |

**A hash is not evidence that anything was hashed.** The last row has a sibling worth knowing
before you meet it, hit the same day in ycleptic: `curl -sL <url> | sha256sum` through a failed
fetch hashes whatever came back instead. Measured 2026-10-01, and the measurement is the point
-- the first two accounts of this, mine and ycleptic's, each named a *host* and they named
opposite ones:

| what you fetch | result |
| :--- | :--- |
| `/packages/source/<l>/<name>/<file>` -- the redirecting form a conda recipe uses; same on `files.pythonhosted.org`, `pypi.org` and `pypi.io` | `404`, **0 bytes**, so the pipe yields `e3b0c44298fc...7852b855`, the SHA256 of the empty string |
| `/packages/<hex>/<hex>/<hex>/<file>` -- the canonical direct URL the JSON API hands you | `404` with a ~370-390 byte S3 `NoSuchKey` XML body that **echoes your key and carries a per-request `RequestId` and `HostId`** -- so the hash is different on every uncached request |
| a **website** path -- `pypi.org/<anything-missing>` | `404` with a ~4.8 KB HTML page, hashing to something with no tell at all |
| a **project page** from a script -- `pypi.org/project/<name>/` | **HTTP 200**, ~3 KB, a "Client Challenge" bot interstitial -- *byte-identical for a real project and a nonexistent one* |

So the host is not the variable; the path shape is. The last row is the one that matters most:
a scripted fetch of a project page returns 200 with a plausible body whether or not the thing
exists, so **`curl -f` exits 0, a status check passes, and a size check passes**.

And row 2 kills the one reassurance left. Repeat the same failing fetch and you get a *byte-
identical* wrong hash -- six for six here -- which reads as a stable, reproducible measurement.
It is not: the response headers say `x-cache: HIT`, `age: 46`. You are re-reading the CDN's copy
of one error. Ask for a key that has never been requested and the size and hash differ every
time (389, 390, 369, 390 bytes; four distinct hashes), because the body embeds a fresh
`RequestId`. **"The hash came out the same twice" is a statement about Fastly, not about the
file.**

Spacing the repeats out does not rescue it, and the numbers are unkind. Polled every 20 s, one
key held `x-cache: HIT` with identical bytes through `age=100` and came back `MISS` with a
*completely different* hash at 120 s. So the window is about two minutes -- which straddles
exactly the interval a person re-checks something at -- and **no header tells you**: these 404s
carry no `cache-control`, no `expires`, no `surrogate-control`, and the first `HIT` arrives with
`age: 0`. There is nothing in the response to calibrate a re-check against. That is why two of us
measuring the same URL disagreed about whether the wrong value was stable: a burst of requests
sits inside one window, two runs minutes apart straddle the boundary, and both readings are
honest.

`e3b0c442...` is the lucky case: it is at least recognizable. Do not learn the constant, do not
learn the host, and do not trust repeatability as a proxy for correctness.

Gate on the content being the thing you asked for. For an archive, require the bytes to *be* an
archive (`file`, or `tar tzf`) before hashing; for an API, require the parsed field you came for.
`pestifer/util/update_check.py` does the latter -- it reads `info.version` out of parsed JSON, so
a challenge page served as 200 raises on `.json()` and is treated as "no answer" rather than as an
answer. Status and size are a cheap first filter, not the check. The same rule covers the table
row above: `sha256sum` reports a missing file on stderr and exits 1, so a pipeline reading only
stdout sees success-shaped emptiness.

**And the inverse, which cost more than any of them: a check that FIRED is not evidence anyone
read it.** Every row above is a check that silently did not run. 2026-10-02 produced the mirror
image. `make_membrane_system._warn_phase_composition_mismatch` warns when a leaflet's declared
phase is at odds with its sterol content; on ex17's lower leaflet (43% CHL1, declared `Ld`) it
says *"expect the leaflet to order (drift toward Lo) during equilibration; consider phase: Lo"*.

That is exactly the slow mode that went on to consume a multi-week investigation -- P-P thickness
climbing 33.9 -> 42.3 A and still rising, `corr(S, P-P) = +0.955`, a premature certification, and
a published APL that is an artifact. The warning is in **every** ex17 build log held, across at
least two release cycles (3.22.1 and 3.22.8 replicates, every gate-off and validation run), and
correctly silent for ex16, which is pure DMPC and whose gate works.

It was overridden by a comment sitting beside the setting: *"the cis double bond can't be ordered
by the trans bias (its Lo ensemble == its fluid Ld one), so declare it Ld: same result."* That
premise is false -- SOPE/SOPS are stearoyl and POPE palmitoyl, so sn-1 is fully saturated and only
sn-2 carries the kink -- and measurably so: building the same patch both ways gives 28.083 A (Ld)
against 29.008 A (Lo), identical across seeds. The ensembles differ. Nobody re-derived the
comment, and the tool's own correct answer sat in the logs for two releases.

Two things follow:

- **Before concluding the tool was silent about something, grep the build logs for it.** The
  diagnostic may already be there, right, and specific.
- **A diagnostic that only reaches a log is one `tail` away from invisible.** Put it where the
  symptom will be read. The gap is still open here: the phase/composition mismatch is raised in
  `make_membrane_system`, while the symptom shows up chunks later in `membrane_equilibrate`'s
  `-membrane.dat`, which carries no trace of it.

A comment that justifies overriding a warning is a claim, and claims in comments rot exactly like
claims in docstrings -- with nothing to fail when they stop being true.

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

## Four defects in one test module, and none of them could fail where anyone would look

Found 2026-10-14 in `tests/unit/test_tasks/test_make_membrane_system.py` -- five `slow` +
`needs_tools` tests that are the **only** tests anywhere that build a bilayer from a composition.
They had been broken for a long time and nothing said so, because of the filter intersection
described above. Worth recording as a set, because they are four different shapes of the same
root condition: *a shared thing that one test mutates and the next one inherits.*

| what was shared | how it failed |
| :--- | :--- |
| `cls.common_patch_relaxation_protocols`, a class-level list of dicts | ycleptic's `dwalk` fills in schema defaults **in place** on the dict it is handed, so after the first test the list carries `temperature`, `pressure`, `addl_paramfiles`, `other_parameters` and lower-cased ensembles; the task then restructures it further. Every later test died in `reconfigure_tasks` with `Attribute 'minimize' invalid ... under 'md'`. |
| `cls.controller`, one per class | `Bilayer.__init__` calls `charmmffcontent.provision(resnames=...)`, which **replaces** the exposed residue set rather than adding to it. The first test narrowed the content to its own lipids; the next got `Cannot find residue POPC in CHARMMFFContent`. |
| `z_dist`, a config key | Exists nowhere -- not in `base.yaml`, not in any Python, not in `bilayer_embed.tcl` (whose arg loop ignores it silently). Two tests passed it and died on schema validation; a third passes it to the Tcl script, where it is simply dropped. |
| `patch_nlipids['lower']`, read as if it sized the quilt | It does not. `make_membrane_system` sizes the lower leaflet from the stress-free ratio `n_upper * apl_upper / apl_lower`, and when embedding sizes **both** from the protein footprint. A 49/49 request correctly builds 50 POPE against 49 POPC, and `assertLessEqual(POPE, 49)` correctly failed. |

The last one is the one to think about. The test was not merely stale: it was asserting against a
config key that **has no effect on the thing it names**. `patch_nlipids['lower']` reaches the two
calibration patches and stops there. Tests now assert against `task.outcome['n_upper']` /
`['n_lower']`, which the task records precisely so the decision is not left to a log line -- and
`n_lower` is a derived number that appears in no config, so without the record there was nothing
to compare a built system to.

Two habits this reinforces, both already in this file in other costumes:

- **A fixture that one test can mutate is not a fixture.** `setUpClass` is the natural place to
  put something expensive, and it is the wrong place for anything a test hands to code that
  normalizes its input in place. Which test breaks depends on alphabetical order, so the symptom
  moves when a test is renamed -- two of the four above only surfaced once an *earlier* test
  stopped failing fast.
- **A config key that the code ignores is worse than one it rejects.** `z_dist` on the Tcl side
  is dropped without a word; `patch_nlipids['lower']` is accepted, validated, and then not used
  for the leaflet it names.

**Still open, and not a test bug.** Once the four above were fixed and the module could run end
to end, `test_makemembranesystem_popc_chl1_psm_chl1` failed 2 of 3 runs: patchB's NPT stage dies
with NAMD's `FATAL ERROR: Periodic cell has become too small for original patch grid!`. It builds
PSM/CHL1 50/50 at `SAPL: 75`, far above where a cholesterol-rich leaflet sits, and the patch
collapses faster than the patch grid can follow.

**And pestifer says so, before the build, in every run** -- for both leaflets:

```
WARNING pestifer.tasks.make_membrane_system: lower_leaflet: phase Ld requested but sterol
fraction is 50% (cholesterol-rich) -- expect the leaflet to order (drift toward Lo) during
equilibration; consider phase: Lo
```

That is the ex17 lesson again, in a test rather than a build: the diagnostic was already there,
already right, already specific, and nobody read it -- here because the module it was warning
about could not run at all. Marked `xfail(strict=False)` with the cause in its `reason`, rather
than retuned: choosing a composition-appropriate `SAPL` (or `phase`, or a longer patch protocol)
is a decision about what the test should build, not a defect in the code under test.

## The working tree IS the installed pestifer on this machine

Verified 2026-10-15, after `coord` raised it: a bare `python3` run from anywhere on panacea
resolves `import pestifer` to `/home/cfa/Git/pestifer/pestifer/__init__.py`, through an
editable `.pth` in the repo venv. **There is no install step between an edit and anything that
imports pestifer.** An uncommitted change is live, including a deliberately broken one.

That matters here because of how often this file asks you to break something on purpose.
"Prove the check can fail -- revert the fix and confirm the test goes red" is the standard
move, and doing it in place publishes the broken state to every process on the host for as
long as it takes to run the check. I did exactly that the day this was written: reverted
example 17's relaxation ladder in the tree, ran the tests against it, restored from a copy,
about 30 seconds live.

**Run negative controls in an export, not in the tree.** The clean-export gate above already
builds one; use it for this too:

```bash
export BR=~/devtests/pestifer/breakroom && rm -rf $BR && mkdir -p $BR
git archive HEAD | tar -x -C $BR && cd $BR
# now break whatever you like and run the check here
```

The cost is one `tar`. The alternative is a window in which someone else's build reads your
half-finished idea as the shipped code -- and they will have no way of knowing, because
nothing about their run looks unusual.

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

## Clearing a piercing is itself an unphysical displacement

Found 2026-10-19 on example 17. `ring_check` resolves a pierced ring by rotating a branch about
a rotatable bond, picking the pose with the **fewest** new heavy-atom clashes. Fewest is not
few. patchB's post-minimize check accepted three poses leaving 5, 2 and 3 new clashes; an
earlier pass accepted one leaving **24**. Nothing relaxed them, and the next stage was
`NVT 410 K, 1 fs, rigidbonds on`.

It died with six RATTLE failures, and every one of them is in or beside a rotated branch:

| atom | is | why |
| :--- | :--- | :--- |
| 448 | `C28 SOPS L1-4` | one end of the rotation bond `445-448` |
| 470 | `C216 SOPS L1-4` | same lipid, inside the rotated branch |
| 19143 | `C218 SOPS L1-179` | rotated twice in that pass |
| 350, 7489, 8475 | `SOPE L1-3`, `CHL1 L1-69`, `CHL1 L1-78` | the neighbors they were pushed into |

The guard had been written as a sandwich -- `ring_check, minimize, ring_check` -- with a comment
reasoning that "a rotation is a small rigid move chosen for fewest new clashes; the dynamics
that follows relaxes it". Dynamics does not get the chance: the first step sees the overlap at
full force. **A check that rotated owes a minimize, and that minimize owes a check.**

The hard part is the terminator, and it is not "add one more minimize". Every minimize needs a
check after it, and every check that rotates needs a minimize after it, so the obvious loop does
not end. It ends by making the last pass *unable to move anything*: `ring_check.resolve: false`
detects and reports only, and anything it finds is fatal. That is the right answer for what
survives two rounds -- a doubly-threaded molecule, which one rigid rotation about one hinge
cannot undo, or a composition whose `SAPL` is far from where that leaflet sits.

Two things to keep:

- **A pose chosen as the best available is not thereby acceptable.** `min(new_clashes)` over a
  set of bad options is still a bad option. The eligibility filter added in `b23cd1e0` helps and
  does not finish the job; the minimize is what finishes it.
- **The sandwich's comment was a claim, and it rotted exactly like the ex17 phase comment
  did.** It was written when `delete` still mattered and reasoned about deletions; by the time
  rotation was the main path it was describing a different mechanism. A comment that justifies
  *not* doing something is the kind to re-derive.

`tests/unit/test_tasks/test_ringcheck.py::TestGuardSplicesTheLoop` pins the spliced order and,
separately, the invariant itself -- no dynamics stage may follow a resolving `ring_check`
without a minimize in between -- so a future splice that gets the order right by accident still
has to satisfy the rule. `TestVerifyOnlyPass` pins that `resolve: false` raises and does not
move atoms, with a negative control that the default *does* reach the rotation path. That
control was initially written as `assertRaises(Exception)` and passed against an
`AttributeError` from an incomplete stub -- this file's own lesson, caught by the fact that the
spy never fired.

## A resolver must not be judged only by the detector it is trying to satisfy

Found 2026-10-22, one release after the slide move above, while drawing a docs figure of the
"resolved" ring -- the figure showed the bond still going through it. `PSFRing.pierced_by`
rejected any bond whose midpoint lay more than **half a bond length** from the ring COM, with
the reason string "both atoms on same side of ring plane". It was a distance gate wearing a
side-test's name, and it hid every off-center threading.

On example 17's replica-1 quilt the 3.29.4 resolver reported 6 of 6 piercings cleared; a
gate-free test found **3 still threaded**, each bond straddling its ring with winding exactly
2π, crossing 1-2 Å off-center. The resolver's move search asks `check_coords` "is it still
pierced?" -- the same detector -- so it searched for, and found, poses the detector could not
see. A small slide or spin is precisely the move that pushes a centered crossing off-center.
3.29.4 was yanked.

The repo-wide lesson is the one this file keeps relearning, in its sharpest form: **an optimizer
will find the blind spots of its own objective.** A search that accepts the first pose the
checker passes is a search for the checker's false negatives. The defense is a second,
independent check -- different algorithm, no shared code -- before believing a resolution. Here
that was a segment-plane intersection with an even-odd ray cast, and it is what confirmed the
corrected 6-to-0 rather than the detector confirming itself.

Also worth knowing: the gate stayed hidden for months because build-produced threadings sit
near the ring center, and the minimizer stretches the threaded bond (2.5-2.9 Å measured, against
1.53 Å), which widens a half-bond-length gate. The detector and a gate-free test agree on
pristine build output. They diverge only after something moves.

`tests/unit/test_psfutils/test_pierced_by_geometry.py` pins the geometry. Note
`test_a_bond_wholly_on_one_side_is_not_pierced`: simply deleting the gate would have been wrong
too -- a bond pointing at the hole from above has a line that threads the ring and a segment
that does not, and only a real straddle test tells them apart.

## A move set is a hypothesis about what can go wrong, and it was wrong here

Found 2026-10-22 from the 3.29.3 sweep, which stopped ex17/rep-01 with three piercings the
resolver could not clear. Both of its remedies are ROTATIONS -- bend the piercing tail
(`try_glycan_pendant`), turn the speared ring (`try_sterol_spin`) -- and that encodes an
assumption nobody had stated: that a piercing is a molecule caught at one point.

CHL1 L1-197 was caught at two. Measured on the state the failing check read:

| ring | pierced by | crossing | ring radius |
| :--- | :--- | :--- | :--- |
| A-ring (6-membered, the 3-OH end) | SOPE 198 `C212-C213` | 0.17 A from center | 2.05 A |
| D-ring (5-membered, the tail end) | SOPS 196 `C215-C216` | 0.18 A from center | 1.79 A |

A sterol's long axis runs near both ring centers, so a spin moves neither hole; bending either
acyl chain leaves the other holding the ring. **903 and 917 candidate rotations, not one of
which un-threaded it.** A rigid translation does: 636 of 1022 clear it, 25 with no hard overlap,
the best at 1.5 A for 5 soft contacts.

Three things worth keeping:

- **The docstring said not to.** `try_sterol_spin` asserted that "a sterol displaced laterally
  would trade one piercing for a clash", and that sentence is why no translation move existed.
  It was never measured. This is the third claim-in-a-comment in a month to be load-bearing and
  false (see the ex17 phase comment and the CCD determinism note), and the pattern is specific:
  a comment that justifies NOT doing something is never exercised by any test, so nothing can
  contradict it.

- **"No rotation cleared it" and "nothing can clear it" are different statements**, and the log
  said the first while everyone read the second -- including me, until I enumerated the poses.
  Before concluding a structure is unresolvable, enumerate what was actually tried.

- **A second failure mode hid behind the first.** CHL1 L1-162 was reported identically and was
  not trapped at all: it had two clearing poses against the pristine structure and lost them to
  the rotations applied before it, because each accepted pose may leave soft clashes. Run first
  it clears; run in detection order it has nothing left. The translation fallback happens to fix
  that too, which is why it was preferred over reordering the search -- but the greedy-allocation
  weakness is still there for anything the fallback cannot reach.

`tests/unit/test_psfutils/test_ring_resolve.py::TestPierceeTranslate` pins the move (including
that it refuses any residue with a bond leaving it -- a rigid translation of one residue inside a
covalent chain would be nonsense), and
`tests/unit/test_tasks/test_lipid_piercing_rotation.py::TestTheSlideIsTheLastResort` pins the
dispatch **behaviorally** rather than by source text, because this file's other dispatch tests
read `inspect.getsource` and a pattern matching a word inside a comment is already on the list
above.

## "Deterministic" in a comment can be true of the algorithm and false of the build

Found 2026-10-21, chasing why example 11 pierced a ring in sweep replica 3 and not in 1 or 2.
`close_loops_ccd` carried the note *"Deterministic: per-loop seed = base + index."* That is
true: `close_one_loop` seeds `np.random.default_rng(seed + i)` and is reproducible **given its
input**. The build is not, because its input is not:

| stage | across two runs of the same command |
| :--- | :--- |
| `00-01-000_psfgen-build.pdb` | identical |
| `00-03-000_md-minimize.coor` -- the minimize feeding `ligate` | **differs** |

Confirmed directly rather than inferred: the same `.namd` on byte-identical input, same binary,
`+p24` twice gives two different `.coor`; `+p1` twice gives identical ones. So a single-threaded
run is the only reproducible one, and that is the escape hatch if a loop-closure bug ever has to
be bisected.

Two things follow, and the second is the expensive one:

- **`--seed` does not reach the loop closure at all.** It sets `namd.seed`; `ligate.ccd.seed` is
  a separate key with its own default. All three sweep replicas of example 11 closed loops from
  the identical seed and still diverged, so "replica-specific" there is a draw, not a property
  of the seed. Anything reasoning about replicas as seed variants is wrong about this stage.

- **A rerun is not a control here.** Having fixed the resolver, I reran example 11 at the same
  seed, got a clean structure, and nearly reported it as proof. It was a different loop
  conformation that happened not to pierce. What settled it was a controlled probe on the ONE
  structure -- same geometry, same mover set, the resolver run on each side of `connect`:

      PRE-connect  (loop still a pendant)   movers=332  CLEARED=True   unresolved=0
      POST-connect (joined at both ends)    movers=331  CLEARED=False  unresolved=1

  When the pipeline between you and the thing you changed is not reproducible, an end-to-end
  rerun cannot be the evidence. Isolate the step.

Do not read "deterministic" in this repository as "rerunning reproduces it" without checking
what feeds the step.

## A schema key with no consumer passes every check pestifer has

Found 2026-10-19, asking why "add one more CHARMM stream file" had no obvious keyword.
`charmmff.user_custom.{rtf,str,prm}` had been declared in `base.yaml`, documented with a worked
example in `usage.rst`, given generated pages under `docs/source/config_ref/`, and read by
`util/citations.py` to decide whether to print the CGenFF-program citations. **No scripter ever
read them.** `scripters/psfgen.py` concatenated `standard.rtf + custom.rtf + standard.str +
custom.str`; `scripters/namd.py` concatenated `standard.prm + standard.str` then `custom.prm +
custom.str`. `user_custom` was in neither, so

```yaml
charmmff:
  user_custom:
    str: [myligand.str]      # validates, builds, loads nothing
```

produced a build with the ligand's topology absent. The key even changed the *citations*
without changing the force field.

Nothing could catch it:

- **`yclept check-spec` passes** -- it checks that declarations are well-formed, not that
  anything consumes them.
- **Every example parses**, because no example used the key.
- **The docs build fine**, and generate a page describing behavior that does not exist.
- **The feature appeared to work** for the common case: a file in a `searchpath` directory is
  picked up automatically *by residue name* (`tasks/psfgen.py:732`), so the dead key was
  bypassed by a live mechanism. It only failed for a file whose residues nothing in the input
  names -- an extra stream loaded for its parameters -- which is the case nobody tried.

This is `z_dist` again (see the membrane-test table) at a different scale: there, a key the Tcl
arg loop silently dropped; here, a whole documented extension point. **A config key that the
code ignores is worse than one it rejects**, and the only thing that finds one is a grep for
the key outside the module that declares it:

```bash
grep -rn "user_custom" --include=*.py --include=*.tcl pestifer/ \
  | grep -vE "charmmffcontent.py|resourcemanager.py|citations.py"
```

Run that against a *working* key as a positive control before believing a quiet result -- the
same grep shape finds `['custom']`'s two real consumers. An empty result that you have not
calibrated is the same success-shaped emptiness as `sha256sum` on a missing file.

Fixed by replacing the three provenance-named groups (`standard` = in the CHARMM release,
`custom` = shipped by pestifer, `user_custom` = yours) with one flat `charmmff.files` list plus
`charmmff.replace_defaults`, since what a user needs to express is *intent*, not identity.
`pestifer/charmmff/fileset.py` holds the list, the two derived views, and the legacy fold.

**The list's order is the load order, and that is load-bearing twice.** Topology order fixes
terminal patches (see the `DEFA` section below). Parameter order decides duplicates, because
`CharmmParamFile.merge` is last-wins -- and that is not hypothetical: merging the shipped set in
its real order versus "all .prm, then all .str" changes exactly **one** record out of 46,527,

```
('NG2O1','CG2R61','CG2R61','NG2S3') n=2    Kchi 1.25  ->  3.1
```

which is `toppar_pestifer_conflict_resolutions.prm` overriding `toppar_all36_carb_imlab.str`,
annotated in the file as *"CGenFF v5.0 value, chosen over imlab's analogy value"*. Tidying the
default list by extension silently reverts a deliberate arbitration. Do not reorder
`fileset.DEFAULT_FILES` for readability.

**And the test for this had to be written one level lower than felt natural.** Asserting that
`fileset.topologies()` contains the user's file does not prove the *scripter* emits a `topology`
line for it -- the same aim-one-level-off error as rows 4 and 5 of the table above.
`TestTheScripterActuallyEmitsAUserFile` drives the real scripters from a config dict; reverting
either scripter's `fileset.resolve` call reddens exactly those three tests and leaves the
nineteen pure-`fileset` tests green.

**The sweep is now a test, and it found eight more.**
`tests/unit/test_core/test_no_dead_schema_keys.py` requires every leaf key in `base.yaml` to be
named somewhere in `pestifer/*.py` or `pestifer/*.tcl`. Run against the tree on the day it was
written it flagged two more live defects:

- **`constraints.consfile`, in seven places** -- `md`, `membrane_equilibrate`,
  `density_equilibrate`, `ligate.steer`, `terminate.package.namd`, and both
  `make_membrane_system` relaxation protocols. `BaseTask.make_constraint_pdb` reads **`consref`**
  (`specs.get('consref','')`), which the schema did not declare, so the documented key did
  nothing *and* the working key could not be set. Renamed to `consref`.
- **`make_membrane_system.embed.q_tolerance`** -- "tolerance for net charge neutralization",
  default `1e-4`, in no Python and no Tcl. Neutralization is `autoionize -neutralize`, which has
  no tolerance argument. Removed.

It is a *name* search, so it proves a mention, not a read; and it is blind to a key composed at
run time, which is what `DYNAMIC_KEYS` is for. That allowlist is the part that rots: five
`membrane_equilibrate.area_*` keys are on it because `_params_for` does
`specs.get(f'{prefix}_{key}')`, and a companion test asserts that f-string is still in the file,
so an allowlist entry cannot outlive its reason and quietly turn the sweep into a rubber stamp.

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

## Three bugs from one idea living at several call sites — stop copying, start extracting

Recorded 2026-09-30, after the third one. The pattern is not "someone forgot"; it is that a
correct fix, applied at the sites known at the time, leaves a *shape* that the next call site
reproduces wrongly:

| the idea | sites | how it failed |
| :--- | :--- | :--- |
| guard on the parsed cell, not the xsc path | 3 | fixed at `make_membrane_system`, missed `RingChecker.check` one layer down (2026-08-27/28) |
| upper-case both sides of an atom-type comparison | 8 | `_bond_key`/`_angle_key` fixed, `_dihedral_key`/`_improper_key`/`_nbfix_key` left, `psf_param_check.py` untouched for a month |
| drop a stream file the release does not ship | 3 | `continuation` and `merge` each had it — **their comments even named the same file** — and `make_membrane_system`'s prebuilt branch, added later, did not (v3.24.2) |

The third is the instructive one, because the duplication was *self-documenting* and still did not
propagate: two call sites carried a comment naming `toppar_water_ions_namd.str` explicitly, and the
third was written without either. Reading the other sites is not what makes a new call site
correct; there is no step at which the author of a new branch is made to read them.

**So the fix for the third instance was not a third copy.** `CHARMMFFContent.stage_streamfiles_locally`
now holds that filter and all three call sites call it. Prefer this whenever you find yourself
about to write the second copy of a behavior — and when you find the second copy already there,
treat it as a defect report about the third site you have not looked at yet.

Two caveats learned while doing it:

- **Extracting widens the blast radius.** Pointing `continuation` and `merge` at a shared helper
  touched two *working* paths to fix a third. That is usually right here, but it is a real cost
  and belongs in the commit message, not hidden in a diff.
- **The helper must be opt-in where the judgment differs.** Dropping a topology template is safe
  where the structure comes from the PSF (`readpsf`, continuation, a prebuilt bilayer) and unsafe
  where the template is what builds the residues. So this is a method callers choose, never
  something `copy_charmmfile_local` does on its own — and its docstring says which is which,
  because the next caller has to make that call.

Sweep with `grep -rn stage_streamfiles_locally` and `grep -rn get_toppar_from_psf`; a new consumer
of a PSF's topology remarks that does neither is the fourth instance.

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
