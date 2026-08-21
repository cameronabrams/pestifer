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
