# BPTI build fixture

A single BPTI (6PTI) build, five files, used by the `needs_tools` integration tests in
`test_merge.py` and `test_desolvate.py`.

| file | from | used by |
| :--- | :--- | :--- |
| `00-01-00_psfgen-build.{psf,pdb}` | the vacuum psfgen build; 1108 atoms, segments A/B/C, three `DISU` patches on A | `test_merge.py` |
| `00-04-00_solvate.{psf,pdb}` | the solvated system, same three `DISU` remarks | `test_desolvate.py` |
| `00-05-00_md-minimize.dcd` | one minimization frame of that solvated system | `test_desolvate.py` |

**Why these are tracked rather than produced.** Both test modules used to read them out of
`scratch/builds/1/`, which is in `.gitignore` -- so they were leftovers of a build someone had
once run on that machine. They passed on a developer box and raised `FileNotFoundError` on every
fresh clone; CI never noticed because both modules are `needs_tools` and skip there. That is the
accumulated-untracked-state failure described in the repository's `CLAUDE.md`, and the fix is for
the fixture to be part of the repository.

`bpti1.yaml` here is the exact config that produced them -- an older example 01, which has since
been rewritten to use `density_equilibrate` in place of five explicit NPT stages. It is kept
because the task indices in the file names (`00-01-00`, `00-04-00`, `00-05-00`) are positions in
*that* task list, so regenerating from today's example 01 would not reproduce these names.
Nothing reads any of this from `scratch/` any more.
