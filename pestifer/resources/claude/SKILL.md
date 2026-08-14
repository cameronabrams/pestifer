---
name: pestifer
description: Build NAMD-ready molecular dynamics systems from PDB/AlphaFold structures using pestifer. Use when asked to prepare, build, or solvate an MD system; when working with pestifer YAML configs; or when a task involves psfgen, CHARMM36 topologies, PSF/PDB generation, membrane embedding, or glycan grafting.
---

# pestifer

Pestifer turns a structure (PDB ID, AlphaFold ID, or local file) into a NAMD-ready system --
PSF, PDB, coordinates, periodic cell, and the CHARMM36 parameters it needs -- from a YAML
config that lists tasks to run in order.

## The working loop

1. Get a config (see *Writing a config* below).
2. **`pestifer build config.yaml --check`** -- validates the YAML against the schema, confirms
   `vmd`/`namd3`/`charmrun`/`catdcd` resolve, statically validates the task pipeline, and prints
   the task plan. Takes about a second and touches nothing. Add `--json` for machine-readable
   output. Exit 0 = would build, exit 1 = would not.
3. Fix whatever it reports and re-check. Schema errors name the valid alternatives.
4. Only then run the build -- **in the background** (see *Builds are slow*).

Never skip step 2. A config error that `--check` finds in one second otherwise surfaces
minutes into a run.

## Writing a config

**Prefer adapting a worked example over composing YAML from scratch.** There are 27 of them
covering the common shapes:

```bash
pestifer show-resources examples     # ID, PDB ID, name, and title for each
pestifer fetch-example 7             # copy example 7's config here without building
```

Pick by similarity to the target: simple globular protein (1-6), glycosylated trimer (7-12,
15), membrane-embedded (16, 17), nucleic acid complex (18), AlphaFold input (20), non-aqueous
solvent (23, 24, 26), fusion construct (27).

**To scaffold from an arbitrary structure**, use `new-system`, which fetches the structure,
reads its header, and writes a config annotated with what it found -- biological assemblies,
per-chain identities, missing loops and tails, engineered mutations, expression tags:

```bash
pestifer new-system 4zmj --inspect
```

Findings are emitted as *commented* YAML stubs; uncomment the ones you want.

**Do not use `pestifer new-system --interactive`.** It prompts on stdin and expects a terminal.
Use `--inspect` and edit the result instead.

A minimal config looks like this:

```yaml
title: what this system is
tasks:
- fetch:
    sourceID: 6pti
    source_format: pdb
- psfgen:
    source:
      biological_assembly: 1
- md:
    ensemble: minimize
- solvate:
- md:
    ensemble: minimize
- md:
    ensemble: NVT
    nsteps: 1000
- density_equilibrate:
- terminate:
    basename: my_system
    package:
      basename: prod_system
      namd:
        ensemble: NPT
```

A task with no value under it runs with all defaults.

## Exploring the schema

Do not guess key names. The schema is self-documenting:

```bash
pestifer config-help --no-interactive                    # top-level keys
pestifer config-help tasks psfgen --no-interactive       # what psfgen accepts
pestifer config-help tasks psfgen mods --no-interactive  # and so on, down the tree
```

Omitting `--no-interactive` starts a prompt loop that will hang a non-interactive session.

## Builds are slow -- run them in the background

Preparation is fast; equilibration is not. A 60,000-atom solvated protein takes about
**11 minutes**, of which ~85% is the final `density_equilibrate` task. Membrane systems and
large glycosylated trimers take substantially longer.

Launch detached and poll the log rather than blocking:

```bash
cd /path/to/clean/dir
setsid nohup pestifer build config.yaml > build.log 2>&1 < /dev/null &
```

Then check progress with `tail build.log`; each task logs when it starts and finishes.

## One clean, empty directory per build

A build writes hundreds of intermediate files and reuses predictable names. Always `mkdir` a
fresh directory and run there. `--check` warns when the working directory is not empty.

## When a build fails

Pestifer writes `.pestifer-manifest.json` recording the last cleanly-completed task:

```bash
pestifer build config.yaml --restart      # resume from where it stopped
pestifer build config.yaml --from psfgen  # resume from a specific task
pestifer build config.yaml --fresh        # ignore the manifest, start over
```

Diagnostics are in `<config-stem>-diagnostics.log`. Every intermediate file -- psfgen scripts,
NAMD logs, per-task structures -- is preserved in `<basename>-artifacts.tar.gz`.

## Reading the output

A `terminate` task with a `package:` block sweeps the working directory when it finishes, so
the run-ready files end up **inside a tarball rather than loose**:

```
prod_system.tar.gz          the deliverable
  prod_system/my_system.psf         topology
  prod_system/my_system.pdb         coordinates (text)
  prod_system/my_system.coor        coordinates (NAMD binary)
  prod_system/my_system.vel         velocities
  prod_system/my_system.xsc         periodic cell
  prod_system/my_system_minimal.prm CHARMM parameters this system needs
  prod_system/prod_system.namd      sample NAMD config
my_system-artifacts.tar.gz  every intermediate file
```

That package directory is self-contained and is what gets copied to a production machine.

## Rendering a built system

`scripts/pestifer-snapshot` (source tree only) renders a headless VMD figure. Solvent is drawn
without hydrogens, so `-solvent points` is the mode that reads correctly; `lines` and
`licorice` degenerate to dots because there are no bonds to draw.

## Other useful subcommands

- `pestifer build-example N` -- fetch and build example N in one step
- `pestifer show-resources resname <RESI>` -- look up a residue name pestifer knows
- `pestifer mdplot`, `pestifer density-profile` -- analyze a finished run
- `pestifer --version`, `pestifer <subcommand> --help`

Full reference: https://pestifer.readthedocs.io/en/latest/
