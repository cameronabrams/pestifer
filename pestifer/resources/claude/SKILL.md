---
name: pestifer
description: Build NAMD-ready molecular dynamics systems from PDB/AlphaFold structures using pestifer. Use when asked to prepare, build, or solvate an MD system; when working with pestifer YAML configs; or when a task involves psfgen, CHARMM36 topologies, PSF/PDB generation, post-translational modifications, membrane embedding, or glycans.
---

# pestifer

Pestifer turns a structure (PDB ID, AlphaFold ID, or local file) into a NAMD-ready system --
PSF, PDB, coordinates, periodic cell, and the CHARMM36 parameters it needs -- from a YAML
config that lists tasks to run in order.

This skill is bundled with pestifer.  `pestifer setup-claude --force` reinstalls it after an
upgrade; `pestifer --version` says which release you are driving.

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

**Prefer adapting a worked example over composing YAML from scratch.** There are 31 of them.
Name one by its shortname (a numeric id also works, but carries no meaning):

```bash
pestifer show-resources examples            # id, PDB ID, shortname and title for each
pestifer fetch-example hiv-protease         # copy that config here without building
```

Pick by similarity to the target:

- simple globular protein: `bpti1`-`bpti4`, `hiv-protease`, `green-mamba-toxin`
- cofactors, ligands, metals: `sperm-whale-myoglobin` (heme), `ferredoxin-fad`,
  `groel-groes-adp`, `insulin-hexamer` (zinc)
- one copy out of a multi-copy asymmetric unit: `vansc-catalytic-domain`
- AlphaFold input: `methylmalonyl-coa-mutase`
- post-translational modifications: `phosphoubiquitin` (installing one), `myristoylated-matrix`
  (one already in the deposit)
- N-glycosylated trimers, loops, cleavage: `hiv-sosip-env-ectodomain1`..`4`,
  `hiv-ad8-env-ectodomain`, `hiv-ae2-env-ectodomain`, `sars-cov2-S-BA2`
- O-glycans: `notch2-egf-o-glycans`
- membrane-embedded: `hiv-mpertm3-membrane1` (one lipid), `hiv-mpertm3-membrane2` (asymmetric)
- protein-DNA: `ecoli-polymerase`
- non-aqueous solvent: `subtilisin-dmso`, `subtilisin-acetone`, `subtilisin-acetonitrile`
- fusion construct: `ubiquitin-gfp-fusion`
- multi-script build, completing a ligand from a second structure: `hiv-env-cd4-17b-liganded`

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

A task with no value under it runs with all defaults.  The tasks available are `fetch`,
`continuation` (start from an existing PSF/PDB instead of a structure), `psfgen`, `ligate`,
`cleave`, `merge`, `manipulate`, `pdb2pqr`, `solvate`, `desolvate`, `make_membrane_system`,
`ring_check`, `md`, `density_equilibrate`, `membrane_equilibrate`, `validate`, `mdplot` and
`terminate`.

## Exploring the schema

Do not guess key names. The schema is self-documenting:

```bash
pestifer config-help tasks --no-interactive              # the task list
pestifer config-help tasks psfgen --no-interactive       # what psfgen accepts
pestifer config-help tasks psfgen mods --no-interactive  # and so on, down the tree
```

Omitting `--no-interactive` starts a prompt loop that will hang a non-interactive session.

## Builds are slow -- run them in the background

Preparation is fast; equilibration is not. Wall times from the 3.22.1 example sweep, on 24 cores:

- 14,000-60,000-atom solvated proteins: 5-14 minutes, most of it `density_equilibrate`
- 250,000-320,000-atom glycosylated trimers and large complexes: 30-50 minutes
- 500,000-660,000 atoms (glycosylated spike, GroEL/GroES): 2-4 hours
- membrane systems: 8-14 hours

Launch detached and poll the log rather than blocking:

```bash
cd /path/to/clean/dir
setsid nohup pestifer build config.yaml > build.log 2>&1 < /dev/null &
```

Then check progress with `tail build.log`; each task logs when it starts and finishes.
`--ncpus N` sets the NAMD core count (default: auto-detect), `--gpu` forces GPU mode, and
`--seed S` sets the NAMD random seed -- give each replica a different one.

## One clean, empty directory per build

A build writes hundreds of intermediate files and reuses predictable names. Always `mkdir` a
fresh directory and run there. `--check` warns when the working directory is not empty.

## When a build fails

Pestifer writes `.pestifer-manifest.json` recording the last cleanly-completed task:

```bash
pestifer build config.yaml --restart      # resume from where it stopped
pestifer build config.yaml --from psfgen  # resume from a specific task (name or index)
pestifer build config.yaml --fresh        # ignore the manifest, start over
```

Diagnostics are in `<config-stem>-diagnostics.log`. Every intermediate file -- psfgen scripts,
NAMD logs, per-task structures -- is preserved in `<basename>-artifacts.tar.gz`.

Two failures stop the build early on purpose, and their messages say what to do:

- **A residue pestifer does not know** (a crystallization additive such as `EDO`, or an
  unparameterized ligand). The error lists each one and prints the `psfgen: source: exclude:`
  lines that drop it. To keep a ligand instead, `pestifer make-ligand-mol2 input.pdb` prepares
  it for CGenFF.
- **Missing force-field parameters.** Before NAMD runs, pestifer checks that the parameter file
  covers every atom type, bond, angle, dihedral, improper and CMAP term in the PSF, and lists all
  that are missing. Usually the fix is to add the stream file that carries them to
  `charmmff: standard: str:` -- a residue's parameters are not always in the file that defines it.

## Checking what you built

Add a `validate` task to test the result rather than assume it. **Select by molecule, never by
chain or segment letter**: letters in a built system are assigned by pestifer, and excluding a
chain frees its letter for the next derived segment (a water or ion segment can inherit it).
Select `protein`, `resname MG`, `glycan`, and so on.

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
run-record.json             what ran, with which versions and seeds
```

That package directory is self-contained and is what gets copied to a production machine.
`pestifer report-methods <run dirs>` drafts a Methods section and bibliography from the
`run-record.json` files.

## Other useful subcommands

- `pestifer build-example <shortname>` -- fetch and build an example in one step
- `pestifer show-resources resname <RESI>` -- look up a residue name pestifer knows
- `pestifer mdplot`, `pestifer density-profile`, `pestifer pressure-profile-ewald` -- analyze a
  finished run
- `pestifer cache status` -- inspect the per-user caches (`clear`, `rebuild`)
- `pestifer --version`, `pestifer <subcommand> --help`

Full reference: https://pestifer.readthedocs.io/en/latest/
