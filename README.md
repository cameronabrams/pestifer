# Pestifer
> NAMD System Preparation Tool

[![PyPI](https://img.shields.io/pypi/v/pestifer.svg)](https://pypi.org/project/pestifer/)
[![PyPI Downloads](https://static.pepy.tech/badge/pestifer)](https://pepy.tech/projects/pestifer)
[![Docs](https://readthedocs.org/projects/pestifer/badge/?version=latest)](https://pestifer.readthedocs.io/en/latest/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.16051498-blue)](https://doi.org/10.5281/zenodo.16051498)

Pestifer is a fully automated simulation-ready MD system preparation tool, requiring as inputs only biomolecular structures (e.g., PDB IDs, PDB files, mmCIF files, alphafold IDs) and a handful of customization parameters, to generate NAMD-compatible input files (PSF, PDB, and xsc).  It is basically a highly functionalized front end for VMD's `psfgen` utility.  It also has a few handy subcommands for working with NAMD output.

## Installation

```bash
pip install pestifer
```

Once installed, the user has access to the main `pestifer` command. 

Pestifer also requires access to the following executables:

1. `namd3` and `charmrun`
2. `vmd` and `catdcd`

Pestifer **includes a mirrored copy of** the [Feb 2026 Charmm36 force field](https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_feb26.tgz), plus a few local corrections to upstream files and some added residue/ligand/ion definitions — see [CHARMM force-field customizations](https://pestifer.readthedocs.io/en/latest/charmmff-customizations.html) for the full list and rationale.

## Documentation

Please visit [readthedocs](https://pestifer.readthedocs.io/en/latest) for full documentation.

## Version History

See the [CHANGELOG](https://github.com/cameronabrams/pestifer/blob/main/CHANGELOG.md) for full details.

## Meta

[https://github.com/cameronabrams](https://github.com/cameronabrams/)

Pestifer is maintained by [Cameron F. Abrams](mailto:cfa22@drexel.edu).

Pestifer is distributed under the MIT license. See ``LICENSE`` for more information.

Pestifer was developed with support from the National Institutes of Health via grants GM100472, AI154071, and AI178833.

## Contributing

Pestifer is developed on GitHub at <https://github.com/cameronabrams/pestifer>.

**Code changes:**

1. Fork, clone, and install editable: `pip install -e .`
2. Branch, make your change, and add tests.
3. Run the suite: `pytest` (or `uv run pytest`).
4. Add a bullet under `## [Unreleased]` in [`CHANGELOG.md`](https://github.com/cameronabrams/pestifer/blob/main/CHANGELOG.md).
5. Push and open a Pull Request.

**Content contributions** — a new example, a PDB-repository entry, or a custom CHARMM residue — use the `modify-package` subcommand, which makes the branch, commits exactly the files it touches, and prints the push / PR steps for you (no manual branching needed). See the [modify-package documentation](https://pestifer.readthedocs.io/en/latest/subs/modify-package.html).

