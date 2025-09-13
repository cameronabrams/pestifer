# Pestifer
> NAMD System Preparation Tool

[![PyPI](https://img.shields.io/pypi/v/pestifer.svg)](https://pypi.org/project/pestifer/)
[![PyPI Downloads](https://static.pepy.tech/badge/pestifer)](https://pepy.tech/projects/pestifer)
[![Docs](https://readthedocs.org/projects/pestifer/badge/?version=latest)](https://pestifer.readthedocs.io/en/latest/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16051498.svg)](https://doi.org/10.5281/zenodo.16051498)

Pestifer is a fully automated simulation-ready MD system preparation tool, requiring as inputs only biomolecular structures (e.g., PDB IDs, PDB files, mmCIF files, alphafold IDs) and a handful of customization parameters, to generate NAMD-compatible input files (PSF, PDB, and xsc).  It is basically a highly functionalized front end for VMD's `psfgen` utility.  It also has a few handy subcommands for working with NAMD output.

## Installation

```bash
pip install pestifer
```

Once installed, the user has access to the main `pestifer` command. 

Pestifer also requires access to the following executables:

1. `namd3` and `charmrun`
2. `vmd` and `catdcd`
3. `packmol`

Pestifer **includes a copy of** the [July 2024 Charmm36 force field](https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_jul24.tgz).

## Documentation

Please visit [readthedocs](https://pestifer.readthedocs.io/en/latest) for full documentation.

## Version History

See the [CHANGELOG](./CHANGELOG.md) for full details.

## Meta

[https://github.com/cameronabrams](https://github.com/cameronabrams/)

Pestifer is maintained by Cameron F. Abrams.

Pestifer is distributed under the MIT license. See ``LICENSE`` for more information.

Pestifer was developed with support from the National Institutes of Health via grants GM100472, AI154071, and AI178833.

## Contributing

1. Fork it (<https://github.com/cameronabrams/pestifer/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

