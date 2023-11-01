# Pestifer
> Automatic NAMD System Input Builder

Pestifer is a Python utility for generating atomistic inputs with appropriate topology and parameter files required for molecular dynamics simulations using NAMD.  It is intended as a fully automated system builder requiring as inputs only the PDB codes of biomolecular structures and a handful of customization parameters.  Pestifer primarily works as a front end for VMD's `psfgen` utility.

## Installation

```bash
pip install pestifer
```

Once installed, the user has access to the main `pestifer` command.

## Release History
* 1.0.7
    * `alpha` crotation for folding a span of residues into an alpha helix
* 1.0.6
    * `cif_residue_map_file` generated to report mapping between CIF-residue numbering and author residue numbering
* 1.0.5
    * enhancements to packaging task
* 1.0.4
    * support for topogromacs added
* 1.0.1
    * Initial version

## Meta

Cameron F. Abrams – cfa22@drexel.edu

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/cameronabrams](https://github.com/cameronabrams/)

[https://github.com/AbramsGroup](https://github.com/AbramsGroup/)

## Contributing

1. Fork it (<https://github.com/AbramsGroup/HTPolyNet/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

