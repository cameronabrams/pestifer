# pestifer

[![tests](https://github.com/cameronabrams/pestifer/actions/workflows/tests.yaml/badge.svg)](https://github.com/cameronabrams/pestifer/actions/workflows/tests.yaml)
[![PyPI](https://img.shields.io/pypi/v/pestifer.svg)](https://pypi.org/project/pestifer/)
[![Python](https://img.shields.io/pypi/pyversions/pestifer)](https://pypi.org/project/pestifer/)
[![License: MIT](https://img.shields.io/pypi/l/pestifer)](https://github.com/cameronabrams/pestifer/blob/main/LICENSE)
[![Docs](https://readthedocs.org/projects/pestifer/badge/?version=latest)](https://pestifer.readthedocs.io/en/latest/)
[![Downloads](https://static.pepy.tech/badge/pestifer)](https://pepy.tech/projects/pestifer)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.16051498-blue)](https://doi.org/10.5281/zenodo.16051498)

**Pestifer** is a versatile system-preparation tool that facilitates the use of the [VMD](https://www.ks.uiuc.edu/Research/vmd/) tool [psfgen](https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/ug.pdf) to generate CHARMM-force-field compatible PSF and PDB files for use in the MD simulation package [NAMD](https://www.ks.uiuc.edu/Research/namd/).  Pestifer automates and extends the standard `psfgen/VMD` workflow using a simple YAML interface to provide a versatile and user-friendly way to set up molecular dynamics simulations.

Pestifer includes the CHARMM36 force field files from the [MacKerell Lab](https://mackerell.umaryland.edu/charmm_ff.shtml) (February 2026 release).

```{note}
Pestifer is under active development.  Pestifer development is supported in part by Grants GM100472, AI154071, and AI178833 from the NIH.
```

## Key capabilities

```{include} ../../README.md
:start-after: '## Key capabilities'
:end-before: '## Installation'
```

## Contents

```{toctree}
:maxdepth: 2

introduction
installation
usage
examples
config_ref
CHARMM FF customizations <charmmff-customizations>
API <api/API>
Tcl sources <tcl/source>
Changelog <changelog>
Roadmap <roadmap>
```
