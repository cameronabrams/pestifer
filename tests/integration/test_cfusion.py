# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Regression test for the Cfusion (C-terminal fusion) coordinate orientation.

The fusion-domain orientation was ported from VMD/Tcl to pure numpy
(``pestifer.util.coord.orient_peptide_fusion``, driven by
``PsfgenScripter.write_cfusion_presegment``).  The donor is parsed from its
source PDB, its fusion-domain residues are selected (in serial/sequence order,
so an in-chain modified residue such as a GFP CRO chromophore stays between its
flanking residues), renumbered onto the base chain's C-terminus, and rigid-body
oriented so the donor N-terminus forms a real peptide bond with the base
carbonyl C.

This mirrors example 27 (ubiquitin + GFP): fuse GFP (1ema, chain A, residues
2-229, chromophore CRO included) onto the C-terminus of ubiquitin (1ubq, chain
A).  The test asserts (a) the build succeeds, (b) the fusion junction is a real
peptide bond (~1.33 A, not the tens-of-angstroms gap of unoriented donor
coordinates), and (c) the CRO chromophore survives, in sequence.
"""

import math
import textwrap
from pathlib import Path

import pytest

from pestifer.core.config import Config
from pestifer.core.controller import Controller

from tests.integration.helpers import assert_psf_sane, parse_psf_pdb

pytestmark = pytest.mark.needs_tools

# ubiquitin is 76 residues; GFP residue 2 (1ema) renumbers onto resid 77, so the
# fusion peptide bond joins base resid 76 (C) to donor resid 77 (N).
_BASE_CTERM = 76
_FUSION_NTERM = 77
_MAX_PEPTIDE_BOND = 2.0   # a real C-N peptide bond is ~1.33 A; an unoriented donor gaps tens of A


@pytest.mark.slow
def test_cfusion_junction_and_chromophore(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = tmp_path / 'cfusion.yaml'
    cfg.write_text(textwrap.dedent("""
        title: Cfusion junction regression (ubiquitin + GFP)
        tasks:
          - fetch:
              sourceID: 1ubq
          - psfgen:
              source:
                biological_assembly: 0
              mods:
                Cfusions:
                  - 1ema:A:2-229,A
    """).lstrip())
    config = Config(userfile=str(cfg)).configure_new()
    report = Controller().configure(config).do_tasks()

    failed = {r['taskname']: r['result'] for r in report.values() if r['result'] != 0}
    assert not failed, f'build task(s) failed: {failed}'

    psf = tmp_path / 'my_system.psf'
    pdb = tmp_path / 'my_system.pdb'
    assert psf.exists() and pdb.exists(), 'expected build outputs not found'
    # the general structural battery.  unminimized: this build stops at psfgen, so psfgen's
    # guessed hydrogen positions are still raw -- heavy-atom connectivity is what matters here.
    assert_psf_sane(psf, pdb, unminimized=True, context='cfusion')

    atoms, _bonds, coords = parse_psf_pdb(psf, pdb)

    # locate the junction backbone atoms by (resid, name) in the single protein chain
    def _find(resid, name):
        ids = [aid for aid, a in atoms.items()
               if a['resid'] == str(resid) and a['name'] == name]
        assert len(ids) == 1, f'expected exactly one resid {resid} atom {name}, found {len(ids)}'
        return coords[ids[0]]

    c_base = _find(_BASE_CTERM, 'C')
    n_fusion = _find(_FUSION_NTERM, 'N')
    d = math.dist(c_base, n_fusion)
    assert d < _MAX_PEPTIDE_BOND, (
        f'fusion junction C({_BASE_CTERM})-N({_FUSION_NTERM}) is {d:.2f} A -- '
        f'donor was not oriented onto the base C-terminus')

    # the GFP chromophore (a HETATM-encoded in-chain modified residue) must survive
    resnames = {a['resname'] for a in atoms.values()}
    assert 'CRO' in resnames, 'GFP CRO chromophore missing from the fused construct'
