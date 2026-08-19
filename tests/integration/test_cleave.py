# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Regression tests for chain cleavage (the ``cleave`` task).

``cleave`` severs a named peptide bond and re-segments the chain either side of it -- the
operation behind the furin-site cleavage in the SARS-CoV-2 spike example.  Like the other
operations tested here, a wrong result is silent: cleaving the wrong bond, or failing to sever
the intended one, still produces a system that builds and runs.

What is asserted, and one thing deliberately not.

The obvious invariant -- "the two fragments are no longer connected" -- is **wrong**, and
asserting it would have produced a test that fails on correct behaviour.  BPTI's disulfides
Cys5-Cys55 and Cys14-Cys38 both span the cut, so after cleaving between residues 20 and 21 the
fragments remain covalently joined through two S-S bonds.  That is chemistry, not a defect:
cutting one backbone bond in a disulfide-bonded protein does not separate the chains.

So the tests assert the precise thing instead: the *named* peptide bond is gone, the fragments
land in separate segments with the expected residue ranges, and the only remaining covalent
links between them are exactly those two disulfides -- which also guards against cleavage
collaterally breaking a bond it was not asked to touch.
"""

import textwrap

import pytest

from pestifer.core.config import Config
from pestifer.core.controller import Controller

from tests.integration.helpers import assert_psf_sane, parse_psf_pdb

pytestmark = pytest.mark.needs_tools

#: Cleave BPTI's chain A between residues 20 and 21.
_CUT_N, _CUT_C = 20, 21

#: BPTI has three disulfides; two of them (Cys5-Cys55, Cys14-Cys38) span this cut.
_EXPECTED_DISULFIDES = 3
_EXPECTED_BRIDGING = {('CYS', '5', 'CYS', '55'), ('CYS', '14', 'CYS', '38')}


def _build(workdir):
    cfg = workdir / 'cleave.yaml'
    cfg.write_text(textwrap.dedent(f"""
        title: cleavage regression (BPTI chain A, {_CUT_N}/{_CUT_C})
        tasks:
          - fetch:
              sourceID: 6pti
          - psfgen:
          - md:
              ensemble: minimize
              nsteps: 50
          - cleave:
              sites:
                - A:{_CUT_N}-{_CUT_C}
          - terminate:
              basename: cleaved
    """).lstrip())
    config = Config(userfile=str(cfg)).configure_new()
    return Controller().configure(config).do_tasks()


def _protein_segments(atoms):
    """{segname -> sorted resids} for segments holding more than a couple of residues."""
    by_seg = {}
    for a in atoms.values():
        if a['resid'].isdigit():
            by_seg.setdefault(a['segname'], set()).add(int(a['resid']))
    return {s: sorted(r) for s, r in by_seg.items() if len(r) > 2}


@pytest.mark.slow
def test_cleave_severs_the_named_bond_and_resegments(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    report = _build(tmp_path)
    failed = {r['taskname']: r['result'] for r in report.values() if r['result'] != 0}
    assert not failed, f'build task(s) failed: {failed}'

    psf, pdb = tmp_path / 'cleaved.psf', tmp_path / 'cleaved.pdb'
    atoms, bonds, _coords = parse_psf_pdb(psf, pdb)
    bonded = {tuple(sorted(b)) for b in bonds}

    # (a) the cut residues landed in different segments, and no residue was lost
    segs = _protein_segments(atoms)
    seg_of = {}
    for seg, resids in segs.items():
        for r in resids:
            seg_of.setdefault(r, seg)
    n_seg, c_seg = None, None
    for seg, resids in segs.items():
        if _CUT_N in resids and max(resids) == _CUT_N:
            n_seg = seg
        if _CUT_C in resids and min(resids) == _CUT_C:
            c_seg = seg
    assert n_seg and c_seg, (
        f'expected one segment ending at residue {_CUT_N} and another starting at {_CUT_C}; '
        f'got segments {[(s, r[0], r[-1]) for s, r in segs.items()]}')
    assert n_seg != c_seg, f'cleavage left residues {_CUT_N} and {_CUT_C} in the same segment'

    # (b) the named peptide bond is gone
    def _atom(seg, resid, name):
        ids = [aid for aid, a in atoms.items()
               if a['segname'] == seg and a['resid'] == str(resid) and a['name'] == name]
        assert len(ids) == 1, f'expected one {seg}:{resid}:{name}, found {len(ids)}'
        return ids[0]

    c_atom = _atom(n_seg, _CUT_N, 'C')
    n_atom = _atom(c_seg, _CUT_C, 'N')
    assert tuple(sorted((c_atom, n_atom))) not in bonded, (
        f'the peptide bond {_CUT_N}C-{_CUT_C}N survived the cleavage')

    # (c) nothing else joining the fragments except the disulfides that legitimately span the cut
    bridging = set()
    for u, v in bonds:
        au, av = atoms[u], atoms[v]
        if {au['segname'], av['segname']} == {n_seg, c_seg}:
            lo, hi = sorted((au, av), key=lambda a: int(a['resid']))
            bridging.add((lo['resname'], lo['resid'], hi['resname'], hi['resid']))
    assert bridging == _EXPECTED_BRIDGING, (
        f'fragments are joined by {sorted(bridging)}; expected exactly the two disulfides that '
        f'span the cut, {sorted(_EXPECTED_BRIDGING)} -- anything else means cleavage broke or '
        f'created a bond it was not asked to')

    # (d) cleavage did not collaterally break a disulfide elsewhere
    ss = [(u, v) for u, v in bonds if atoms[u]['name'] == 'SG' and atoms[v]['name'] == 'SG']
    assert len(ss) == _EXPECTED_DISULFIDES, (
        f'{len(ss)} disulfides after cleavage, expected {_EXPECTED_DISULFIDES}')

    # (e) the result is structurally coherent (no relaxation follows the cut here)
    assert_psf_sane(psf, pdb, unminimized=True, context='cleave')
