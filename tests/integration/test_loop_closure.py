# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Regression tests for loop modeling (the ``ligate`` task).

Residues that are present in a structure's SEQRES but absent from its coordinates -- disordered
loops -- are modelled in by growing them from the resolved residue before the gap, then closed
onto the downstream anchor by cyclic coordinate descent (``pestifer.psfutil.loop_ccd``, driven by
``LigateTask.close_loops_ccd``).

A bad closure is silent.  The loop is *present* either way; it is simply strained, knotted, or
not actually joined to the protein.  Nothing in the build fails, and the defect only surfaces
later as a NAMD blow-up or, worse, as a plausible-looking trajectory of a wrong structure.  So
these tests assert on geometry:

* every modelled residue is actually built;
* the chain is continuous across the closure -- the junction peptide bond is a real bond, not the
  tens-of-angstroms gap of an unclosed loop;
* after the minimization that always follows ``ligate`` in a real config, the system satisfies
  the general structural battery, including no atom overlaps.

The last point is deliberately checked *after* minimization.  CCD closes a loop geometrically and
leaves H-H contacts behind -- verified here at 0.32-0.70 A immediately post-closure, all of them
inside the modelled loops -- which the following minimization relieves.  Asserting on the
unminimized structure would be asserting on an intermediate no user consumes.

4ZMJ (HIV-1 Env SOSIP trimer) is used because it carries three real gaps of very different sizes
(9, 11 and 21 residues), including one spanning insertion codes.  The 21-residue loop is the
interesting case: short loops close almost regardless of method.
"""

import math
import textwrap

import pytest

from pestifer.core.config import Config
from pestifer.core.controller import Controller

from tests.integration.helpers import assert_psf_sane, parse_psf_pdb

pytestmark = pytest.mark.needs_tools

#: Gaps in 4ZMJ that the modeller must fill: (segname, first resid, last resid).
_MODELLED_LOOPS = [('G', 400, 410), ('B', 548, 568)]

#: A peptide C-N bond is ~1.33 A.  An unclosed loop leaves its terminus tens of angstroms from
#: the anchor, so anything under 2 A means the chain is genuinely joined.
_MAX_PEPTIDE_BOND = 2.0


def _build(workdir):
    """fetch -> psfgen -> minimize -> ligate -> minimize, as a real loop-modelling config does."""
    cfg = workdir / 'loop.yaml'
    cfg.write_text(textwrap.dedent("""
        title: loop closure regression (4zmj)
        tasks:
          - fetch:
              sourceID: 4zmj
          - psfgen:
              source:
                biological_assembly: 0
          - md:
              ensemble: minimize
              nsteps: 100
          - ligate:
              method: ccd
              ccd:
                refine: 200
                ensemble: 3
          - md:
              ensemble: minimize
              nsteps: 500
          - terminate:
              basename: ligated
    """).lstrip())
    config = Config(userfile=str(cfg)).configure_new()
    return Controller().configure(config).do_tasks()


@pytest.mark.slow
def test_modelled_loops_are_built_and_closed(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    report = _build(tmp_path)
    failed = {r['taskname']: r['result'] for r in report.values() if r['result'] != 0}
    assert not failed, f'build task(s) failed: {failed}'

    psf, pdb = tmp_path / 'ligated.psf', tmp_path / 'ligated.pdb'
    atoms, bonds, coords = parse_psf_pdb(psf, pdb)

    # (a) the modelled residues exist at all
    for seg, lo, hi in _MODELLED_LOOPS:
        present = {int(a['resid']) for a in atoms.values()
                   if a['segname'] == seg and a['resid'].isdigit() and lo <= int(a['resid']) <= hi}
        assert present == set(range(lo, hi + 1)), (
            f'segment {seg}: modelled loop {lo}-{hi} is incomplete; missing '
            f'{sorted(set(range(lo, hi + 1)) - present)}')

    # (b) the backbone is continuous through each modelled loop and across both junctions.
    #     Walk C(i)->N(i+1) from the residue before the gap to the one after it.
    by_key = {}
    for aid, a in atoms.items():
        if a['resid'].isdigit():
            by_key[(a['segname'], int(a['resid']), a['name'])] = aid
    bonded = {tuple(sorted(b)) for b in bonds}

    for seg, lo, hi in _MODELLED_LOOPS:
        checked = 0
        for r in range(lo - 1, hi + 1):
            c, n = by_key.get((seg, r, 'C')), by_key.get((seg, r + 1, 'N'))
            if c is None or n is None:
                continue
            assert tuple(sorted((c, n))) in bonded, (
                f'segment {seg}: no peptide bond between residue {r} C and {r + 1} N -- '
                f'the modelled loop is not joined to the chain')
            d = math.dist(coords[c], coords[n])
            assert d < _MAX_PEPTIDE_BOND, (
                f'segment {seg}: peptide bond {r}C-{r + 1}N is {d:.1f} A -- loop not closed')
            checked += 1
        # guard: if the residue naming ever changes, this test must not silently check nothing
        assert checked >= (hi - lo), (
            f'segment {seg}: only {checked} backbone bonds checked across loop {lo}-{hi}; '
            f'this test would not detect an unclosed loop')

    # (c) the finished system is structurally coherent
    assert_psf_sane(psf, pdb, context='loop closure')


@pytest.mark.slow
def test_ligate_on_a_gapless_structure_is_a_no_op(tmp_path, monkeypatch):
    """A ligate task on a structure with no gaps must succeed, not fail.

    Regression: the bypass path returned ``None`` rather than 0, and once a failed task became
    fatal, any config carrying a ligate task -- copying an Env example onto a complete structure
    is the obvious way to get one -- aborted the whole build.
    """
    monkeypatch.chdir(tmp_path)
    cfg = tmp_path / 'nogap.yaml'
    cfg.write_text(textwrap.dedent("""
        title: ligate no-op regression (BPTI has no gaps)
        tasks:
          - fetch:
              sourceID: 6pti
          - psfgen:
          - ligate:
              method: ccd
          - terminate:
              basename: nogap
    """).lstrip())
    config = Config(userfile=str(cfg)).configure_new()
    report = Controller().configure(config).do_tasks()

    results = {r['taskname']: r['result'] for r in report.values()}
    assert results.get('ligate') == 0, (
        f'ligate on a gapless structure returned {results.get("ligate")!r}; a structure with no '
        f'loops is a valid input, not a build failure')
    assert not [t for t, r in results.items() if r != 0], f'build task(s) failed: {results}'
