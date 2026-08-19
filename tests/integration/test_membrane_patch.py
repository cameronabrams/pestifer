# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Regression test for grid-based bilayer construction.

``make_membrane_system`` is the most intricate code in the package and the least covered by unit
tests.  It has produced silent geometric defects before -- the grid/psfgen split stranding atoms
at the origin, and an orthohexagonal lattice search that returned badly-shaped quilts.  Both
build and run; neither announces itself.

The invariants asserted here are the packer's own claims, stated in
``Bilayer.write_grid_pdb``: lipids are dropped on a per-leaflet lattice, oriented with head
groups pointing *away* from the midplane and tails toward it.  So:

* both leaflets are populated, with the counts the composition asked for;
* every lipid's phosphate sits farther from the midplane than that lipid's own centroid -- i.e.
  no lipid is upside down, which a sign error in the orientation step would produce wholesale;
* the two headgroup planes are separated by a physically sensible distance for the lipid;
* and the general structural battery passes, which is what catches origin-stranded atoms.

Marked ``expensive``: this builds a real bilayer and takes ~5.5 minutes, four times the rest of
the integration suite combined.  It is therefore excluded from the release gate and run
deliberately::

    pytest tests/integration -m expensive --runslow

A deliberately small patch (2x2, minimize-only relaxation) is used.  The point is to exercise
grid packing and the psfgen split, not to equilibrate anything -- the equilibration path is
better tested with synthetic series at unit level than by waiting for real MD.
"""

import collections
import textwrap

import numpy as np
import pytest

from pestifer.core.config import Config
from pestifer.core.controller import Controller

from tests.integration.helpers import assert_psf_sane, parse_psf_pdb

pytestmark = pytest.mark.needs_tools

_LIPID = 'DMPC'

#: DMPC phosphate-to-phosphate bilayer thickness is ~35-40 A.  A generous window still excludes
#: an interdigitated collapse (too small) and two leaflets that never met (too large).
_MIN_PP_SEPARATION, _MAX_PP_SEPARATION = 25.0, 55.0


def _build(workdir):
    cfg = workdir / 'patch.yaml'
    cfg.write_text(textwrap.dedent(f"""
        title: minimal bilayer patch (structure only)
        tasks:
          - make_membrane_system:
              bilayer:
                packer: grid
                SAPL: 60.0
                npatch: [2, 2]
                composition:
                  lower_leaflet:
                    - name: {_LIPID}
                      frac: 1.00
                  upper_leaflet:
                    - name: {_LIPID}
                      frac: 1.00
                seed: 270272
                relaxation_protocols:
                  patch:
                    - md:
                        ensemble: minimize
                        nsteps: 100
                  quilt:
                    - md:
                        ensemble: minimize
                        nsteps: 100
          - terminate:
              basename: patch
    """).lstrip())
    config = Config(userfile=str(cfg)).configure_new()
    return Controller().configure(config).do_tasks()


@pytest.mark.slow
@pytest.mark.expensive
def test_grid_bilayer_is_a_correctly_oriented_bilayer(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    report = _build(tmp_path)
    failed = {r['taskname']: r['result'] for r in report.values() if r['result'] != 0}
    assert not failed, f'build task(s) failed: {failed}'

    psf, pdb = tmp_path / 'patch.psf', tmp_path / 'patch.pdb'
    atoms, _bonds, coords = parse_psf_pdb(psf, pdb)

    # per-lipid centroid and phosphate z
    centroid_z, phosphate_z = collections.defaultdict(list), {}
    for aid, a in atoms.items():
        if a['resname'] != _LIPID:
            continue
        key = (a['segname'], a['resid'])
        centroid_z[key].append(coords[aid][2])
        if a['name'] == 'P':
            phosphate_z[key] = coords[aid][2]

    assert centroid_z, f'no {_LIPID} lipids in the built patch'
    cz = {k: float(np.mean(v)) for k, v in centroid_z.items()}
    assert set(cz) == set(phosphate_z), 'some lipids have no phosphate atom'

    midplane = 0.5 * (min(cz.values()) + max(cz.values()))
    lower = [k for k, z in cz.items() if z < midplane]
    upper = [k for k, z in cz.items() if z >= midplane]

    # (a) both leaflets exist and are equally populated, as a 1.00/1.00 composition asks
    assert lower and upper, (
        f'lipids did not separate into two leaflets about z={midplane:.1f}: '
        f'{len(lower)} below, {len(upper)} above')
    assert len(lower) == len(upper), (
        f'leaflets are unequal: {len(lower)} lower vs {len(upper)} upper, from a composition '
        f'that requested the same species at frac 1.00 in each')

    # (b) no lipid is upside down: its head group is farther from the midplane than its centroid
    inverted = []
    for leaflet, sign in ((lower, -1.0), (upper, +1.0)):
        for k in leaflet:
            if sign * (phosphate_z[k] - cz[k]) <= 0:
                inverted.append((k, round(phosphate_z[k], 1), round(cz[k], 1)))
    assert not inverted, (
        f'{len(inverted)} lipid(s) point the wrong way -- phosphate is nearer the midplane than '
        f'the lipid centroid, e.g. {inverted[:5]} (resid, P z, centroid z)')

    # (c) the two headgroup planes sit a physically sensible distance apart
    sep = float(np.mean([phosphate_z[k] for k in upper]) -
                np.mean([phosphate_z[k] for k in lower]))
    assert _MIN_PP_SEPARATION < sep < _MAX_PP_SEPARATION, (
        f'phosphate-to-phosphate separation is {sep:.1f} A, outside the plausible '
        f'{_MIN_PP_SEPARATION}-{_MAX_PP_SEPARATION} A window for a {_LIPID} bilayer')

    # (d) the general battery -- this is what catches atoms stranded at the origin by the
    #     grid/psfgen split, the defect that motivated testing this path at all
    assert_psf_sane(psf, pdb, unminimized=True, context='membrane patch')
