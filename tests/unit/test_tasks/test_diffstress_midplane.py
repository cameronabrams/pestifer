"""The differential-stress diagnostic splits the leaflets where the bilayer IS, not mid-cell.

Aimed at the call site, ``MakeMembraneSystemTask._leaflet_stress_frames``, on synthetic files of
the shapes the diagnostic pass writes (no NAMD): a production log with step-numbered profile
records, a DCD with the cell, an xst and a PSF.  The bilayer sits off-center -- midplane at slab
3.5 of 10 -- and all the excess tension is in slab 4, just ABOVE the midplane.  Charged to the
leaflets about the measured midplane that is upper-leaflet tension, Dgamma > 0; split at the
central slab, as the diagnostic used to, slab 4 counts as LOWER and Dgamma comes out negative.
That a one-slab misplacement flips the answer is why this matters (example 17: 8-15 mN/m per slab).
"""
import struct

import numpy as np

from pestifer.tasks.make_membrane_system import MakeMembraneSystemTask, _leaflet_stress_summary

NSLABS, C, ORIGIN_Z, CENTER = 10, 60.0, 30.0, 21.0      # slab 0 starts at 0; 6 A slabs
STEPS = [500 * (k + 1) for k in range(24)]


def _rec(b):
    return struct.pack('<i', len(b)) + b + struct.pack('<i', len(b))


def _write(tmp):
    rng = np.random.default_rng(2)
    n = 200
    z = np.r_[CENTER + rng.uniform(2, 18, n), CENTER - rng.uniform(2, 18, n)]
    natoms = len(z)
    with open(tmp / 'p.psf', 'w') as f:
        f.write(f'PSF\n\n{natoms:8d} !NATOM\n')
        for i in range(natoms):
            f.write(f'{i + 1:8d} MEMB {i // 10 + 1:4d} DMPC C21  CTL2  0.000000 12.0110 0\n')
    icntrl = [0] * 20
    icntrl[0], icntrl[1], icntrl[2], icntrl[10], icntrl[19] = len(STEPS), STEPS[0], 500, 1, 24
    out = _rec(b'CORD' + struct.pack('<20i', *icntrl)) + _rec(struct.pack('<i', 1) + b' ' * 80) \
        + _rec(struct.pack('<i', natoms))
    xyz = np.zeros((natoms, 3), dtype='<f4')
    xyz[:, 2] = z
    for _ in STEPS:
        out += _rec(struct.pack('<6d', 50.0, 90.0, 50.0, 90.0, 90.0, C))
        out += b''.join(_rec(xyz[:, d].tobytes()) for d in range(3))
    (tmp / 'p.dcd').write_bytes(out)
    prof = np.zeros((NSLABS, 3))
    prof[4, 2] = 100.0                              # Pzz - PL = 100 bar in slab 4 only
    with open(tmp / 'p.log', 'w') as f:
        for s in [0] + STEPS:
            f.write(f'PRESSUREPROFILE: {s} ' + ' '.join(map(str, prof.ravel())) + '\n')
    with open(tmp / 'p.xst', 'w') as f:
        f.write('#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z\n')
        for s in [0] + STEPS:
            f.write(f'{s} 50 0 0 0 50 0 0 0 {C} 25 25 {ORIGIN_Z}\n')
    (tmp / 'p.namd').write_text('pressureProfile on\n')
    return {k: str(tmp / f'p.{k}') for k in ('namd', 'log', 'dcd', 'xst', 'psf')}


def test_the_split_is_at_the_measured_midplane(tmp_path):
    files = _write(tmp_path)
    task = MakeMembraneSystemTask.__new__(MakeMembraneSystemTask)
    frames = task._leaflet_stress_frames(files, reciprocal=False)
    res = _leaflet_stress_summary(frames['lower'], frames['upper'])
    slab_tension = 0.01 * 100.0 * (C / NSLABS)       # 6 mN/m
    assert abs(res['upper'] - slab_tension) < 1e-6
    assert abs(res['lower']) < 1e-6
    assert res['dgamma'] > 0, 'slab 4 was charged to the lower leaflet: split at mid-cell?'
