# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""How much of a built membrane's area-per-lipid is inherited from the calibration patch?

`make_membrane_system` sizes a leaflet as ``n = target_area / apl``, so where the calibration
patch stops sets the lipid count arithmetically.  Whether that matters depends on something this
script measures rather than assumes: whether the quilt equilibration afterwards relieves a bad
count by moving its own lateral area, or simply keeps what it was given.

Measured over 35 example-17 builds (2026-10-15), it keeps it -- the regression of the final quilt
APL on the calibration APL has slope +0.858, and the quilt's area barely responds to the lipid
count it was handed (r = +0.195).  So the reported APL is largely a stopping artifact, and the
spread across those builds is 40.42-49.49 A^2.

Parses build logs only; does not import pestifer, so it runs over an archive from any version.
Compares builds at ONE target area (the same protein footprint), since n and apl are coupled
through it.

    python3 -I scripts/measure_apl_seed_propagation.py \
        $(grep -rl "Stress-free leaflet counts" ~/devtests/pestifer --include=*.log)
"""
import re, os, sys
import numpy as np
CNT  = re.compile(r'Stress-free leaflet counts: upper (\d+), lower (\d+) \(target equilibrium box area (\d+)')
VERD = re.compile(r'membrane_equilibrate-patchA: (CONVERGED at step (\d+)|CEILING)')
PA   = re.compile(r'membrane_equilibrate-patchA: \[stage 2\].*?APL~([0-9.]+)')
QA   = re.compile(r'membrane_equilibrate-quilt: \[stage 2\].*?area=([0-9.e+]+)')
seen, rows = set(), []
for p in sys.argv[1:]:
    t = open(p, errors='replace').read()
    c, v, pa, qa = CNT.search(t), VERD.search(t), PA.findall(t), QA.findall(t)
    if not (c and v and pa and qa):
        continue
    n_u, area, q = int(c.group(1)), float(c.group(3)), float(qa[-1])
    if abs(area - 13983) > 1:            # ex17's protein footprint only
        continue
    key = (round(float(pa[-1]), 2), n_u, q)
    if key in seen:                       # the same build logged twice (build.log + diagnostics)
        continue
    seen.add(key)
    rows.append((float(pa[-1]), 'CEILING' if v.group(1) == 'CEILING' else 'converged', n_u, q, q/n_u))
print(f'{"calib APL":>9} {"stopped":>9} {"n_upper":>7} {"quilt area":>10} {"FINAL APL":>9} {"final-calib":>11}')
for apl, verdict, n_u, q, fin in sorted(rows):
    print(f'{apl:9.2f} {verdict:>9} {n_u:7d} {q:10.0f} {fin:9.2f} {fin-apl:+11.2f}')
f = np.array([r[4] for r in rows])
print(f'\nFINAL APL across {len(f)} ex17 builds: mean {f.mean():.2f}  sd {f.std(ddof=1):.2f}  '
      f'range {f.min():.2f}-{f.max():.2f}')
c = np.array([r[0] for r in rows])
print(f'calibration APL          : mean {c.mean():.2f}  sd {c.std(ddof=1):.2f}  range {c.min():.2f}-{c.max():.2f}')
print(f'corr(calibration, final) = {np.corrcoef(c, f)[0,1]:+.3f}')

n = np.array([r[2] for r in rows], dtype=float)
q = np.array([r[3] for r in rows], dtype=float)
print()
print('--- does the quilt adjust its area to the lipid count it was given? ---')
print(f'corr(n_upper, quilt_area) = {np.corrcoef(n, q)[0,1]:+.3f}   '
      f'(+1 = fully compensating -> quilt owns the APL; 0 = no response -> APL is the seed)')
print(f'quilt_area: mean {q.mean():.0f}  sd {q.std(ddof=1):.0f}  ({100*q.std(ddof=1)/q.mean():.1f}% of mean)')
print(f'n_upper   : mean {n.mean():.0f}  sd {n.std(ddof=1):.0f}  ({100*n.std(ddof=1)/n.mean():.1f}% of mean)')
slope = np.polyfit(c, f, 1)[0]
print(f'\nregression final_APL on calibration_APL: slope {slope:+.3f} '
      f'(1.0 = the quilt keeps whatever it was seeded with; 0.0 = it finds its own)')
