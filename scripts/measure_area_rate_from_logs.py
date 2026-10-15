# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Recompute `MembraneEquilibrateTask.area_rate` for every stage-2 series in a set of build logs,
grouped by the verdict the convergence gate gave it.

Parses the logs only -- it does not import pestifer -- so it can be run over an archive of builds
from any version.  The statistic is identical to `area_rate`: a linear fit over the late half of
the stage-2 area series, expressed as a percentage of the mean area per 100,000 steps.

    python3 -I scripts/measure_area_rate_from_logs.py $(find ~/devtests/pestifer -name "*build.log")
"""
import re, sys, os
import numpy as np

STAGE2 = re.compile(r'(membrane_equilibrate-(\w+)): \[stage 2\] chunk \d+ @ step (\d+) .*?area=([0-9.e+]+)')
VERDICT = re.compile(r'(membrane_equilibrate-(\w+)): (CONVERGED at step \d+|CEILING[^\n]*)')

series, verdicts = {}, {}
for path in sys.argv[1:]:
    try:
        text = open(path, errors='replace').read()
    except OSError:
        continue
    for m in STAGE2.finditer(text):
        key = (path, m.group(2))
        series.setdefault(key, []).append((float(m.group(3)), float(m.group(4))))
    for m in VERDICT.finditer(text):
        verdicts[(path, m.group(2))] = m.group(3).split(':')[0].split(' at ')[0]

rows = []
for key, pts in series.items():
    if len(pts) < 8:
        continue
    t, a = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
    h = len(t) // 2
    t, a = t[h:], a[h:]                       # late half of stage 2
    if t[-1] == t[0]:
        continue
    slope = np.polyfit(t, a, 1)[0]            # area per step
    rate = 100.0 * slope * 1e5 / a.mean()     # % of mean area per 100k steps
    rows.append((verdicts.get(key, 'none'), abs(rate), key[1], os.path.relpath(key[0])))

for v in ('CONVERGED', 'CEILING'):
    r = sorted(x[1] for x in rows if x[0] == v)
    if r:
        print(f'{v:10s} n={len(r):3d}  median={np.median(r):.3f}  max={max(r):.3f}  min={min(r):.3f}')
print(f'{"(no verdict)":10s} n={len([x for x in rows if x[0]=="none"])}')
print()
for row in sorted(rows, key=lambda x: -x[1])[:6]:
    print(f'  {row[0]:10s} {row[1]:6.3f} %/100k  {row[2]:8s} {row[3][:70]}')

