# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""How much of a built membrane's area-per-lipid is inherited from the calibration patch?

`make_membrane_system` sizes a leaflet as ``n = target_area / apl``, so where the calibration
patch stops sets the lipid count arithmetically.  Whether that propagates into the finished
system depends on whether the quilt equilibration afterwards relieves a bad count by moving its
own lateral area.  This measures it, stratified by calibration protocol.

**Stratify, and read the stratified numbers, not the pooled one.**  Builds divide into two
protocols -- a gated calibration (patchA CONVERGED, the production path) and a gate-off one
(patchA CEILING) -- which sit at systematically different calibration APLs.  Pooling them makes
a between-protocol difference look like a within-protocol slope.  pestifer-sweep hit exactly
that on an overlapping corpus: a pooled r of +0.948 that became -0.636 and +0.621 inside the two
groups, from two clumps with a gap between them.  This script reports pooled AND per-stratum, and
prints the overlap, so the artifact is visible rather than latent.

What 35 example-17 builds at one protein footprint (13983 A^2) gave, 2026-10-15:

    regression of final quilt APL on calibration APL
      pooled  n=35   +0.858                              <- inflated by stratification
      gated   n=23   +0.599   95% CI [-0.027, +1.226]  p=0.060
      CEILING n=12   +0.845   95% CI [+0.084, +1.606]  p=0.033
    reference values: +1.03 if the quilt is inert and keeps its seed, 0.00 if it finds its own area

    does the quilt's area respond to the lipid count it was handed?
      gated    +10.9 A^2/lipid  p=0.40      full compensation would need ~47
      CEILING   -9.8 A^2/lipid  p=0.55      full compensation would need ~44

So: no evidence the quilt compensates, and **not enough power to exclude that it does** -- the
production stratum's slope CI spans from fully compensating to fully inert.  The claim this
script was written to support ("the quilt keeps 86% of its seed") does not survive stratification
at this sample size, and should be quoted as the per-stratum range or not at all.

What does NOT depend on the mechanism, and is the durable result: across those 35 builds the
final APL ranges **40.42-49.49 A^2** (mean 46.24, sd 2.37) for one composition at one footprint.

Parses build logs only; does not import pestifer, so it runs over an archive from any version.

    python3 -I scripts/measure_apl_seed_propagation.py \
        $(grep -rl "Stress-free leaflet counts" ~/devtests/pestifer --include=*.log)
"""
import re, os, sys
import numpy as np
from scipy import stats as st

CNT  = re.compile(r'Stress-free leaflet counts: upper (\d+), lower (\d+) \(target equilibrium box area (\d+)')
VERD = re.compile(r'membrane_equilibrate-patchA: (CONVERGED at step (\d+)|CEILING)')
PA   = re.compile(r'membrane_equilibrate-patchA: \[stage 2\].*?APL~([0-9.]+)')
QA   = re.compile(r'membrane_equilibrate-quilt: \[stage 2\].*?area=([0-9.e+]+)')


def collect(paths, footprint=13983.0):
    """(calib APL, protocol, n_upper, quilt area, final APL) per build, deduplicated.

    One build often writes both `build.log` and `pestifer_diagnostics.log`; counting it twice
    would halve the standard errors below for free.
    """
    seen, rows = set(), []
    for p in paths:
        try:
            t = open(p, errors='replace').read()
        except OSError:
            continue
        c, v, pa, qa = CNT.search(t), VERD.search(t), PA.findall(t), QA.findall(t)
        if not (c and v and pa and qa):
            continue
        n_u, area, q = int(c.group(1)), float(c.group(3)), float(qa[-1])
        if abs(area - footprint) > 1:      # one footprint only: n and apl are coupled through it
            continue
        key = (round(float(pa[-1]), 2), n_u, q)
        if key in seen:
            continue
        seen.add(key)
        rows.append((float(pa[-1]), 'CEILING' if v.group(1) == 'CEILING' else 'gated',
                     n_u, q, q / n_u))
    return rows


def report(sub, label, footprint=13983.0):
    if len(sub) < 3:
        print(f'{label:18s} n={len(sub):2d}  (too few)')
        return
    c = np.array([r[0] for r in sub]); f = np.array([r[4] for r in sub])
    n = np.array([r[2] for r in sub], float); q = np.array([r[3] for r in sub], float)
    lr, rq = st.linregress(c, f), st.linregress(n, q)
    tcrit = st.t.ppf(0.975, len(c) - 2)
    lo, hi = lr.slope - tcrit * lr.stderr, lr.slope + tcrit * lr.stderr
    print(f'{label:18s} n={len(c):2d}  calibAPL {c.min():.2f}-{c.max():.2f}   '
          f'finalAPL {f.min():.2f}-{f.max():.2f}')
    print(f'{"":18s}   final~calib slope {lr.slope:+.3f}  95% CI [{lo:+.3f},{hi:+.3f}]  '
          f'p={lr.pvalue:.3f}   (inert {q.mean()/footprint:+.2f}, compensating +0.00)')
    print(f'{"":18s}   quilt area ~ n_upper: {rq.slope:+.1f} A^2/lipid  p={rq.pvalue:.2f}   '
          f'(full compensation needs ~{f.mean():.0f})')


if __name__ == '__main__':
    rows = collect(sys.argv[1:])
    if not rows:
        sys.exit('no builds matched')
    report(rows, 'POOLED')
    report([r for r in rows if r[1] == 'CEILING'], 'gate-off/CEILING')
    report([r for r in rows if r[1] == 'gated'], 'gated/CONVERGED')
    ce = sorted(r[0] for r in rows if r[1] == 'CEILING')
    ga = sorted(r[0] for r in rows if r[1] == 'gated')
    if ce and ga:
        lo, hi = max(min(ce), min(ga)), min(max(ce), max(ga))
        if lo <= hi:
            print(f'\ncalibAPL overlap between strata: {lo:.2f}-{hi:.2f} '
                  f'(CEILING n={sum(lo<=x<=hi for x in ce)}, gated n={sum(lo<=x<=hi for x in ga)})')
        else:
            print(f'\nNO OVERLAP between strata ({max(ce):.2f} then {min(ga):.2f}) -- '
                  f'the pooled slope is a between-protocol difference, not a within-protocol one')
