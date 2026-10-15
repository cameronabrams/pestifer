# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Does the calibration protocol decide the finished membrane's area-per-lipid?

`make_membrane_system` sizes a leaflet as ``n = target_area / apl``, so where the calibration
patch stops sets the lipid count arithmetically.  This measures what reaches the finished system,
over a set of build logs, split by calibration protocol.

**Read the reported APL out of the log; do not derive it.**  pestifer logs
``quilt: <L> lipids/leaflet after relaxation (APL <A>)``, where ``<L>`` is the mean of the two
leaflets' POST-relaxation counts -- not the seeded ``n_upper``, and not the seeded mean.  The
leaflets end asymmetric (seeded ``n_upper - n_lower`` averages +28 here, and the ending counts
differ again from the seeded ones), so ``quilt_area / n_upper`` misses the logged value by a mean
of 3.00 and up to 3.74 A^2 -- about the size of the entire effect anyone is arguing about.  An
earlier version of this script used that derived form and every number it printed was wrong.

What 13 example-17 builds at one protein footprint (13983 A^2) gave, 2026-10-15:

    protocol            n   calibration APL    FINAL APL (logged)
    gate-off/CEILING    6     41.25-42.19       46.20 +- 0.34   [45.83-46.63]
    gated/CONVERGED     7     44.47-46.96       49.24 +- 0.70   [47.95-49.96]
                                    difference  +3.04 A^2, Welch t=10.16, p<0.0001

So the finished membrane is **reproducible within a protocol** -- sd 0.34 and 0.70 A^2 across
independent builds -- and the two protocols give two tightly clustered, cleanly separated
answers 3 A^2 apart.  The variation is not noise; it is the stopping rule.

**Two things this does NOT establish.**  First, which answer is right: that needs an APL for this
composition from outside the pipeline.  Second, causation -- the two groups have NO overlap in
calibration APL (42.19 then 44.47), and the gate-off set was a separate deliberate experiment, so
"protocol" carries every other difference between those runs with it.  A within-protocol
regression cannot separate them either: the slopes are -0.534 (CEILING, p=0.29) and +0.423
(gated, p=0.18) -- opposite signs, neither significant.  That sign flip across a gap is the same
artifact pestifer-sweep documented on 2026-10-15; it is why this script prints the gap.

Parses build logs only; does not import pestifer, so it runs over an archive from any version.

    python3 -I scripts/measure_apl_seed_propagation.py \
        $(grep -rl "Stress-free leaflet counts" ~/devtests/pestifer --include=*.log)
"""
import re, sys
import numpy as np
from scipy import stats as st

CNT  = re.compile(r'Stress-free leaflet counts: upper (\d+), lower (\d+) \(target equilibrium box area (\d+)')
VERD = re.compile(r'membrane_equilibrate-patchA: (CONVERGED at step (\d+)|CEILING)')
PA   = re.compile(r'membrane_equilibrate-patchA: \[stage 2\].*?APL~([0-9.]+)')
QLOG = re.compile(r'quilt: ([0-9.]+) lipids/leaflet after relaxation \(APL ([0-9.]+)')


def collect(paths, footprint=13983.0):
    """One record per build: calibration APL, protocol, seeded counts, and the LOGGED quilt APL.

    Deduplicated -- a build commonly writes both `build.log` and `pestifer_diagnostics.log`, and
    counting it twice would shrink every standard error below for free.
    """
    seen, rows = set(), []
    for p in paths:
        try:
            t = open(p, errors='replace').read()
        except OSError:
            continue
        c, v, pa, ql = CNT.search(t), VERD.search(t), PA.findall(t), QLOG.findall(t)
        if not (c and v and pa and ql):
            continue
        n_u, n_l, area = int(c.group(1)), int(c.group(2)), float(c.group(3))
        if abs(area - footprint) > 1:     # one footprint only: n and apl are coupled through it
            continue
        apl = float(ql[-1][1])
        key = (round(float(pa[-1]), 2), n_u, apl)
        if key in seen:
            continue
        seen.add(key)
        rows.append(dict(calib=float(pa[-1]), n_u=n_u, n_l=n_l, apl=apl,
                         proto='CEILING' if v.group(1) == 'CEILING' else 'gated'))
    return rows


def main(paths):
    rows = collect(paths)
    if not rows:
        sys.exit('no builds matched')
    print(f'{len(rows)} builds with a logged quilt APL at this footprint\n')
    print(f'{"calib":>6} {"protocol":>9} {"seed n_u/n_l":>13} {"FINAL APL":>10}')
    for r in sorted(rows, key=lambda x: x['calib']):
        print(f'{r["calib"]:6.2f} {r["proto"]:>9} {r["n_u"]:6d}/{r["n_l"]:<6d} {r["apl"]:10.2f}')
    groups = {k: np.array([r['apl'] for r in rows if r['proto'] == k]) for k in ('CEILING', 'gated')}
    calibs = {k: np.array([r['calib'] for r in rows if r['proto'] == k]) for k in ('CEILING', 'gated')}
    print()
    for k in ('CEILING', 'gated'):
        g, c = groups[k], calibs[k]
        if g.size < 2:
            print(f'{k:8s} n={g.size} too few'); continue
        print(f'{k:8s} n={g.size:2d}  calib {c.min():.2f}-{c.max():.2f}  '
              f'final APL {g.mean():.2f} +- {g.std(ddof=1):.2f}  [{g.min():.2f}-{g.max():.2f}]')
    if min(g.size for g in groups.values()) >= 2:
        t = st.ttest_ind(groups['gated'], groups['CEILING'], equal_var=False)
        print(f'difference {groups["gated"].mean()-groups["CEILING"].mean():+.2f} A^2  '
              f'Welch t={t.statistic:.2f}  p={t.pvalue:.4f}')
        lo, hi = calibs['CEILING'].max(), calibs['gated'].min()
        if lo < hi:
            print(f'\ncalibration APL GAP between protocols: {lo:.2f} .. {hi:.2f} -- no overlap, so '
                  f'this is a group difference under confounding, NOT a causal slope')
    for k in ('CEILING', 'gated'):
        c, g = calibs[k], groups[k]
        if g.size >= 3:
            lr = st.linregress(c, g)
            print(f'  within {k:8s}: slope {lr.slope:+.3f}  p={lr.pvalue:.3f}')


if __name__ == '__main__':
    main(sys.argv[1:])
