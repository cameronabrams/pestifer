"""Generate the figures for docs/source/pierced-rings.rst.

Run from the repository root:

    uv run python docs/source/figs/pierced-rings/make_figures.py

Every figure is drawn from data or generated geometry, not by hand.  The real-structure figure
reads ``ex17_chl1_197.json``, extracted from the state example 17 replica 1 stopped on in the
3.29.3 sweep, with the "after" coordinates from the 3.29.5 resolver.  Nothing here needs VMD or
NAMD.

The piercing criterion used for annotation is the one ``PSFRing.pierced_by`` applies since
3.29.5: the bond SEGMENT crosses the ring's mean plane, and it crosses INSIDE the ring.  An
earlier draft of this script tested where the infinite LINE through the bond met the plane, and
called a cleared ring still pierced -- which is how the detector's own blind spot was found.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon
import numpy as np

HERE = Path(__file__).parent
DPI = 200

STEROL = '#4c72b0'
CHAIN_A = '#55a868'
CHAIN_D = '#c44e52'
PIERCE = '#d62728'
CLEAR = '#2ca02c'
RINGFILL = '#9ecae1'
GRAY = '#888888'
LIGHT = '#dddddd'

plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11, 'figure.dpi': DPI,
                     'savefig.bbox': 'tight', 'savefig.pad_inches': 0.08})


# ----------------------------------------------------------------------------------------------
# geometry
# ----------------------------------------------------------------------------------------------

def ring_frame(R):
    """Centroid, in-plane basis (u, v) and unit normal n of a ring's atoms."""
    c = R.mean(axis=0)
    _u, _s, vt = np.linalg.svd(R - c)
    return c, vt[0], vt[1], vt[2]


def point_in_polygon(p, poly):
    """Even-odd ray cast in 2-D."""
    inside = False
    for (x1, y1), (x2, y2) in zip(poly, np.roll(poly, -1, axis=0)):
        if (y1 > p[1]) != (y2 > p[1]) and p[0] < x1 + (p[1] - y1) * (x2 - x1) / (y2 - y1):
            inside = not inside
    return inside


def segment_crossing(R, A, B):
    """Where the bond SEGMENT AB crosses the ring's mean plane, in the ring's 2-D frame.

    Returns (poly, crossing_or_None, pierced).  ``None`` means both bond atoms are on the same
    side of the plane, so the segment does not cross it at all."""
    c, u, v, n = ring_frame(R)
    poly = np.stack([(R - c) @ u, (R - c) @ v], axis=1)
    hA, hB = (A - c) @ n, (B - c) @ n
    if hA * hB >= 0:
        return poly, None, False
    X = A + hA / (hA - hB) * (B - A)
    x2 = np.array([(X - c) @ u, (X - c) @ v])
    return poly, x2, point_in_polygon(x2, poly)


def hexagon(r=1.0, center=(0, 0), rot=0.0, n=6):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False) + rot
    return np.stack([center[0] + r * np.cos(t), center[1] + r * np.sin(t)], axis=1)


# ----------------------------------------------------------------------------------------------
# 1. the real double threading
# ----------------------------------------------------------------------------------------------

def fig_real_case():
    d = json.loads((HERE / 'ex17_chl1_197.json').read_text())
    by = {(a['resid'], a['name']): a for a in d['atoms']}
    serial = {a['serial']: a for a in d['atoms']}
    ster = np.array([a['before'] for a in d['atoms'] if a['resid'] == 197 and a['heavy']])

    cS = ster.mean(axis=0)
    _u, _s, vt = np.linalg.svd(ster - cS)
    long_ax, plate_ax, normal = vt
    tilt = np.deg2rad(55)
    vx, vy = long_ax, np.cos(tilt) * plate_ax + np.sin(tilt) * normal

    def proj(X):
        X = np.atleast_2d(X) - cS
        return np.stack([X @ vx, X @ vy], axis=1)

    fig = plt.figure(figsize=(10.5, 6.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.3, 1], hspace=0.42, wspace=0.18)

    ax = fig.add_subplot(gs[0, :])
    color = {197: STEROL, 198: CHAIN_A, 196: CHAIN_D}
    for i, j in d['bonds']:
        a, b = serial[i], serial[j]
        if a['heavy'] and b['heavy']:
            P = proj(np.array([a['before'], b['before']]))
            ax.plot(P[:, 0], P[:, 1], '-', color=color[a['resid']], lw=1.4,
                    alpha=0.95 if a['resid'] == 197 else 0.5, zorder=2)
    for key in ('A', 'D'):
        R = np.array([by[(197, nm)]['before'] for nm in d['rings'][key]])
        ax.add_patch(Polygon(proj(R), closed=True, fc=RINGFILL, ec=STEROL, lw=1.6, alpha=0.75,
                             zorder=3))
    for pc in d['piercings']:
        A = by[(pc['piercer_resid'], pc['bond'][0])]['before']
        B = by[(pc['piercer_resid'], pc['bond'][1])]['before']
        P = proj(np.array([A, B]))
        ax.plot(P[:, 0], P[:, 1], '-', color=PIERCE, lw=4, solid_capstyle='round', zorder=4)
    centers = {k: proj(np.array([by[(197, nm)]['before'] for nm in d['rings'][k]])).mean(axis=0)
               for k in ('A', 'D')}
    mid = (centers['A'] + centers['D']) / 2
    for key, label in (('A', 'A-ring (3-OH end)'), ('D', 'D-ring (tail end)')):
        p = centers[key]
        out = np.sign(p[0] - mid[0]) or 1.0
        ax.annotate(label, p, xytext=(p[0] + 4.5 * out, p[1] + 5.0), ha='center',
                    arrowprops=dict(arrowstyle='-', color=GRAY, lw=0.8))
    for c, lab, lw in ((STEROL, 'CHL1 197', 2), (CHAIN_A, 'SOPE 198 (threads the A-ring)', 2),
                       (CHAIN_D, 'SOPS 196 (threads the D-ring)', 2), (PIERCE, 'the piercing bonds', 4)):
        ax.plot([], [], '-', color=c, lw=lw, label=lab)
    ax.legend(loc='lower center', ncol=4, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.17))
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Example 17, replica 1: one cholesterol threaded at both ends by two different '
                 'acyl chains (heavy atoms)')

    results = []
    col = 0
    for pc in d['piercings']:
        names = d['rings'][pc['ring']]
        for stage in ('before', 'after'):
            R = np.array([by[(197, nm)][stage] for nm in names])
            A = np.array(by[(pc['piercer_resid'], pc['bond'][0])][stage])
            B = np.array(by[(pc['piercer_resid'], pc['bond'][1])][stage])
            poly, x2, pierced = segment_crossing(R, A, B)
            results.append((pc['ring'], stage, x2, pierced, poly))
            a = fig.add_subplot(gs[1, col])
            col += 1
            a.add_patch(Polygon(poly, closed=True, fc=RINGFILL, ec=STEROL, lw=1.6, alpha=0.8))
            a.plot(*np.vstack([poly, poly[:1]]).T, 'o', color=STEROL, ms=3.5)
            c0 = poly.mean(axis=0)
            for p, nm in zip(poly, names):
                dv = (p - c0) / (np.linalg.norm(p - c0) + 1e-9)
                a.text(*(p + 0.5 * dv), nm, ha='center', va='center', fontsize=6.5, color=GRAY)
            if x2 is None:
                note = 'bond no longer crosses\nthe ring plane'
            else:
                a.plot(*x2, marker='X', ms=8, color=PIERCE if pierced else CLEAR, mec='k',
                       mew=0.6, zorder=5)
                where = 'inside' if pierced else 'outside'
                note = f'bond crosses {np.linalg.norm(x2):.2f} Å\nfrom center, {where} the ring'
            a.text(0.5, -0.05, note, transform=a.transAxes, ha='center', va='top', fontsize=8)
            a.set_title(f"{pc['ring']}-ring, {'as the build stopped' if stage == 'before' else 'after the resolver'}",
                        fontsize=9.5)
            a.set_xlim(-3.4, 3.4)
            a.set_ylim(-3.0, 3.0)
            a.set_aspect('equal')
            a.axis('off')
    fig.savefig(HERE / 'real-double-threading.png')
    plt.close(fig)
    return results


# ----------------------------------------------------------------------------------------------
# 2. the detector's former blind spot
# ----------------------------------------------------------------------------------------------

def fig_gate():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.3), gridspec_kw=dict(width_ratios=[1, 1.15]))

    poly = hexagon(1.40, rot=np.pi / 6)
    a1.add_patch(Polygon(poly, closed=True, fc=RINGFILL, ec=STEROL, lw=1.8, alpha=0.8))
    a1.add_patch(Circle((0, 0), 0.765, fill=False, ls='--', ec='k', lw=1.1))
    a1.text(0, -0.86, 'old gate: ½ bond length', ha='center', va='top', fontsize=8)
    cases = [((0.15, 0.25), PIERCE, 'centered:\ndetected before and after 3.29.5'),
             ((1.0, -0.25), PIERCE, 'off-center, inside the ring:\nMISSED before 3.29.5, detected now'),
             ((2.1, 0.75), CLEAR, 'outside the ring:\nnot a piercing')]
    for (x, y), c, lab in cases:
        a1.plot(x, y, marker='X', ms=13, color=c, mec='k', mew=0.6, zorder=5)
    a1.annotate(cases[0][2], cases[0][0], xytext=(-2.9, 2.0), fontsize=8,
                arrowprops=dict(arrowstyle='-', color=GRAY, lw=0.8))
    a1.annotate(cases[1][2], cases[1][0], xytext=(0.2, -2.6), fontsize=8,
                arrowprops=dict(arrowstyle='-', color=GRAY, lw=0.8))
    a1.annotate(cases[2][2], cases[2][0], xytext=(1.0, 2.0), fontsize=8,
                arrowprops=dict(arrowstyle='-', color=GRAY, lw=0.8))
    a1.set_xlim(-3.1, 3.4)
    a1.set_ylim(-3.0, 2.8)
    a1.set_aspect('equal')
    a1.axis('off')
    a1.set_title('(a) Where the bond crosses the ring plane (face-on)')

    # (b) side view: what the straddle test adds
    a2.plot([-1.4, 1.4], [0, 0], '-', color=STEROL, lw=5, solid_capstyle='round')
    a2.text(1.55, 0, 'ring\n(edge-on)', va='center', fontsize=8, color=STEROL)
    a2.plot([-0.35, -0.35], [-0.8, 0.8], '-', color=PIERCE, lw=4, solid_capstyle='round')
    a2.text(-0.35, -1.0, 'straddles the plane:\npierced', ha='center', va='top', fontsize=8)
    a2.plot([0.55, 0.55], [0.35, 1.9], '-', color=CLEAR, lw=4, solid_capstyle='round')
    a2.plot([0.55, 0.55], [-1.3, 0.35], ':', color=GRAY, lw=1.2)
    a2.text(0.75, 1.55, 'its LINE would pass\nthrough the hole, but the\nbond itself never reaches\nthe plane: not pierced',
            fontsize=8, va='center')
    a2.set_xlim(-2.0, 3.6)
    a2.set_ylim(-2.0, 2.3)
    a2.set_aspect('equal')
    a2.axis('off')
    a2.set_title('(b) Why the replacement is a straddle test (side view)')
    fig.savefig(HERE / 'detector-gate.png')
    plt.close(fig)


# ----------------------------------------------------------------------------------------------
# 3. the three moves
# ----------------------------------------------------------------------------------------------

def _plate(ax, center=(0, 0), alpha=0.85, ls='-', fc=RINGFILL):
    """A cartoon sterol ring plate: three hexagons and a pentagon, fused in a row."""
    cx, cy = center
    rings = [hexagon(0.6, (cx - 1.56, cy), np.pi / 6), hexagon(0.6, (cx - 0.52, cy), np.pi / 6),
             hexagon(0.6, (cx + 0.52, cy), np.pi / 6), hexagon(0.52, (cx + 1.46, cy), 0, n=5)]
    for r in rings:
        ax.add_patch(Polygon(r, closed=True, fc=fc, ec=STEROL, lw=1.4, alpha=alpha, ls=ls))
    return rings


def fig_moves():
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.9))

    # (a) rotate the piercing tail about an upstream bridge
    a = axs[0]
    a.add_patch(Polygon(hexagon(0.8, (0, 0), np.pi / 6), closed=True, fc=RINGFILL, ec=STEROL, lw=1.6))
    chain = np.array([[-2.6, 1.6], [-1.9, 1.2], [-1.2, 1.5], [-0.5, 1.0], [0.0, 0.0],
                      [0.3, -1.0], [1.0, -1.4], [1.7, -1.1]])
    a.plot(chain[:, 0], chain[:, 1], '-o', color=CHAIN_A, ms=3, lw=1.8)
    a.plot(chain[3:5, 0], chain[3:5, 1], '-', color=PIERCE, lw=4)
    hinge = chain[1:3]
    a.plot(hinge[:, 0], hinge[:, 1], '-', color='k', lw=3)
    a.text(-1.55, 1.75, 'hinge\n(a bridge)', ha='center', fontsize=8)
    th = np.deg2rad(-75)
    piv = chain[2]
    rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    moved = (chain[2:] - piv) @ rot.T + piv
    a.plot(moved[:, 0], moved[:, 1], '--o', color=CHAIN_A, ms=3, lw=1.4, alpha=0.6)
    a.add_patch(FancyArrowPatch((0.0, 0.5), (-1.5, -0.9), connectionstyle='arc3,rad=-0.4',
                                arrowstyle='-|>', mutation_scale=14, color=GRAY))
    a.set_title('(a) Bend the piercing tail\nabout an upstream bond')

    # (b) spin the pierced molecule about its long axis
    a = axs[1]
    _plate(a)
    a.plot([-2.6, 2.4], [0, 0], '--', color='k', lw=1)
    a.text(2.45, 0.08, 'long\naxis', fontsize=8, va='bottom')
    a.plot(-1.56, 0.05, marker='X', ms=11, color=PIERCE, mec='k', mew=0.6)
    a.add_patch(FancyArrowPatch((-0.3, 1.2), (0.6, 1.2), connectionstyle='arc3,rad=0.9',
                                arrowstyle='-|>', mutation_scale=14, color=GRAY))
    a.text(0, -1.4, 'turns the ring atoms; a crossing near\nthe axis barely moves relative to them',
           ha='center', fontsize=8)
    a.set_title('(b) Spin the pierced molecule\nabout its own long axis')

    # (c) slide the pierced molecule
    a = axs[2]
    _plate(a, (0, 0), alpha=0.35, ls='--', fc=LIGHT)
    _plate(a, (0.0, -1.45))
    for x in (-1.56, 1.46):
        a.plot(x, 0.02, marker='X', ms=11, color=CLEAR, mec='k', mew=0.6, zorder=5)
    a.add_patch(FancyArrowPatch((2.3, 0.0), (2.3, -1.45), arrowstyle='-|>', mutation_scale=14,
                                color=GRAY))
    a.text(0, 0.95, 'two crossings, at opposite ends:\nno rotation frees both, one slide does',
           ha='center', fontsize=8)
    a.set_title('(c) Slide the whole molecule\n(last resort; free molecules only)')

    for a in axs:
        a.set_xlim(-2.9, 2.9)
        a.set_ylim(-2.3, 2.3)
        a.set_aspect('equal')
        a.axis('off')
    fig.savefig(HERE / 'resolution-moves.png')
    plt.close(fig)


# ----------------------------------------------------------------------------------------------
# 4. where the checks sit around a minimization
# ----------------------------------------------------------------------------------------------

def fig_loop():
    steps = [('a displacement', LIGHT),
             ('ring_check', RINGFILL), ('minimize', '#f2f2f2'), ('ring_check', RINGFILL),
             ('minimize', '#f2f2f2'), ('ring_check\nresolve: false', '#fde0dd'),
             ('dynamics', LIGHT)]
    fig, ax = plt.subplots(figsize=(11.5, 2.6))
    w, h, gap = 1.35, 0.95, 0.32
    xs = []
    for i, (lab, fc) in enumerate(steps):
        x = i * (w + gap)
        xs.append(x)
        ax.add_patch(FancyBboxPatch((x, 0), w, h, boxstyle='round,pad=0.04', fc=fc, ec='k', lw=1))
        ax.text(x + w / 2, h / 2, lab, ha='center', va='center', fontsize=8.5)
        if i:
            ax.add_patch(FancyArrowPatch((xs[i - 1] + w + 0.03, h / 2), (x - 0.03, h / 2),
                                         arrowstyle='-|>', mutation_scale=12, color='k'))
    notes = {0: 'placement, embed,\ngraft, loop closure', 1: 'clear what the\ndisplacement made', 3: 'did the minimize\ncomplete a\nnear-threading?',
             4: 'relax whatever that\ncheck moved', 5: 'detect only:\nanything left is fatal'}
    for i, txt in notes.items():
        ax.text(xs[i] + w / 2, -0.12, txt, ha='center', va='top', fontsize=7.5, color='#444444')
    ax.set_xlim(-0.2, xs[-1] + w + 0.2)
    ax.set_ylim(-0.75, h + 0.15)
    ax.axis('off')
    fig.savefig(HERE / 'resolution-loop.png')
    plt.close(fig)


# ----------------------------------------------------------------------------------------------
# 5. ligate: resolve before connect
# ----------------------------------------------------------------------------------------------

def fig_ligate():
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.7))
    t = np.linspace(0, 1, 9)
    loop = np.stack([-1.0 + 2.0 * t, 1.6 * np.sin(np.pi * t) + 0.0], axis=1)
    for k, a in enumerate(axs):
        a.plot([-3.2, -1.0], [0, 0], '-', color=GRAY, lw=6, solid_capstyle='round')
        a.plot([1.25, 3.2], [0, 0], '-', color=GRAY, lw=6, solid_capstyle='round')
        a.text(-2.1, -0.35, 'resolved chain', ha='center', va='top', fontsize=8)
        a.text(2.2, -0.35, 'downstream anchor', ha='center', va='top', fontsize=8)
        a.add_patch(Polygon(hexagon(0.42, (0.15, 1.15), np.pi / 6), closed=True, fc=RINGFILL,
                            ec=STEROL, lw=1.4))
        if k == 0:
            a.plot(loop[:, 0], loop[:, 1], '-o', color=CHAIN_A, ms=3, lw=2)
            a.plot(loop[-1, 0], loop[-1, 1], 'o', ms=8, mfc='white', mec=CHAIN_A, mew=1.5)
            a.text(1.0, 0.25, 'free end', fontsize=8, ha='left')
            a.plot(loop[2:4, 0], loop[2:4, 1], '-', color='k', lw=3)
            a.add_patch(FancyArrowPatch((0.9, 1.5), (0.2, 2.4), connectionstyle='arc3,rad=0.4',
                                        arrowstyle='-|>', mutation_scale=13, color=GRAY))
            a.set_title('(a) Before connect: the loop is a pendant.\n'
                        'A hinge swings only the loop.  CLEARED.')
        else:
            a.plot(loop[:, 0], loop[:, 1], '-o', color=CHAIN_A, ms=3, lw=2)
            a.plot([loop[-1, 0], 1.25], [loop[-1, 1], 0], '-', color=CHAIN_A, lw=2)
            a.plot(loop[2:4, 0], loop[2:4, 1], '-', color='k', lw=3)
            a.add_patch(Polygon([[-0.45, 0.9], [3.3, 0.9], [3.3, -0.9], [-0.45, -0.9]],
                                closed=True, fc='#fde0dd', ec='none', alpha=0.6, zorder=0))
            a.text(1.4, -1.05, 'a hinge here would have to move\neverything downstream, so none qualifies',
                   ha='center', va='top', fontsize=8)
            a.set_title('(b) After connect: joined at both ends.\n'
                        'No hinge has a small moving side.  NOT cleared.')
        a.set_xlim(-3.4, 3.5)
        a.set_ylim(-2.1, 2.7)
        a.set_aspect('equal')
        a.axis('off')
    fig.savefig(HERE / 'ligate-before-connect.png')
    plt.close(fig)


if __name__ == '__main__':
    for ring, stage, x2, pierced, poly in fig_real_case():
        where = 'no crossing' if x2 is None else f'{np.linalg.norm(x2):.2f} A from center'
        print(f'{ring}-ring {stage:6s}: {where:24s} -> {"PIERCED" if pierced else "clear"}')
    fig_gate()
    fig_moves()
    fig_loop()
    fig_ligate()
    print('wrote', sorted(p.name for p in HERE.glob('*.png')))
