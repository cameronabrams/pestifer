"""Build the lipid gallery: one page per entry of the shipped lipid collection, plus the index.

Every number on a page is read from the entry (its ``info.yaml`` and PSF) or computed from its
conformers by pestifer's own geometry code; nothing is transcribed.  Run it from the repository
root with the docs extras available (RDKit for the 2D drawings, VMD and ImageMagick for the
contact sheets)::

    uv run --extra ligand-paramgen python docs/source/lipid_gallery/build_gallery.py

``--pages-only`` rewrites the ``.rst`` files without rendering, for a change to page text.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import math
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO))

from pestifer.charmmff.make_pdb_collection import acyl_chain_rows, chain_extensions  # noqa: E402

COLLECTION = REPO / 'pestifer' / 'resources' / 'charmmff' / 'feb26' / 'pdbrepository' / 'lipid.tgz'
LO = '__Lo'

# Contact-sheet geometry.  A conformer is rendered straight at tile size; lipidsnap.tcl views it
# down the membrane normal, so on screen Y is the conformer's z and X is x cos15 - y sin15, and
# VMD centers it on the mean of its heavy atoms, the ones drawn (within 0.1 A on POPC's ten
# tiles; the all-atom mean is 1.2 A off, which drove tails into the tile's bottom edge).  Measured on POPC's ten fluid conformers at tile
# height 332: a heavy-atom extent of d A draws as PX_PER_A_PER_SCALE * scale * d pixels, plus
# STICK_PX for the stick radius (both fits within one pixel).
TILE_W, TILE_H, COLUMNS = 188, 332, 5
PX_PER_A_PER_SCALE = 10.6 / 0.0954
STICK_PX = 4.0
FILL = 0.90          # the farthest atom of any conformer of the entry reaches 90% of the tile
TILT = math.radians(15.0)

GROUPS = [
    ('top_all36_lipid.rtf', 'Glycerophospholipids'),
    ('toppar_all36_lipid_miscellaneous.str', 'Further glycerophospholipids'),
    ('toppar_all36_lipid_inositol.str', 'Phosphoinositides'),
    ('toppar_all36_lipid_cardiolipin.str', 'Cardiolipins'),
    ('toppar_all36_lipid_sphingo.str', 'Sphingolipids'),
    ('toppar_all36_lipid_cholesterol.str', 'Sterols'),
    ('toppar_all36_lipid_detergent.str', 'Detergents'),
]

PROVENANCE = [   # (generation key, row label, formatter)
    ('sampler', 'Sampler', lambda v: f'``{v}``'),
    ('phase', 'Phase', lambda v: 'none (fluid)' if v is None else str(v)),
    ('torsion_bias', 'Trans bias', lambda v: f'{v}'),
    ('chain_order', 'Chain order', lambda v: f'{v:.3f}'),
    ('mc_seed', 'MC seed', lambda v: f'{v}'),
    ('n_chains', 'Chains', lambda v: f'{v}'),
    ('cylinder_apl', 'Cylinder APL', lambda v: f'{v} Å²'),
    ('cylinder_inflation', 'Cylinder inflation', lambda v: f'{v}'),
    ('mc_n_equil', 'Equilibration steps', lambda v: f'{v}'),
    ('mc_n_decorr', 'Decorrelation steps', lambda v: f'{v}'),
    ('mc_max_angle', 'Max MC angle', lambda v: f'{math.degrees(v):.0f}°'),
    ('mc_radius_scale', 'Radius scale', lambda v: f'{v}'),
    ('axial_bias', 'Axial bias', lambda v: f'{v}'),
]


def read_pdb(path):
    lines = [ln for ln in open(path) if ln.startswith(('ATOM', 'HETATM'))]
    names = [ln[12:16].strip() for ln in lines]
    xyz = np.array([[float(ln[30:38]), float(ln[38:46]), float(ln[46:54])] for ln in lines])
    return names, xyz


def natom(psf):
    for ln in open(psf):
        if '!NATOM' in ln:
            return int(ln.split()[0])
    raise ValueError(f'{psf}: no !NATOM record')


class Ensemble:
    """One directory of the collection: a base entry, or its ``__Lo`` sibling."""

    def __init__(self, root: Path, name: str, psf: Path):
        self.name, self.dir, self.psf = name, root / name, psf
        self.info = yaml.safe_load(open(self.dir / 'info.yaml'))
        self.conformers = sorted(self.info['conformers'], key=lambda c: c['pdb'])
        self.generation = self.info.get('generation') or {}
        self.rows = []
        for c in self.conformers:
            names, xyz = read_pdb(self.dir / c['pdb'])
            chains = acyl_chain_rows(names)
            ext = min(chain_extensions(xyz, chains, np.array([0.0, 0.0, 1.0]))) if chains else None
            self.rows.append((c['pdb'][:-4], c['head-tail-length'], c['max-internal-length'], ext))

    def half_extents(self):
        """Largest on-screen half-height and half-width of any conformer, from its heavy-atom mean."""
        hh = hw = 0.0
        for c in self.conformers:
            names, xyz = read_pdb(self.dir / c['pdb'])
            heavy = xyz[[not n.startswith('H') for n in names]]
            heavy = heavy - heavy.mean(axis=0)
            hh = max(hh, np.abs(heavy[:, 2]).max())
            hw = max(hw, np.abs(heavy[:, 0] * math.cos(TILT) - heavy[:, 1] * math.sin(TILT)).max())
        return hh, hw


class Entry:
    def __init__(self, root: Path, resname: str):
        self.resname = resname
        self.slug = resname.lower()
        psf = root / resname / f'{resname}-init.psf'
        self.fluid = Ensemble(root, resname, psf)
        self.lo = Ensemble(root, resname + LO, psf) if (root / (resname + LO)).is_dir() else None
        info = self.fluid.info
        self.charge = float(info['charge'])
        self.natom = natom(psf)
        self.topfile = info['defined-in']
        self.synonym = (info.get('synonym') or '').lstrip('! ').strip()

    @property
    def ensembles(self):
        return [self.fluid] + ([self.lo] if self.lo else [])

    def scale(self):
        hh = max(e.half_extents()[0] for e in self.ensembles)
        hw = max(e.half_extents()[1] for e in self.ensembles)
        s_h = (FILL * TILE_H / 2 - STICK_PX) / (PX_PER_A_PER_SCALE * hh)
        s_w = (FILL * TILE_W / 2 - STICK_PX) / (PX_PER_A_PER_SCALE * hw)
        return min(s_h, s_w)


def render_sheet(ens: Ensemble, scale: float, out_png: Path, work: Path):
    """One tile per conformer, in conformer order, tiled five across, saved as a palette PNG."""
    tiles = []
    for c in ens.conformers:
        tile = work / f'{ens.name}-{c["pdb"][:-4]}.png'
        with open(HERE / 'lipidsnap.tcl') as tcl:
            subprocess.run(['vmd', '-dispdev', 'text', '-args', '-psf', str(ens.psf),
                            '-coor', str(ens.dir / c['pdb']), '-o', str(tile),
                            '-size', f'{TILE_W}x{TILE_H}', '-scale', f'{scale:.5f}'],
                           stdin=tcl, capture_output=True, check=True)
        if not tile.exists():
            raise RuntimeError(f'VMD wrote no image for {ens.name} {c["pdb"]}')
        tiles.append(str(tile))
    cols = min(COLUMNS, len(tiles))
    sheet = work / f'{ens.name}-sheet.png'
    subprocess.run(['magick', 'montage', *tiles, '-tile', f'{cols}x{math.ceil(len(tiles) / cols)}',
                    '-geometry', '+0+0', '-background', 'white', str(sheet)], check=True)
    subprocess.run(['magick', str(sheet), '-depth', '8', '-colors', '64',
                    '-define', 'png:compression-level=9', str(out_png)], check=True)
    check_sheet(out_png, len(tiles))


def check_sheet(png: Path, n: int):
    """Every tile drawn and clear of its border: a blank tile is VMD's empty-scene failure (a
    script fed with ``-e`` instead of on stdin), and a clipped one means the scale rule missed."""
    from PIL import Image
    a = np.asarray(Image.open(png).convert('L'))
    cols = min(COLUMNS, n)
    if a.shape != (math.ceil(n / cols) * TILE_H, cols * TILE_W):
        raise RuntimeError(f'{png.name}: sheet is {a.shape[::-1]} for {n} conformers')
    for k in range(n):
        ink = a[(k // cols) * TILE_H:(k // cols + 1) * TILE_H, (k % cols) * TILE_W:(k % cols + 1) * TILE_W] < 245
        if ink.sum() < 200:
            raise RuntimeError(f'{png.name}: tile {k} is blank')
        ys, xs = np.nonzero(ink)
        if min(ys.min(), xs.min(), TILE_H - 1 - ys.max(), TILE_W - 1 - xs.max()) < 2:
            raise RuntimeError(f'{png.name}: tile {k} reaches the edge of its tile')


def table(rows, header):
    widths = [max(len(header[i]), *(len(r[i]) for r in rows)) for i in range(len(header))]
    rule = '  '.join('=' * w for w in widths)
    fmt = lambda r: '  '.join(c.ljust(w) for c, w in zip(r, widths)).rstrip()
    return '\n'.join([rule, fmt(header), rule, *map(fmt, rows), rule])


def conformer_table(ens: Ensemble):
    chains = any(r[3] is not None for r in ens.rows)
    header = ['Conformer', 'Head-tail (Å)', 'Max internal (Å)'] + (['Min chain extent'] if chains else [])
    rows = [[f'``{n}``', f'{ht:.2f}', f'{mi:.2f}'] + ([f'{x:.2f}'] if chains else [])
            for n, ht, mi, x in ens.rows]
    return table(rows, header), chains


def page(e: Entry, has_structure: bool) -> str:
    out = [f'.. _lipid_gallery {e.slug}:', '', e.resname, '=' * len(e.resname), '']
    if e.synonym:
        out += [f'**{e.synonym}**', '']
    nconf = len(e.fluid.conformers)
    conf = f'{nconf} fluid' + (f', {len(e.lo.conformers)} liquid-ordered (``{e.lo.name}``)' if e.lo else '')
    out += [_identity(e, conf), '']
    if has_structure:
        out += ['Structure', '---------', '', f'.. image:: img/{e.slug}-structure.svg', '   :width: 340px',
                f'   :alt: Two-dimensional structure of {e.resname}', '']
    out += ['Conformers', '----------', '']
    if e.fluid.generation.get('sampler') == 'single':
        out += [f'A sterol is dominated by its rigid fused rings, so the collection carries one '
                f'minimized structure for it rather than a sampled ensemble.', '']
    elif e.lo:
        out += [f'The grid packer draws a conformer per lipid from an ensemble, so what the ensemble '
                f'contains is what a built bilayer is made of. {e.resname} ships two: a fluid one, '
                f'used by default, and a liquid-ordered one, used for a leaflet declared '
                f'``phase: Lo``. All conformers on this page are rendered at one shared scale.', '']
    else:
        out += [f'The grid packer draws a conformer per lipid from this ensemble, so what the '
                f'ensemble contains is what a built bilayer is made of. A leaflet declared '
                f'``phase: Lo`` generates an ordered ensemble when a build first needs it. All '
                f'conformers are rendered at one shared scale.', '']
    any_chains = False
    for ens, title, img in ((e.fluid, f'Fluid (``{e.fluid.name}``)', f'{e.slug}-conformers.png'),
                            (e.lo, f'Liquid-ordered (``{e.lo.name}``)' if e.lo else '', f'{e.slug}-lo-conformers.png')):
        if ens is None:
            continue
        if e.lo:
            out += [title, '~' * len(title), '']
        tab, chains = conformer_table(ens)
        any_chains |= chains
        what = 'liquid-ordered ' if ens is e.lo else ''
        out += [f'.. image:: img/{img}', '   :width: 100%' if len(ens.conformers) > 1 else f'   :width: {TILE_W}px',
                f'   :alt: {len(ens.conformers)} {what}{e.resname} conformer(s) at a shared scale', '', tab, '']
    if any_chains:
        out += ['*Min chain extent* is how far the least-extended acyl chain descends below its own '
                'attachment carbon, as a fraction of that chain\'s all-trans contour; an all-trans '
                'chain scores 1.0, and a conformer under 0.25 is rejected at generation as folded.', '']
    out += ['Provenance', '----------', '',
            'How ' + ('each ensemble was' if e.lo else 'the ensemble was') +
            ' generated, read from the entry\'s own ``info.yaml``:', '']
    ens = e.ensembles
    keys = [k for k in PROVENANCE if any(k[0] in x.generation for x in ens)]
    rows = [[label] + [fmt(x.generation[k]) if k in x.generation else '—' for x in ens] for k, label, fmt in keys]
    header = ['Setting'] + [f'``{x.name}``' for x in ens]
    out += [table(rows, header), '']
    return '\n'.join(out)


def _identity(e: Entry, conf: str) -> str:
    rows = [['CHARMM resname', f'``{e.resname}``'], ['Net charge', f'{e.charge:.2f}'],
            ['Atoms', str(e.natom)], ['Topology file', f'``{e.topfile}``'], ['Conformers', conf]]
    w0 = max(len(r[0]) for r in rows)
    w1 = max(len(r[1]) for r in rows)
    rule = f'{"=" * w0}  {"=" * w1}'
    return '\n'.join([rule, *(f'{a.ljust(w0)}  {b}'.rstrip() for a, b in rows), rule])


def index_tables(entries) -> str:
    out = []
    by_file = {}
    for e in entries:
        by_file.setdefault(e.topfile, []).append(e)
    known = [f for f, _ in GROUPS]
    for f, title in GROUPS + [(f, f) for f in sorted(by_file) if f not in known]:
        if f not in by_file:
            continue
        out += [title, '-' * len(title), '', f'Defined in ``{f}``.', '']
        rows = [[f':ref:`{e.resname} <lipid_gallery {e.slug}>`', e.synonym or '—', f'{e.charge:.0f}',
                 'yes' if e.lo else '']
                for e in sorted(by_file[f], key=lambda e: e.resname.lower())]
        out += [table(rows, ['Lipid', 'Name in the topology', 'Charge', 'Lo ensemble']), '']
    return '\n'.join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--collection', default=str(COLLECTION))
    ap.add_argument('--out', default=str(HERE))
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--pages-only', action='store_true')
    ap.add_argument('--only', nargs='*', help='resnames to build (default: every entry)')
    args = ap.parse_args()
    out = Path(args.out)
    (out / 'img').mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='lipid-gallery-'))
    try:
        with tarfile.open(args.collection) as tf:
            tf.extractall(work, filter='data')
        root = next(p for p in work.iterdir() if p.is_dir())
        names = sorted((d.name for d in root.iterdir() if d.is_dir() and LO not in d.name), key=str.lower)
        entries = [Entry(root, n) for n in names]
        build = [e for e in entries if not args.only or e.resname in args.only]
        drawn = {}
        if not args.pages_only:
            import lipid2d
            tiles = work / 'tiles'
            tiles.mkdir()

            def one(e: Entry):
                s = e.scale()
                render_sheet(e.fluid, s, out / 'img' / f'{e.slug}-conformers.png', tiles)
                if e.lo:
                    render_sheet(e.lo, s, out / 'img' / f'{e.slug}-lo-conformers.png', tiles)
                return e.resname, s
            with cf.ThreadPoolExecutor(args.jobs) as pool:
                for name, s in pool.map(one, build):
                    print(f'{name}: scale {s:.4f}', flush=True)
            for e in build:
                svg = out / 'img' / f'{e.slug}-structure.svg'
                try:
                    lipid2d.draw(str(e.fluid.dir / e.fluid.conformers[0]['pdb']), str(svg))
                    svg.with_suffix('.png').unlink(missing_ok=True)   # lipid2d's preview; pages use the SVG
                except Exception as exc:
                    print(f'{e.resname}: no 2D drawing ({type(exc).__name__}: {exc})', flush=True)
                    svg.unlink(missing_ok=True)
        for e in build:
            has = (out / 'img' / f'{e.slug}-structure.svg').exists()
            (out / f'{e.slug}.rst').write_text(page(e, has))
        text = (out / 'index.rst').read_text()
        start, end = '.. gallery-entries-start', '.. gallery-entries-end'
        head, rest = text.split(start, 1)
        tail = rest.split(end, 1)[1]
        toc = ['.. toctree::', '   :hidden:', '', *(f'   {e.slug}' for e in entries)]
        (out / 'index.rst').write_text(f'{head}{start}\n\n' + '\n'.join(toc) + '\n\n'
                                       + index_tables(entries) + f'\n{end}{tail}')
        print(f'{len(build)} page(s) written; index lists {len(entries)}')
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    main()
