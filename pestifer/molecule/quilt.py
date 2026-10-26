# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Build a full membrane by tiling relaxed patches and cropping them to a target box.

Each leaflet of the quilt is taken from a relaxed, periodic patch: tiled on that patch's own
lateral cell, and cropped to the target box by each lipid's headgroup reference atom.  Nothing is
placed, so the quilt starts at the patches' packing.  See ``docs/design/patch-quilting.md``.

Two kinds of contact in the result were never made by any relaxed structure, and they are the
only places a tiled quilt can carry a threaded ring or a near-coincident atom pair:

- the **crop seams**, where a lipid near one crop edge meets the periodic image of a lipid near
  the opposite edge (absent when the box holds a whole number of tiles in that direction);
- the **midplane of an asymmetric build**, where two patches' leaflets meet for the first time.

:func:`find_conflicts` scans the whole assembled quilt for both, rather than trying to enumerate
the seam pairs: every in-tile contact is one the patch already had, so a clean patch yields
findings only at the seams, and a finding anywhere else is a report about the patch.
:func:`resolve_conflicts` removes lipids to clear them -- never moves them.

The module is numpy over arrays read from the patch's PSF and coordinates; the only file it
writes is a lipids-only PDB in the same layout the grid packer writes, so the downstream psfgen
and solvation path is unchanged.
"""
import logging
import os

from dataclasses import dataclass, field

import networkx as nx
import numpy as np

from scipy.spatial import cKDTree

logger = logging.getLogger(__name__)

# Reference atoms, in order of preference.  A phospholipid's phosphorus, else a sterol's 3-OH
# oxygen; anything carrying neither falls back to its centroid.  Never the center of mass for
# leaflet membership: a sterol that has sunk toward the midplane is still in the leaflet its
# hydroxyl points out of.
REFERENCE_ATOM_NAMES = ('P', 'O3')

# A lipid whose reference atom sits this close to the midplane is reported, not silently
# assigned: it is mid flip-flop, and the side it lands on is arbitrary.
MIDPLANE_AMBIGUITY = 2.0

# Atoms of two different lipid copies closer than this are a clash, and one of the two is removed.
# The same floor the grid packer uses (`corruption_floor` in Bilayer.write_grid_pdb): below it,
# coincident atoms can wreck the first minimize; above it, an overlap is ordinary work for that
# minimize, as it is for every gridded membrane.  Measured on an ex17 patch tiled 2.5 x 2.5 with
# a 2 A gap, 0.8 A removed 47 lipids where 0.3 A removed 20 -- every extra removal a seam lipid
# that the minimize would have moved a fraction of an angstrom.
CLASH_DISTANCE = 0.3

# Candidate radius for a piercing test, bond midpoint to ring centroid; PSFRing.pierced_by's own
# gate is 3.5 A, so any value at or above it finds the same piercings.
PIERCE_SEARCH = 4.0

# Rings larger than this are not checked (the same bound RingChecker uses).
MAX_RING_SIZE = 7


@dataclass
class PatchLipid:
    """One lipid of a relaxed patch: its atoms, as rows into the patch arrays."""
    rows: np.ndarray            # atom rows in the patch arrays
    ref: int                    # row of the reference atom, or -1 to use the centroid
    resname: str
    leaflet: str = ''           # 'upper' | 'lower'


@dataclass
class Patch:
    """The lipids of one relaxed, periodic patch, wrapped and registered.

    After :meth:`register` every lipid is whole, its lateral centroid lies in ``[0, ax) x [0, ay)``,
    and the bilayer midplane is at ``z = 0``.
    """
    names: np.ndarray           # per-atom atom names
    resnames: np.ndarray        # per-atom residue names
    xyz: np.ndarray             # (N, 3)
    bonds: np.ndarray           # (M, 2) atom rows, intra-lipid
    lipids: list
    cell: np.ndarray            # (3,) orthorhombic edges
    origin: np.ndarray          # (3,) cell center, NAMD convention
    rings: dict = field(default_factory=dict)   # lipid index -> list of row arrays
    lipid_bonds: dict = field(default_factory=dict)   # lipid index -> (k, 2) rows
    midplane_ambiguous: list = field(default_factory=list)

    def ref_point(self, i):
        lip = self.lipids[i]
        if lip.ref >= 0:
            return self.xyz[lip.ref]
        return self.xyz[lip.rows].mean(axis=0)

    def centroid(self, i):
        return self.xyz[self.lipids[i].rows].mean(axis=0)

    def density(self):
        """Lipids per square angstrom of a leaflet of this patch, taken as the MEAN of its two
        leaflets' counts: a symmetric patch is tensionless as a whole, not per leaflet, so a
        101/99 patch samples one preferred area per lipid (the convention
        ``make_membrane_system`` calibrates by)."""
        n = 0.5 * len(self.lipids)
        return n / (self.cell[0] * self.cell[1])

    def register(self):
        """Wrap each lipid whole into the cell by its lateral centroid, put the cell's lower-left
        corner at the xy origin, assign leaflets, and shift the midplane to ``z = 0``.

        Wrapped by the centroid because the crop is by the centroid: wrapped by any other point,
        a lipid whose centroid falls just outside the cell loses a copy at the edge of every
        whole-tile box (a 2 x 2 box kept 380 of 400 per leaflet that way).
        """
        lo = self.origin - 0.5 * self.cell
        for i, lip in enumerate(self.lipids):
            r = self.centroid(i)
            shift = np.zeros(3)
            shift[:2] = -lo[:2] - self.cell[:2] * np.floor((r[:2] - lo[:2]) / self.cell[:2])
            self.xyz[lip.rows] += shift
        refs = np.array([self.ref_point(i) for i in range(len(self.lipids))])
        # midplane: start from the mean reference z, split, then take the mean of the two
        # leaflets' means -- a leaflet with more sterols (refs sitting deeper) must not drag it
        z0 = refs[:, 2].mean()
        up = refs[:, 2] > z0
        if up.all() or not up.any():
            raise ValueError('patch reference atoms do not split into two leaflets')
        zm = 0.5 * (refs[up, 2].mean() + refs[~up, 2].mean())
        self.xyz[:, 2] -= zm
        refs[:, 2] -= zm
        self.midplane_ambiguous = []
        for i, lip in enumerate(self.lipids):
            lip.leaflet = 'upper' if refs[i, 2] > 0.0 else 'lower'
            if abs(refs[i, 2]) < MIDPLANE_AMBIGUITY:
                self.midplane_ambiguous.append((i, lip.resname, float(refs[i, 2])))
        self._index_topology()
        return self

    def _index_topology(self):
        """Per-lipid bond lists and rings, from the bond graph."""
        owner = np.full(len(self.xyz), -1, dtype=int)
        for i, lip in enumerate(self.lipids):
            owner[lip.rows] = i
        self.lipid_bonds = {i: [] for i in range(len(self.lipids))}
        for a, b in self.bonds:
            if owner[a] >= 0 and owner[a] == owner[b]:
                self.lipid_bonds[owner[a]].append((a, b))
        self.lipid_bonds = {i: np.array(v, dtype=int).reshape(-1, 2)
                            for i, v in self.lipid_bonds.items()}
        # rings are a property of the residue type, so find them once per resname and map by
        # atom name; a residue whose atom names do not match its type's template gets its own
        template = {}
        self.rings = {}
        for i, lip in enumerate(self.lipids):
            key = (lip.resname, tuple(self.names[lip.rows]))
            if key not in template:
                local = {r: k for k, r in enumerate(lip.rows)}
                G = nx.Graph()
                G.add_edges_from((local[a], local[b]) for a, b in self.lipid_bonds[i])
                template[key] = [c for c in nx.minimum_cycle_basis(G) if len(c) <= MAX_RING_SIZE]
                # minimum_cycle_basis returns node SETS; order each ring around its cycle so
                # the winding test sees a polygon, not a scatter
                template[key] = [_ordered_cycle(G, c) for c in template[key]]
            self.rings[i] = [lip.rows[np.array(c)] for c in template[key]]


def _ordered_cycle(G, nodes):
    """The nodes of a ring, in order around it.  A ring from a minimum cycle basis has no
    chords (a chord would close a smaller ring), so its induced subgraph is the cycle itself."""
    return [u for u, _v in nx.find_cycle(G.subgraph(nodes))]


def read_patch(psf, coor, xsc, lipid_resnames=None):
    """Read a relaxed patch's lipids from its PSF, coordinates (PDB or NAMD ``.coor``), and xsc.

    ``lipid_resnames`` overrides which residues count as lipid; by default it is everything
    :func:`~pestifer.util.densityprofile.classify_species` calls a lipid.
    """
    from ..util.density_convergence import _parse_psf_atoms, _read_coor_xyz
    from ..util.densityprofile import classify_species
    from ..util.util import cell_from_xsc
    seg, resid, resn, name, _mass = _parse_psf_atoms(psf)
    xyz = _read_coor_xyz(coor, len(name))
    box, origin = cell_from_xsc(xsc)
    if box is None:
        raise ValueError(f'{xsc}: no periodic cell; a patch to tile must be periodic')
    if lipid_resnames is None:
        is_lipid = classify_species(resn) == 'lipid'
    else:
        is_lipid = np.isin(resn, list(lipid_resnames))
    rows = np.flatnonzero(is_lipid)
    if rows.size == 0:
        raise ValueError(f'{psf}: no lipid atoms')
    keep = {int(r): k for k, r in enumerate(rows)}
    bonds = _read_psf_bonds(psf)
    sel = np.array([a in keep and b in keep for a, b in bonds], dtype=bool) if len(bonds) else np.zeros(0, bool)
    bonds = np.array([(keep[a], keep[b]) for a, b in bonds[sel]], dtype=int).reshape(-1, 2)
    names, resnames, xyz = name[rows], resn[rows], xyz[rows].astype(float)
    lipids = []
    key = [f'{s}:{r}' for s, r in zip(seg[rows], resid[rows])]
    # residues in file order
    first = {}
    for k, kk in enumerate(key):
        first.setdefault(kk, []).append(k)
    for kk, ks in first.items():
        ks = np.array(ks, dtype=int)
        ref = -1
        for nm in REFERENCE_ATOM_NAMES:
            hit = ks[names[ks] == nm]
            if hit.size:
                ref = int(hit[0])
                break
        lipids.append(PatchLipid(rows=ks, ref=ref, resname=str(resnames[ks[0]])))
    return Patch(names=names, resnames=resnames, xyz=xyz, bonds=bonds, lipids=lipids,
                 cell=np.diagonal(box).astype(float).copy(), origin=np.asarray(origin, float).copy())


def _read_psf_bonds(psf):
    """``(M, 2)`` zero-based atom rows of every bond in a PSF."""
    with open(psf) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines) and '!NBOND' not in lines[i]:
        i += 1
    if i == len(lines):
        return np.zeros((0, 2), dtype=int)
    nbond = int(lines[i].split()[0])
    vals = []
    j = i + 1
    while len(vals) < 2 * nbond and j < len(lines):
        vals.extend(int(v) for v in lines[j].split())
        j += 1
    return (np.array(vals[:2 * nbond], dtype=int).reshape(-1, 2) - 1)


@dataclass
class Copy:
    """One tiled copy of a patch lipid."""
    patch: Patch
    lipid: int                  # index into patch.lipids
    shift: np.ndarray           # (3,) translation applied to the patch coordinates
    leaflet: str

    @property
    def xyz(self):
        return self.patch.xyz[self.patch.lipids[self.lipid].rows] + self.shift

    @property
    def centroid(self):
        return self.patch.centroid(self.lipid) + self.shift


def _is_whole(L, a, tol=1e-3):
    n = L / a
    return abs(n - round(n)) < tol and round(n) >= 1


def seams_of(patches, Lx, Ly):
    """``(x, y)``: whether the box has a crop seam in each direction for ANY of ``patches``.

    A direction has no seam only if the box holds a whole number of every patch's tiles there.
    The gap is then applied to every leaflet in a direction that has a seam for any of them, so
    that both leaflets lose the same fraction of their area to it -- a gap paid by one leaflet
    only is a built-in error in the leaflet ratio.
    """
    return (not all(_is_whole(Lx, p.cell[0]) for p in patches),
            not all(_is_whole(Ly, p.cell[1]) for p in patches))


def tile_leaflet(patch: Patch, leaflet: str, Lx: float, Ly: float, gap: float, seams):
    """Tile one leaflet of a registered patch over ``[0, Lx) x [0, Ly)`` and crop it.

    A lipid is kept if its lateral **centroid** lies in ``[g/2, L - g/2)`` in each direction,
    with ``g = gap`` in a direction that has a seam and 0 in one that does not.  The centroid,
    not the headgroup: a lipid's body reaches 5-8 A sideways from its phosphorus, so a
    headgroup crop leaves bodies overlapping across the seam even with a gap.  Measured on an
    ex17 patch, a 2.5 x 2.5-tile box cropped by phosphorus at a 4 A gap left 21 clashes and 19
    threaded rings; cropped by centroid, 7 and 2.
    """
    ax, ay = patch.cell[0], patch.cell[1]
    gx = gap if seams[0] else 0.0
    gy = gap if seams[1] else 0.0
    nx_, ny_ = int(np.ceil(Lx / ax - 1e-9)), int(np.ceil(Ly / ay - 1e-9))
    members = [i for i, lip in enumerate(patch.lipids) if lip.leaflet == leaflet]
    cen = np.array([patch.centroid(i) for i in members]).reshape(-1, 3)
    out = []
    for ix in range(nx_):
        for iy in range(ny_):
            s = np.array([ix * ax, iy * ay, 0.0])
            p = cen + s
            ok = ((p[:, 0] >= gx / 2) & (p[:, 0] < Lx - gx / 2) &
                  (p[:, 1] >= gy / 2) & (p[:, 1] < Ly - gy / 2))
            out.extend(Copy(patch, members[k], s, leaflet) for k in np.flatnonzero(ok))
    return out


def _edge_distance(copy: Copy, Lx, Ly, seams):
    """How far a copy's centroid is from the nearest box edge that is a seam -- or from any box
    edge when there are no seams, since the box edge is then where a trim does least harm."""
    r = copy.centroid
    dx, dy = min(r[0], Lx - r[0]), min(r[1], Ly - r[1])
    if not (seams[0] or seams[1]):
        return min(dx, dy)
    return min(dx if seams[0] else np.inf, dy if seams[1] else np.inf)


def _pierced(ring_xyz, bond_xyz):
    """The same straddle + winding test as :meth:`PSFRing.pierced_by`, on raw arrays."""
    from ..psfutil.psfbond import PSFBond
    from ..psfutil.psfring import PSFRing
    ring = PSFRing(list(range(len(ring_xyz))))
    ring.P = np.asarray(ring_xyz, float)
    ring.calculate_stuff()
    bond = PSFBond([0, 1])
    bond.P = np.asarray(bond_xyz, float)
    bond.calculate_stuff()
    return bool(ring.pierced_by(bond)['pierced'])


@dataclass
class Conflict:
    kind: str                   # 'pierce' | 'clash'
    a: int                      # copy index: the piercer, or one side of a clash
    b: int                      # copy index: the piercee, or the other side
    detail: str = ''


def find_conflicts(copies, Lx, Ly, clash_distance=CLASH_DISTANCE, pierce_search=PIERCE_SEARCH):
    """Every threaded ring and every near-coincident atom pair between two different copies.

    Periodic in x and y over the target box; not periodic in z (the quilt holds lipids only).
    """
    if not copies:
        return []
    xyz_parts, owner_parts = [], []
    for ci, c in enumerate(copies):
        x = c.xyz
        xyz_parts.append(x)
        owner_parts.append(np.full(len(x), ci, dtype=int))
    xyz = np.concatenate(xyz_parts)
    owner = np.concatenate(owner_parts)
    zlo = xyz[:, 2].min() - 1.0
    zspan = xyz[:, 2].max() - zlo + 10.0 * max(clash_distance, pierce_search)
    edges = np.array([Lx, Ly, zspan])

    def boxed(p):
        q = np.array(p, float, copy=True)
        q[..., 2] -= zlo
        q[..., 0] %= Lx
        q[..., 1] %= Ly
        return q

    conflicts = []
    # clashes
    tree = cKDTree(boxed(xyz), boxsize=edges)
    for i, j in tree.query_pairs(clash_distance, output_type='ndarray'):
        if owner[i] != owner[j]:
            conflicts.append(Conflict('clash', int(owner[i]), int(owner[j]),
                                      f'{np.linalg.norm(_mic(xyz[j] - xyz[i], Lx, Ly)):.2f} A'))
    # piercings: every ring against every bond of another copy whose midpoint is near it
    offsets = np.cumsum([0] + [len(p) for p in xyz_parts])
    bond_rows, bond_owner = [], []
    ring_rows, ring_owner = [], []
    for ci, c in enumerate(copies):
        lip = c.patch.lipids[c.lipid]
        local = {r: k for k, r in enumerate(lip.rows)}
        b = c.patch.lipid_bonds.get(c.lipid, np.zeros((0, 2), int))
        if len(b):
            bond_rows.append(np.array([[local[x] + offsets[ci], local[y] + offsets[ci]] for x, y in b]))
            bond_owner.append(np.full(len(b), ci))
        for ring in c.patch.rings.get(c.lipid, []):
            ring_rows.append(np.array([local[r] + offsets[ci] for r in ring]))
            ring_owner.append(ci)
    if ring_rows and bond_rows:
        bond_rows = np.concatenate(bond_rows)
        bond_owner = np.concatenate(bond_owner)
        mid = xyz[bond_rows].mean(axis=1)
        btree = cKDTree(boxed(mid), boxsize=edges)
        coms = np.array([xyz[r].mean(axis=0) for r in ring_rows])
        for k, hits in enumerate(btree.query_ball_point(boxed(coms), pierce_search)):
            rc = ring_owner[k]
            rxyz = xyz[ring_rows[k]]
            for bi in hits:
                if bond_owner[bi] == rc:
                    continue
                bxyz = xyz[bond_rows[bi]]
                # bring the bond to the ring's image
                bxyz = bxyz + _mic(bxyz.mean(axis=0) - coms[k], Lx, Ly) - (bxyz.mean(axis=0) - coms[k])
                if _pierced(rxyz, bxyz):
                    conflicts.append(Conflict('pierce', int(bond_owner[bi]), int(rc)))
    return conflicts


def _mic(d, Lx, Ly):
    d = np.array(d, float, copy=True)
    d[..., 0] -= Lx * np.round(d[..., 0] / Lx)
    d[..., 1] -= Ly * np.round(d[..., 1] / Ly)
    return d


def resolve_conflicts(copies, conflicts, Lx, Ly, seams):
    """Choose copies to remove so that no conflict survives.

    Of the two copies in a conflict, the one whose centroid is closer to a seam goes -- it is the
    lipid that only just made it into the box.  If they are equally close (only possible with no
    seam to tell them apart), the piercer goes for a piercing and the later copy for a clash.
    Piercings are taken first, so a lipid removed for threading a ring is reported as such even
    when it also clashed.  Returns the set of removed copy indices.
    """
    removed = set()
    for c in sorted(conflicts, key=lambda c: c.kind != 'pierce'):
        if c.a in removed or c.b in removed:
            continue
        da = _edge_distance(copies[c.a], Lx, Ly, seams)
        db = _edge_distance(copies[c.b], Lx, Ly, seams)
        if da < db:
            removed.add(c.a)
        elif db < da:
            removed.add(c.b)
        elif c.kind == 'pierce':
            removed.add(c.a)
        else:
            removed.add(max(c.a, c.b))
    return removed


@dataclass
class Quilt:
    copies: list
    Lx: float
    Ly: float
    Lz: float
    zcenter: float
    seams: tuple                # (x, y): whether each direction has a crop seam
    removed: list               # (leaflet, resname, reason, detail) of each removed copy
    stress_free: dict           # leaflet -> lipids the box holds at that leaflet's patch density
    fill: float                 # fraction of its stress-free count each leaflet carries
    n_upper: int
    n_lower: int

    @property
    def box(self):
        return np.diag([self.Lx, self.Ly, self.Lz])

    @property
    def origin(self):
        return np.array([0.5 * self.Lx, 0.5 * self.Ly, self.zcenter])


def _assemble(upper, lower, Lx, Ly, gap, clash_distance):
    seams = seams_of([upper, lower], Lx, Ly)
    copies = tile_leaflet(upper, 'upper', Lx, Ly, gap, seams) + \
        tile_leaflet(lower, 'lower', Lx, Ly, gap, seams)
    conflicts = find_conflicts(copies, Lx, Ly, clash_distance=clash_distance)
    gone = resolve_conflicts(copies, conflicts, Lx, Ly, seams)
    why = {}
    for c in sorted(conflicts, key=lambda c: c.kind != 'pierce'):
        for k in (c.a, c.b):
            if k in gone and k not in why:
                why[k] = (c.kind, c.detail)
    # Equalize.  Crop noise and conflict removals leave the two leaflets at different fractions
    # of their stress-free counts -- 1-5% apart on an ex17 patch tiled against ITSELF -- and that
    # difference is a built-in differential stress, the error this whole build path exists to
    # avoid.  So the leaflet that came out fuller is trimmed, nearest the box edges first, to the
    # other's fraction.  The ratio is then the patches' stress-free ratio to within one lipid.
    sf = {'upper': upper.density() * Lx * Ly, 'lower': lower.density() * Lx * Ly}
    alive = {lf: [k for k, c in enumerate(copies) if c.leaflet == lf and k not in gone]
             for lf in ('upper', 'lower')}
    fill = min(1.0, min(len(alive[lf]) / sf[lf] for lf in alive))
    for lf, ks in alive.items():
        excess = len(ks) - int(round(fill * sf[lf]))
        if excess > 0:
            ks = sorted(ks, key=lambda k: _edge_distance(copies[k], Lx, Ly, seams))
            for k in ks[:excess]:
                gone.add(k)
                why[k] = ('ratio', '')
    removed = [(copies[k].leaflet, copies[k].patch.lipids[copies[k].lipid].resname,
                why[k][0], why[k][1]) for k in sorted(gone)]
    kept = [c for k, c in enumerate(copies) if k not in gone]
    return kept, seams, removed, sf, fill


def build_quilt(upper: Patch, lower: Patch, Lx: float, Ly: float, gap: float,
                clash_distance: float = CLASH_DISTANCE, expand: bool = False):
    """Tile ``upper``'s upper leaflet and ``lower``'s lower leaflet over ``Lx x Ly``, clear every
    seam conflict by removal, equalize the leaflets' fill, and return the :class:`Quilt`.  Both
    patches must be registered.

    Every lipid removed, and the gap, leave lateral room that the relaxation closes, so the quilt
    relaxes to about ``fill`` of the box area.  With ``expand`` the box is enlarged once by
    ``1/sqrt(fill)`` and the quilt rebuilt, so that it relaxes to about ``Lx x Ly`` -- for an
    embedding, where the box was sized to clear a protein.

    The box height is the mean of the two patches' heights, centered so that the upper chamber
    is half of ``upper``'s height and the lower half of ``lower``'s: the water each leaflet's
    patch carried per lipid, which is what the solvation step will refill.
    """
    kept, seams, removed, sf, fill = _assemble(upper, lower, Lx, Ly, gap, clash_distance)
    if expand and fill < 1.0:
        s = 1.0 / np.sqrt(fill)
        Lx, Ly = s * Lx, s * Ly
        kept, seams, removed, sf, fill = _assemble(upper, lower, Lx, Ly, gap, clash_distance)
    hu, hl = upper.cell[2], lower.cell[2]
    return Quilt(copies=kept, Lx=Lx, Ly=Ly, Lz=0.5 * (hu + hl), zcenter=0.25 * (hu - hl),
                 seams=seams, removed=removed, stress_free=sf, fill=fill,
                 n_upper=sum(1 for c in kept if c.leaflet == 'upper'),
                 n_lower=sum(1 for c in kept if c.leaflet == 'lower'))


def write_quilt_pdb(quilt: Quilt, path, aliases=None):
    """Write the quilt's lipids as a topology-free coordinate PDB, in the grid packer's layout:
    one unique resid per lipid, chain ``A``, residue name in the four-column field exactly as
    the conformer PDBs the patches were gridded from carry it, coordinates in columns 31-54.

    ``aliases`` maps a residue name longer than four characters to the placeholder psfgen will
    map back (:meth:`Bilayer.pdb_resname_aliases`); without one, a long name is cut to four
    characters, which psfgen then reads as a different residue or none.

    Frame: the box spans ``[0, Lx) x [0, Ly)`` with the midplane at ``z = 0``, so the xsc written
    alongside uses :attr:`Quilt.origin`.
    """
    if len(quilt.copies) > 9999:
        # bilayer_patch.tcl groups lipids by the 4-column resid, so a wrap would merge two lipids
        raise ValueError(f'{len(quilt.copies)} lipids exceed the 9999 a single quilt PDB can number')
    serial = 0
    with open(path, 'w') as f:
        for resid, c in enumerate(quilt.copies, start=1):
            lip = c.patch.lipids[c.lipid]
            rn = (aliases or {}).get(str(lip.resname), str(lip.resname)[:4])
            for nm, (x, y, z) in zip(c.patch.names[lip.rows], c.xyz):
                serial += 1
                nm = str(nm)
                name_field = nm if len(nm) == 4 else f' {nm:<3}'
                f.write(f'ATOM  {serial % 100000:5d} {name_field:<4} {rn:<4}A{resid:4d}    '
                        f'{x:8.3f}{y:8.3f}{z:8.3f}  0.00  0.00\n')
        f.write('END\n')
    return os.fspath(path)
