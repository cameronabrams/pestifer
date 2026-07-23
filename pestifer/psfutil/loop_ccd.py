# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Physics-based loop closure by cyclic coordinate descent (CCD).

Offline, deterministic, and free of VMD/NAMD: everything here operates on numpy
coordinate arrays. This is the geometry core of P1 of the loop-modeling redesign
(``docs/design/loop-modeling.md``): it closes a built (``guesscoord``) loop onto its
downstream anchor by distributing the closure across the loop's *own* backbone
dihedrals, rather than dragging one end across space (the steered-MD artifact this
replaces).

The core is deliberately protein-agnostic -- it moves index-addressed points about
index-addressed bond axes toward a fixed target -- so it can be unit-tested on synthetic
geometry. The task layer is responsible for mapping a protein loop's backbone atoms onto
these arrays.

Reference: A. A. Canutescu & R. L. Dunbrack Jr., "Cyclic coordinate descent: A robotics
algorithm for protein loop closure", Protein Science 12:963-972 (2003).
"""
import numpy as np

from ..util.coord import rotate_points_about_axis


def dihedral_deg(p1, p2, p3, p4):
    """
    Signed dihedral (degrees, (-180, 180]) of the four points p1-p2-p3-p4 (each (3,)).

    Sign follows the IUPAC / VMD ``measure dihed`` convention (verified against VMD), so
    values are directly comparable to standard Ramachandran phi/psi (e.g. alpha ~ -63/-43).
    """
    b1 = np.asarray(p2, float) - np.asarray(p1, float)
    b2 = np.asarray(p3, float) - np.asarray(p2, float)
    b3 = np.asarray(p4, float) - np.asarray(p3, float)
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return float(np.degrees(np.arctan2(-y, x)))   # -y matches VMD's sign


def set_dihedral(coords, i, j, k, l, target_deg, moving_mask):
    """
    Rotate the atoms flagged by ``moving_mask`` about the j-k bond so that the dihedral
    ``coords[i]-coords[j]-coords[k]-coords[l]`` equals ``target_deg``.

    ``coords`` is (M, 3); ``i, j, k, l`` are indices into it; ``moving_mask`` is a length-M
    boolean array of the atoms downstream of the j-k bond (must include ``l``, must exclude
    ``i, j, k``). Returns a rotated copy; the input is not modified.
    """
    X = np.array(coords, dtype=float)
    current = dihedral_deg(X[i], X[j], X[k], X[l])
    delta = target_deg - current
    mask = np.asarray(moving_mask, dtype=bool)
    # With the VMD-convention dihedral above, a right-handed rotation of the downstream atoms
    # about the j->k axis increases the dihedral, so rotate by +delta about (X[k]-X[j]).
    X[mask] = rotate_points_about_axis(X[mask], X[j], X[k] - X[j], delta)
    return X

# Analytic Ramachandran basins for the P1 initial guess: name -> (phi_deg, psi_deg, weight).
# Zero-data (no derived library yet -- that is P2). Centers are the conventional means of the
# major coil-accessible basins; weights are a coil-like mixture (helix/sheet/PPII dominant,
# left-handed rare).
RAMACHANDRAN_BASINS = {
    'alphaR': (-63.0, -43.0, 0.42),   # right-handed alpha
    'beta':   (-135.0, 135.0, 0.33),  # beta / extended
    'ppII':   (-75.0, 145.0, 0.20),   # polyproline II
    'alphaL': (60.0, 45.0, 0.05),     # left-handed alpha (mostly Gly)
}


def sample_backbone_dihedrals(rng, n, jitter_deg=10.0):
    """
    Draw ``n`` (phi, psi) pairs, in degrees, from the :data:`RAMACHANDRAN_BASINS` mixture.

    Deterministic for a fixed ``rng`` (a ``numpy.random.Generator``): a build seeded the
    same way samples the same starting backbone. Each residue picks a basin by weight, then
    a Gaussian jitter (std ``jitter_deg``) is added to each angle. Returns an ``(n, 2)``
    array of ``[phi, psi]`` in degrees (not wrapped; callers that need [-180, 180) should
    wrap).
    """
    names = list(RAMACHANDRAN_BASINS)
    weights = np.array([RAMACHANDRAN_BASINS[k][2] for k in names], dtype=float)
    weights /= weights.sum()
    picks = rng.choice(len(names), size=n, p=weights)
    centers = np.array([[RAMACHANDRAN_BASINS[names[p]][0],
                         RAMACHANDRAN_BASINS[names[p]][1]] for p in picks], dtype=float)
    return centers + rng.normal(0.0, jitter_deg, size=(n, 2))


_BACKBONE = ('N', 'CA', 'C', 'O')


def _resid_key(line):
    """Residue identifier from a PDB ATOM/HETATM line, matching how
    :meth:`Molecule.protein_loop_gaps` reports resids: the plain ``int`` sequence number
    (cols 23-26) when there is no insertion code, else the ``'<seq><icode>'`` string (e.g.
    ``'185A'``) so insertion-coded loops (HIV Env V1/V2, etc.) key consistently."""
    resseq = int(line[22:26])
    icode = line[26].strip()
    return f'{resseq}{icode}' if icode else resseq


def backbone_from_pdb(pdb_path, chainID=None, segname=None):
    """
    Parse backbone atoms (N, CA, C, O) from a PDB into ``{resid: {atomname: xyz}}``.

    Selects records by ``segname`` (cols 73-76) when given, else by ``chainID`` (col 22),
    else all. resid is the sequence number (int) or ``'<seq><icode>'`` when an insertion code
    is present (see :func:`_resid_key`). Later records win on a duplicate (resid, atomname) --
    callers should pass an altloc-free / single-model PDB.
    """
    out = {}
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            name = line[12:16].strip()
            if name not in _BACKBONE:
                continue
            if segname is not None and line[72:76].strip() != segname:
                continue
            if chainID is not None and line[21:22].strip() != chainID:
                continue
            try:
                resid = _resid_key(line)
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            out.setdefault(resid, {})[name] = xyz
    return out


def build_loop_ccd_problem(backbone, loop_resids, n_anchor_resid):
    """
    Assemble the CCD arrays for closing one loop, from a per-residue backbone map.

    Parameters
    ----------
    backbone : dict[int, dict[str, ndarray]]
        ``{resid: {atomname: xyz}}`` for at least the atoms named in ``N, CA, C, O`` of the
        loop residues plus the N-terminal anchor's C (needed as the phi reference of the first
        loop residue).
    loop_resids : sequence[int]
        The loop residues, N->C order (the movable residues).
    n_anchor_resid : int
        The resolved residue just before the loop (fixed); its C anchors the first phi.

    Returns
    -------
    dict
        A mapping with keys::

            coords       : (M, 3) backbone atoms, order = [n_anchor:C] + per loop residue N,CA,C,O
            row          : {(resid, atomname): row index into coords}
            bonds        : list of (a_row, b_row) rotatable backbone bonds, N->C
                           (each loop residue contributes phi = N-CA and psi = CA-C)
            moving_masks : list of bool arrays (len M), atoms downstream of each bond
            end_idx      : rows of the last loop residue's N, CA, C (the end effector)

    The caller supplies the closure ``target`` (the desired positions of ``end_idx``): the
    native positions for a delete-and-rebuild test, or the downstream anchor's geometry for a
    real build.
    """
    order = []          # (resid, atomname) in coords order
    coords = []
    coords.append(backbone[n_anchor_resid]['C']); order.append((n_anchor_resid, 'C'))
    for r in loop_resids:
        for an in _BACKBONE:
            if an in backbone[r]:
                coords.append(backbone[r][an]); order.append((r, an))
    coords = np.asarray(coords, dtype=float)
    row = {key: i for i, key in enumerate(order)}
    M = len(coords)

    def downstream_mask(after_row):
        m = np.zeros(M, dtype=bool)
        m[after_row + 1:] = True
        return m

    bonds, masks = [], []
    for r in loop_resids:
        n, ca, c = row[(r, 'N')], row[(r, 'CA')], row[(r, 'C')]
        # phi: rotate about N-CA -> moves everything after CA (C, O, downstream residues)
        bonds.append((n, ca)); masks.append(downstream_mask(ca))
        # psi: rotate about CA-C -> moves everything after C (O, downstream residues)
        bonds.append((ca, c)); masks.append(downstream_mask(c))
    last = loop_resids[-1]
    end_idx = [row[(last, 'N')], row[(last, 'CA')], row[(last, 'C')]]
    return {'coords': coords, 'row': row, 'bonds': bonds,
            'moving_masks': masks, 'end_idx': end_idx}


def loop_atoms_from_pdb(pdb_path, resids, segname=None, chainID=None):
    """
    Parse ALL atoms of the given ``resids`` from a PDB into an ordered structure.

    Returns ``(order, coords, serials)`` where ``order[k] = (resid, atomname)``,
    ``coords`` is (M, 3), and ``serials[k]`` is the PDB serial of row k (for writing
    coordinates back with :func:`pestifer.util.coord.pdb_replace_coords`). Atoms are kept in
    file order within each residue and in ``resids`` order across residues.
    """
    want = set(resids)
    per = {r: [] for r in resids}
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            if segname is not None and line[72:76].strip() != segname:
                continue
            if chainID is not None and line[21:22].strip() != chainID:
                continue
            try:
                resid = _resid_key(line)
            except ValueError:
                continue
            if resid not in want:
                continue
            name = line[12:16].strip()
            serial = int(line[6:11])
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            per[resid].append((name, serial, xyz))
    order, coords, serials = [], [], []
    for r in resids:
        for (name, serial, xyz) in per[r]:
            order.append((r, name)); coords.append(xyz); serials.append(serial)
    return order, np.asarray(coords, dtype=float), serials


# residue-p atoms that stay fixed under a phi (N-CA) rotation: the N-side backbone
_PHI_FIXED = {'N', 'HN', 'H', 'HT1', 'HT2', 'HT3', 'CA'}
# residue-p atoms that move under a psi (CA-C) rotation: only the carbonyl/terminal oxygens
_PSI_MOVERS = {'O', 'OT1', 'OT2', 'OXT', 'OT'}


def build_loop_problem(order, coords, loop_resids):
    """
    Build the full-atom CCD arrays for a loop, with topology-correct moving masks.

    Each loop residue contributes two rotatable backbone bonds: phi (N-CA) and psi (CA-C),
    swept N->C. A phi rotation carries the residue's sidechain and carbonyl plus all
    downstream residues; a psi rotation carries only the carbonyl O and all downstream
    residues (the sidechain stays on the CA side). So sidechain atoms move rigidly with the
    backbone -- essential, since CCD only rotates about backbone bonds.

    Parameters
    ----------
    order : list[tuple[int, str]]
        ``(resid, atomname)`` per row, from :func:`loop_atoms_from_pdb`.
    coords : (M, 3) array
    loop_resids : sequence[int]  (N->C order)

    Returns
    -------
    dict with keys ``coords`` (copy), ``row`` {(resid,atomname):i}, ``bonds``, ``moving_masks``.
    Callers pick ``end_idx``/``target`` (e.g. last residue's CA/C/O onto an anchor target).
    """
    M = len(order)
    row = {}
    for i, key in enumerate(order):
        row[key] = i           # later duplicate names (altloc) win; caller should pass single-altloc
    pos = {r: p for p, r in enumerate(loop_resids)}  # loop position of each resid

    def resid_pos(r):
        return pos.get(r, len(loop_resids))          # non-loop (shouldn't occur) sorts last

    bonds, masks = [], []
    for p, r in enumerate(loop_resids):
        n, ca, c = row[(r, 'N')], row[(r, 'CA')], row[(r, 'C')]
        phi_mask = np.zeros(M, dtype=bool)
        psi_mask = np.zeros(M, dtype=bool)
        for i, (rr, name) in enumerate(order):
            downstream = resid_pos(rr) > p
            if downstream:
                phi_mask[i] = True
                psi_mask[i] = True
            elif rr == r:
                if name not in _PHI_FIXED:
                    phi_mask[i] = True
                if name in _PSI_MOVERS:
                    psi_mask[i] = True
        bonds.append((n, ca)); masks.append(phi_mask)
        bonds.append((ca, c)); masks.append(psi_mask)
    return {'coords': np.array(coords, dtype=float), 'row': row,
            'bonds': bonds, 'moving_masks': masks}


def apply_backbone_dihedrals(coords, prob, loop_resids, phipsi, prev_C, next_N):
    """
    Set each loop residue's (phi, psi) to sampled values, in N->C order, on the full-atom
    ``prob`` from :func:`build_loop_problem`.

    Rotates about each residue's N-CA (phi) and CA-C (psi) bonds using that bond's
    topology-correct moving mask, so whole residues (sidechains included) follow. References
    for the terminal torsions come from outside the loop: ``prev_C`` (the N-anchor's C, for
    the first residue's phi) and ``next_N`` (the C-anchor's N, for the last residue's psi).

    Parameters
    ----------
    coords : (M, 3)  -- not modified; a copy is returned
    prob : dict from build_loop_problem (provides ``row`` and ``moving_masks``)
    loop_resids : sequence[int]  (N->C)
    phipsi : (len(loop_resids), 2) array of (phi, psi) in degrees
    prev_C, next_N : (3,) reference positions

    Returns
    -------
    np.ndarray (M, 3)
    """
    X = np.array(coords, dtype=float)
    r = prob['row']
    masks = prob['moving_masks']
    for p, resid in enumerate(loop_resids):
        N, CA, C = r[(resid, 'N')], r[(resid, 'CA')], r[(resid, 'C')]
        pC = X[r[(loop_resids[p - 1], 'C')]] if p > 0 else np.asarray(prev_C, float)
        nN = X[r[(loop_resids[p + 1], 'N')]] if p < len(loop_resids) - 1 else np.asarray(next_N, float)
        # phi = prevC-N-CA-C : rotate downstream of N-CA
        cur = dihedral_deg(pC, X[N], X[CA], X[C])
        X[masks[2 * p]] = rotate_points_about_axis(X[masks[2 * p]], X[N], X[CA] - X[N],
                                                   phipsi[p, 0] - cur)
        # psi = N-CA-C-nextN : rotate downstream of CA-C. The LAST residue's psi is defined by
        # the fixed C-anchor N (not a loop atom), so it cannot be set by rotating loop atoms --
        # it is determined by CCD closure. Skip it.
        if p < len(loop_resids) - 1:
            cur = dihedral_deg(X[N], X[CA], X[C], nN)
            X[masks[2 * p + 1]] = rotate_points_about_axis(X[masks[2 * p + 1]], X[CA], X[C] - X[CA],
                                                           phipsi[p, 1] - cur)
    return X


def place_atom_nerf(a, b, c, bond, angle_deg, dihedral_deg):
    """
    NeRF placement: return the position of atom D given three reference atoms a-b-c and the
    internal coordinates bond (``|c-D|``), angle (b-c-D, degrees), and dihedral (a-b-c-D, degrees).
    """
    a = np.asarray(a, float); b = np.asarray(b, float); c = np.asarray(c, float)
    ang = np.radians(angle_deg); dih = np.radians(dihedral_deg)
    bc = c - b; bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    nn = np.linalg.norm(n)
    n = n / nn if nn > 1e-9 else np.array([0.0, 0.0, 1.0])
    m = np.cross(n, bc)
    # D in the local frame (bc, m, n) about c
    d2 = np.array([-bond * np.cos(ang),
                   bond * np.sin(ang) * np.cos(dih),
                   bond * np.sin(ang) * np.sin(dih)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


# Standard trans-peptide internal coordinates for the loop_last -> anchor junction.
_PEP = dict(C_N=1.33, CA_C=1.52, C_O=1.23,
            CA_N_C=121.7, N_C_CA=116.6, N_C_O=122.7)


def anchor_closure_target(anchor_N, anchor_CA, anchor_C, phi_deg=-120.0):
    """
    Build ideal target positions for the loop's C-terminal residue backbone atoms
    ``(CA, C, O)`` such that that residue forms a good trans peptide bond to the downstream
    anchor (anchor:N).

    Constructed from the fixed anchor N/CA/C by NeRF placement using standard peptide
    geometry; ``phi_deg`` is the assumed C(prev)-N-CA-C torsion (the one free parameter --
    the subsequent psfgen guesscoord+regenerate in ``connect`` refines the exact junction, so
    a generic extended value is fine). Returns ``(target_CA, target_C, target_O)`` to be
    superimposed on the loop's last-residue ``(CA, C, O)`` by CCD.
    """
    # prev C bonded to anchor N: bond C-N, angle CA-N-C, torsion C-CA-N-(prevC)=phi
    C = place_atom_nerf(anchor_C, anchor_CA, anchor_N, _PEP['C_N'], _PEP['CA_N_C'], phi_deg)
    # prev CA from prev C: omega = 180 (trans)
    CA = place_atom_nerf(anchor_CA, anchor_N, C, _PEP['CA_C'], _PEP['N_C_CA'], 180.0)
    # prev O from prev C: in the peptide plane, opposite CA (torsion 0 wrt anchor N)
    O = place_atom_nerf(anchor_CA, anchor_N, C, _PEP['C_O'], _PEP['N_C_O'], 0.0)
    return np.array([CA, C, O])


def optimal_ccd_angle(moving, target, pivot, axis, weights=None):
    """
    Closed-form rotation angle (degrees, right-handed about ``axis``) that minimizes
    ``sum_k w_k |R(theta) moving_k - target_k|^2`` for a rotation about the line through
    ``pivot`` along ``axis``.

    Because each moving atom's component perpendicular to the axis rotates in-plane, the
    objective's theta-dependence is ``A cos(theta) + B sin(theta)`` (up to constants), whose
    maximizer is ``theta* = atan2(B, A)`` with ``A = sum_k w_k (r_k . f_k)`` and
    ``B = sum_k w_k ((k_hat x r_k) . f_k)``, where ``r_k`` / ``f_k`` are the perpendicular
    components (relative to the axis) of ``moving_k`` / ``target_k`` about ``pivot``.

    Parameters
    ----------
    moving, target : array_like, shape (K, 3)
        Current and desired positions of the K end-effector atoms.
    pivot : array_like, shape (3,)
    axis : array_like, shape (3,)  (need not be unit length; a zero axis yields 0.0)
    weights : array_like, shape (K,), optional

    Returns
    -------
    float
        The optimal angle in degrees.
    """
    moving = np.asarray(moving, dtype=float)
    target = np.asarray(target, dtype=float)
    pivot = np.asarray(pivot, dtype=float)
    axis = np.asarray(axis, dtype=float)
    n = np.linalg.norm(axis)
    if n == 0.0 or len(moving) == 0:
        return 0.0
    k = axis / n
    w = np.ones(len(moving)) if weights is None else np.asarray(weights, dtype=float)
    r = moving - pivot
    r = r - np.outer(r @ k, k)          # perpendicular component of moving atoms
    f = target - pivot
    f = f - np.outer(f @ k, k)          # perpendicular component of targets
    A = float(np.sum(w * np.sum(r * f, axis=1)))
    kxr = np.cross(np.broadcast_to(k, r.shape), r)
    B = float(np.sum(w * np.sum(kxr * f, axis=1)))
    return float(np.degrees(np.arctan2(B, A)))


def end_rmsd(coords, end_idx, target):
    """RMSD (Angstrom) between the end-effector atoms ``coords[end_idx]`` and ``target``."""
    d = np.asarray(coords, dtype=float)[end_idx] - np.asarray(target, dtype=float)
    return float(np.sqrt(np.mean(np.sum(d * d, axis=1))))


def ccd_close(coords, bonds, moving_masks, end_idx, target,
              weights=None, max_iters=500, tol=0.08, stall_tol=1e-6):
    """
    Close a loop onto a fixed target by cyclic coordinate descent.

    Sweeps the rotatable ``bonds`` in order (N->C is the natural order for a loop); for each
    bond it applies :func:`optimal_ccd_angle` to bring the end-effector atoms onto ``target``,
    rotating only the atoms flagged by that bond's mask. Repeats until the end RMSD falls
    below ``tol`` or progress stalls.

    Parameters
    ----------
    coords : array_like, shape (M, 3)
        Loop atom coordinates (not modified; a copy is returned).
    bonds : list[tuple[int, int]]
        ``(a, b)`` index pairs into ``coords``; each defines a rotation axis (the line
        through ``coords[a]`` along ``coords[b] - coords[a]``).
    moving_masks : list[array_like[bool]], each length M
        For each bond, which atoms rotate with it. Every end-effector atom must be movable by
        every bond it is downstream of; typically a bond's mask is "all atoms C-terminal to
        it", which includes ``end_idx``.
    end_idx : array_like[int]
        Indices of the end-effector atoms whose landing on ``target`` defines closure.
    target : array_like, shape (len(end_idx), 3)
    weights : array_like, optional
        Per-end-atom weights for the least-squares objective.
    max_iters : int
    tol : float
        Convergence RMSD in Angstrom.
    stall_tol : float
        Stop early if an iteration improves the RMSD by less than this.

    Returns
    -------
    (coords_closed, final_rmsd, iterations) : (np.ndarray, float, int)
    """
    X = np.array(coords, dtype=float)
    target = np.asarray(target, dtype=float)
    end_idx = np.asarray(end_idx, dtype=int)
    last = end_rmsd(X, end_idx, target)
    for it in range(1, max_iters + 1):
        for (a, b), mask in zip(bonds, moving_masks):
            pivot = X[a]
            axis = X[b] - X[a]
            theta = optimal_ccd_angle(X[end_idx], target, pivot, axis, weights)
            mask = np.asarray(mask, dtype=bool)
            X[mask] = rotate_points_about_axis(X[mask], pivot, axis, theta)
        r = end_rmsd(X, end_idx, target)
        if r < tol or (last - r) < stall_tol:
            return X, r, it
        last = r
    return X, end_rmsd(X, end_idx, target), max_iters


def _heavy_mask(order):
    """Boolean mask (len == len(order)) of heavy (non-hydrogen) atoms."""
    return np.array([not name.startswith('H') for (_r, name) in order], dtype=bool)


def heavy_env_coords_from_pdb(pdb_path, exclude_serials=(), segname=None, chainID=None):
    """
    Heavy-atom coordinates of the *environment* -- every atom in ``pdb_path`` except those
    whose PDB serial is in ``exclude_serials`` (typically the loop being scored). Used as the
    frozen backdrop for :func:`loop_clash_report`'s loop-vs-environment check.
    """
    excl = set(int(s) for s in exclude_serials)
    pts = []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            if segname is not None and line[72:76].strip() != segname:
                continue
            if chainID is not None and line[21:22].strip() != chainID:
                continue
            name = line[12:16].strip()
            if name.startswith('H'):
                continue
            if int(line[6:11]) in excl:
                continue
            pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.asarray(pts, dtype=float)


def loop_clash_report(order, coords, loop_resids, env_coords=None,
                      soft_cutoff=2.0, deep_cutoff=1.6, ca_thread=3.0):
    """
    Steric diagnostic for a (closed) loop, over heavy atoms only.

    Distinguishes *soft* overlaps (van-der-Waals bumps that minimization/MD can relax)
    from *topological* failures -- interpenetration or chain threading -- which a local
    optimizer can **never** fix (separating superimposed atoms would require pulling one
    strand back through another). A loop flagged ``topological`` is not a valid MD starting
    structure regardless of how much it is minimized.

    Checks two clash populations: intra-loop (non-adjacent residues, ``|dresid| >= 2``)
    and, if ``env_coords`` is given, loop-vs-environment (the rest of the structure). The
    minimum non-adjacent Ca-Ca distance is a threading proxy (native coil rarely brings
    non-sequential Ca closer than ~4.5 A; well below ``ca_thread`` means strands crossing).

    Parameters
    ----------
    order : list[tuple[int, str]]
        ``(resid, atomname)`` per row of ``coords`` (from :func:`loop_atoms_from_pdb`).
    coords : (M, 3) array
    loop_resids : sequence[int]
    env_coords : (K, 3) array, optional
        Heavy-atom coordinates of every atom NOT in the loop (the frozen environment).
    soft_cutoff, deep_cutoff : float
        Heavy-atom separations (A) defining a soft contact and a deep (interpenetrating)
        contact. ``deep_cutoff`` ~1.6 A is well inside any real heavy-atom vdW contact.
    ca_thread : float
        Non-adjacent Ca-Ca distance (A) below which the backbone is treated as threaded.

    Returns
    -------
    dict with keys ``n_soft``, ``n_deep``, ``n_env_soft``, ``n_env_deep``, ``worst``
    (smallest heavy-atom separation seen, A), ``min_ca`` (A), and ``topological`` (bool:
    any deep contact, or ``min_ca < ca_thread``).
    """
    from scipy.spatial import cKDTree
    X = np.asarray(coords, dtype=float)
    heavy = _heavy_mask(order)
    hX = X[heavy]
    hres = [order[i][0] for i in range(len(order)) if heavy[i]]
    # sequence position of each loop residue (handles insertion codes / mixed int-str resids):
    # "adjacent" means near in chain order, not near in resid number.
    pos = {r: i for i, r in enumerate(loop_resids)}
    report = dict(n_soft=0, n_deep=0, n_env_soft=0, n_env_deep=0, worst=9.99,
                  min_ca=9.99, topological=False)
    if len(hX) == 0:
        return report
    # intra-loop, non-adjacent residues
    tree = cKDTree(hX)
    for a, b in tree.query_pairs(soft_cutoff):
        if abs(pos[hres[a]] - pos[hres[b]]) < 2:
            continue
        d = float(np.linalg.norm(hX[a] - hX[b]))
        report['n_soft'] += 1
        report['worst'] = min(report['worst'], d)
        if d < deep_cutoff:
            report['n_deep'] += 1
    # loop vs environment
    if env_coords is not None and len(env_coords):
        etree = cKDTree(np.asarray(env_coords, dtype=float))
        for k, hits in enumerate(etree.query_ball_point(hX, soft_cutoff)):
            for h in hits:
                d = float(np.linalg.norm(hX[k] - etree.data[h]))
                report['n_env_soft'] += 1
                report['worst'] = min(report['worst'], d)
                if d < deep_cutoff:
                    report['n_env_deep'] += 1
    # min non-adjacent Ca-Ca (threading proxy)
    ca_rows = [i for i, (_r, name) in enumerate(order) if name == 'CA']
    ca_res = [order[i][0] for i in ca_rows]
    ca_xyz = X[ca_rows]
    for i in range(len(ca_rows)):
        for j in range(i + 2, len(ca_rows)):
            if abs(pos[ca_res[i]] - pos[ca_res[j]]) < 2:
                continue
            report['min_ca'] = min(report['min_ca'], float(np.linalg.norm(ca_xyz[i] - ca_xyz[j])))
    report['topological'] = bool(report['n_deep'] or report['n_env_deep']
                                 or report['min_ca'] < ca_thread)
    return report


def _clash_score(report, weight_deep=100.0):
    """Scalar penalty from a :func:`loop_clash_report`: deep (interpenetrating) overlaps
    dominate soft ones. Lower is better; 0 means clash-free."""
    deep = report['n_deep'] + report['n_env_deep']
    soft = report['n_soft'] + report['n_env_soft']
    return weight_deep * deep + soft


def refine_declash_ccd(coords, prob, order, loop_resids, end_idx, target, rng, clash_env_fn,
                       n_iters=250, perturb_deg=25.0, close_tol=0.15, close_iters=400,
                       weight_deep=100.0, T0=8.0, Tmin=0.4, clash_fn=None):
    """
    Iteratively declash a *closed* compact loop by perturb-and-reclose simulated annealing.

    Random sampling rarely lands a crowded loop in the small clash-free region of
    conformation space; iterative descent walks there. Each step perturbs one backbone
    dihedral by a small random angle, re-closes the end onto ``target`` with CCD (a short move
    because the loop stays compact), and scores the result with :func:`loop_clash_report`
    against ``clash_env_fn(coords)`` -- a callback returning the environment coordinates for a
    trial pose, which lets the caller fold in the loop's own symmetry images (so the mutually
    interacting copies of an axial loop declash *together*, radially outward). Moves that lower
    the clash score are accepted; worse moves are accepted with a Metropolis probability under
    a geometrically annealed temperature, to escape local minima. Deterministic for a fixed
    ``rng``.

    Parameters
    ----------
    coords : (M, 3)   -- a closed starting pose (e.g. from :func:`ccd_close`)
    prob : dict from :func:`build_loop_problem` (``bonds``, ``moving_masks``)
    order : list[(resid, atomname)] per row (for :func:`loop_clash_report`)
    loop_resids, end_idx, target : as for :func:`ccd_close`
    rng : numpy.random.Generator
    clash_env_fn : callable (M,3) array -> (K,3) environment coords for scoring that pose
    n_iters, perturb_deg, close_tol, close_iters, weight_deep, T0, Tmin : tuning

    Returns
    -------
    (best_coords, best_report) : (np.ndarray, dict)
    """
    bonds = prob['bonds']
    masks = prob['moving_masks']
    nb = len(bonds)

    def evaluate(X):
        rep = loop_clash_report(order, X, loop_resids, env_coords=clash_env_fn(X))
        return _clash_score(rep, weight_deep), rep

    cur = np.array(coords, dtype=float)
    cur_s, cur_rep = evaluate(cur)
    best, best_s, best_rep = np.array(cur), cur_s, cur_rep
    if best_s == 0:
        return best, best_rep
    early = max(1, nb // 3)              # leading third of the bonds = longest levers
    for it in range(n_iters):
        T = T0 * (Tmin / T0) ** (it / max(1, n_iters - 1))
        # Two move types: an occasional LARGE swing about a leading (long-lever) bond, which
        # rotates the whole loop bulk (e.g. radially outward, away from an assembly axis), and
        # otherwise a small local perturbation that refines shape. The swing supplies the big
        # concerted motion a crowded loop needs to escape; the local move polishes.
        if rng.random() < 0.35:
            b = int(rng.integers(early))
            delta = float(rng.normal(0.0, 90.0))
        else:
            b = int(rng.integers(nb))
            delta = float(rng.normal(0.0, perturb_deg))
        a, c = bonds[b]
        mask = np.asarray(masks[b], dtype=bool)
        trial = np.array(cur)
        trial[mask] = rotate_points_about_axis(trial[mask], trial[a], trial[c] - trial[a], delta)
        if clash_fn is not None:          # guarded re-close: refuse to reintroduce clashes
            trial, rmsd, _ = ccd_close_guarded(trial, bonds, masks, end_idx, target, clash_fn,
                                               tol=close_tol, max_iters=close_iters)
        else:
            trial, rmsd, _ = ccd_close(trial, bonds, masks, end_idx, target,
                                       tol=close_tol, max_iters=close_iters)
        if rmsd > 0.6:                    # could not re-close cleanly -> reject
            continue
        s, rep = evaluate(trial)
        if s <= cur_s or rng.random() < np.exp(-(s - cur_s) / max(T, 1e-6)):
            cur, cur_s, cur_rep = trial, s, rep
            if s < best_s:
                best, best_s, best_rep = np.array(trial), s, rep
                if best_s == 0:
                    break
    return best, best_rep


def ccd_close_guarded(coords, bonds, moving_masks, end_idx, target, clash_fn,
                      tol=0.15, max_iters=500, stall_tol=1e-6):
    """
    CCD closure that keeps the environment present and refuses to clash into it.

    Identical sweep to :func:`ccd_close`, but at each rotatable bond the closure-optimal
    rotation is committed **only if it does not increase** the environment clash count
    (``clash_fn(coords) -> number of deep clashes``); otherwise that bond is left alone this
    pass. So the loop walks its end toward ``target`` without ever being driven into the
    protein background (or, when ``clash_fn`` includes them, its own symmetry images). If a
    clash-free path to the anchor does not exist, the end simply will not fully close -- the
    returned ``rmsd`` reports that honestly rather than forcing an interpenetrating pose.

    Returns ``(coords, final_rmsd, final_clash)``.
    """
    X = np.array(coords, dtype=float)
    target = np.asarray(target, dtype=float)
    end_idx = np.asarray(end_idx, dtype=int)
    cur_clash = clash_fn(X)
    last = end_rmsd(X, end_idx, target)
    for it in range(1, max_iters + 1):
        for (a, b), mask in zip(bonds, moving_masks):
            theta = optimal_ccd_angle(X[end_idx], target, X[a], X[b] - X[a])
            mask = np.asarray(mask, dtype=bool)
            trial = X.copy()
            trial[mask] = rotate_points_about_axis(trial[mask], X[a], X[b] - X[a], theta)
            tc = clash_fn(trial)
            if tc <= cur_clash:
                X = trial
                cur_clash = tc
        r = end_rmsd(X, end_idx, target)
        if r < tol or abs(last - r) < stall_tol:
            break
        last = r
    return X, end_rmsd(X, end_idx, target), cur_clash


def close_one_loop(src_pdb, segname, loop_resids, n_anchor_resid, c_anchor_resid,
                   all_loop_serials, seed, ensemble=10, refine=250, guard=False,
                   max_iters=2000, tol=0.1, extra_env=None, progress_queue=None):
    """
    Close one interior loop against the resolved structure and (optionally) a set of
    already-closed loops -- the standalone, picklable unit the ligate task runs in parallel
    across loops and re-runs sequentially for the rare threaded pair.

    Reads the loop and its environment from ``src_pdb``: the environment is every heavy atom
    NOT in ``all_loop_serials`` (all modeled loops are throwaway) and NOT in this loop's two
    flanking anchors (bonded junctions), plus ``extra_env`` (heavy coords of loops already
    closed, passed during sequential repair). Seeds each loop residue from Ramachandran basins,
    closes onto the anchor by CCD (guarded if ``guard``), polishes by iterative declash
    (``refine`` iterations), and keeps the least-clashing of ``ensemble`` seeds. Deterministic
    for a fixed ``seed``.

    ``progress_queue`` (optional) is a picklable queue -- a ``multiprocessing.Manager().Queue`` when
    called across a process pool -- onto which one item is pushed per ensemble candidate actually
    run, so the ligate driver can render a determinate progress bar with a time-remaining estimate.

    Returns a dict: ``segname``, ``loop`` (resids), ``serials``, ``order``, ``closed`` (M,3),
    ``rep`` (:func:`loop_clash_report`), ``heavy`` (heavy-atom coords), ``ca`` (Ca coords).
    """
    from scipy.spatial import cKDTree
    from scipy.spatial.distance import cdist
    rng = np.random.default_rng(seed)
    bb = backbone_from_pdb(src_pdb, segname=segname)
    order, coords, serials = loop_atoms_from_pdb(src_pdb, loop_resids, segname=segname)
    prob = build_loop_problem(order, coords, loop_resids)
    r = prob['row']
    last = loop_resids[-1]
    end_idx = [r[(last, 'CA')], r[(last, 'C')]]
    target = anchor_closure_target(bb[c_anchor_resid]['N'], bb[c_anchor_resid]['CA'],
                                   bb[c_anchor_resid]['C'])[[0, 1]]
    _ao, _ac, anch_ser = loop_atoms_from_pdb(src_pdb, [n_anchor_resid, c_anchor_resid], segname=segname)
    env = heavy_env_coords_from_pdb(src_pdb, exclude_serials=list(all_loop_serials) + list(anch_ser))
    if extra_env is not None and len(extra_env):
        extra_env = np.asarray(extra_env, dtype=float)
        env = np.vstack([env, extra_env]) if len(env) else extra_env

    heavy_mask = np.array([not nm.startswith('H') for (_rr, nm) in order])
    loop_pos = {rr: i for i, rr in enumerate(loop_resids)}
    hpos = np.array([loop_pos[order[i][0]] for i in range(len(order)) if heavy_mask[i]])
    _intra = np.triu(np.abs(hpos[:, None] - hpos[None, :]) >= 2, k=1)
    env_tree = cKDTree(env) if len(env) else None

    def clash_env_fn(_pose):
        return env

    def clash_fn(X):
        hx = X[heavy_mask]
        deep = 0
        if env_tree is not None:
            deep += int(env_tree.query_ball_point(hx, 1.6, return_length=True).sum())
        deep += int(((cdist(hx, hx) < 1.6) & _intra).sum())
        return deep

    best = None
    for cand in range(max(1, ensemble)):
        phipsi = sample_backbone_dihedrals(rng, len(loop_resids))
        start = apply_backbone_dihedrals(prob['coords'], prob, loop_resids, phipsi,
                                         prev_C=bb[n_anchor_resid]['C'], next_N=bb[c_anchor_resid]['N'])
        if guard:
            closed, _rmsd, _nc = ccd_close_guarded(start, prob['bonds'], prob['moving_masks'],
                                                   end_idx, target, clash_fn,
                                                   tol=max(tol, 0.15), max_iters=max_iters)
        else:
            closed, _rmsd, _it = ccd_close(start, prob['bonds'], prob['moving_masks'],
                                           end_idx, target, tol=tol, max_iters=max_iters)
        if refine > 0:
            closed, rep = refine_declash_ccd(closed, prob, order, loop_resids, end_idx, target,
                                             rng, clash_env_fn, n_iters=refine, close_iters=120,
                                             clash_fn=clash_fn if guard else None)
        else:
            rep = loop_clash_report(order, closed, loop_resids, env_coords=env)
        key = (rep['n_deep'] + rep['n_env_deep'], rep['n_soft'] + rep['n_env_soft'],
               end_rmsd(closed, end_idx, target))
        if best is None or key < best[0]:
            best = (key, closed, rep)
        if progress_queue is not None:
            # one tick per candidate actually run; the driver trues-up on loop completion
            # since an early break (below) means fewer than `ensemble` ticks are emitted.
            try:
                progress_queue.put_nowait(1)
            except Exception:
                pass
        if key[0] == 0 and key[1] == 0:
            break
    _key, closed, rep = best
    ca = closed[[i for i, (_rr, nm) in enumerate(order) if nm == 'CA']]
    return dict(segname=segname, loop=loop_resids, serials=serials, order=order,
                closed=closed, rep=rep, heavy=closed[heavy_mask], ca=ca)
