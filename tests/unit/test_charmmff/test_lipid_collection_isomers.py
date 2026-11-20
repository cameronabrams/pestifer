# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Double-bond isomers in the lipid conformer collection.

Found 2026-11-20 comparing example 17 against a CHARMM-GUI build of the same composition: about
two-thirds of pestifer's oleoyl C9=C10 bonds were TRANS, already in the first minimized structure
and never restored by MD (``X CEL1 CEL1 X`` has multiplicity 2, so trans is a minimum too).  The
athermal-MC sampler pivoted every all-carbon tail bridge, the double bond included.  The shipped
``lipid.tgz`` had a flipped double bond in 180 of its 196 unsaturated entries.

Each entry's ``<RESI>-init.pdb`` is built by psfgen from CHARMM's internal coordinates, so it
carries the isomer CHARMM defines (cis for oleoyl, trans for the sphingosine C4=C5).  That is the
reference every conformer is checked against.

The double bonds are found here from BOND LENGTH in the init structure, not from the
coordination rule the sampler uses (``athermal_mc._is_multiple_bond``), so the check does not share
the code it is checking.
"""
import glob
import os
import shutil
import tarfile
import tempfile

import numpy as np
import pytest

import pestifer
from pestifer.charmmff.athermal_mc import _element_of_mass, build_lipid_mc, run_mc

_DOUBLE_BOND_MAX = 1.42   # A; C=C is ~1.34, a C-C single bond ~1.53


def _collection_tarball():
    root = os.path.dirname(pestifer.__file__)
    hits = sorted(glob.glob(os.path.join(root, 'resources', 'charmmff', '*', 'pdbrepository', 'lipid.tgz')))
    assert hits, 'no shipped lipid.tgz found'
    return hits[-1]


@pytest.fixture(scope='module')
def collection():
    """The shipped collection, extracted once into a private temporary directory."""
    tmp = tempfile.mkdtemp(prefix='lipid-isomers-')
    with tarfile.open(_collection_tarball()) as tf:
        tf.extractall(tmp, filter='data')
    yield os.path.join(tmp, 'lipid')
    shutil.rmtree(tmp, ignore_errors=True)


def read_pdb(path):
    names, xyz = [], []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(('ATOM', 'HETATM')):
                names.append(ln[12:16].strip())
                xyz.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return names, np.array(xyz)


def read_psf(path):
    """(names, masses, bonds as 0-based index pairs) from an X-PLOR PSF."""
    with open(path) as fh:
        lines = fh.readlines()
    i = next(k for k, ln in enumerate(lines) if '!NATOM' in ln)
    n = int(lines[i].split()[0])
    rows = [lines[i + 1 + k].split() for k in range(n)]
    names = [r[4] for r in rows]
    masses = [float(r[7]) for r in rows]
    j = next(k for k, ln in enumerate(lines) if '!NBOND' in ln)
    nb = int(lines[j].split()[0])
    vals, k = [], j + 1
    while len(vals) < 2 * nb:
        vals += [int(x) - 1 for x in lines[k].split()]
        k += 1
    return names, masses, [tuple(p) for p in np.array(vals).reshape(-1, 2)]


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def double_bond_quartets(names, masses, bonds, xyz, *also):
    """Name quartets ``(i, a, b, j)`` spanning each C=C: a C-C bond short in ``xyz`` and in every
    structure in ``also``.

    One structure is not enough.  psfgen's IC build sometimes draws a single bond short -- it
    did for ``C1S-C2S`` in CER3E and ``C12-C13`` in two cardiolipin arms whose topology has them
    saturated -- and the minimized conformers then show it at 1.53 A.  A C=C is short in both.
    """
    heavy = {k for k, m in enumerate(masses) if _element_of_mass(m) != 'H'}
    nbr = {k: [] for k in heavy}
    for a, b in bonds:
        if a in heavy and b in heavy:
            nbr[a].append(b)
            nbr[b].append(a)
    out = []
    for a, b in bonds:
        if a in heavy and b in heavy and _element_of_mass(masses[a]) == 'C' \
                and _element_of_mass(masses[b]) == 'C' \
                and all(np.linalg.norm(x[a] - x[b]) < _DOUBLE_BOND_MAX for x in (xyz, *also)):
            ia = sorted(k for k in nbr[a] if k != b)
            jb = sorted(k for k in nbr[b] if k != a)
            if ia and jb:
                out.append(tuple(names[k] for k in (ia[0], a, b, jb[0])))
    return out


def is_trans(quartet, names, xyz):
    idx = {n: k for k, n in enumerate(names)}
    return abs(dihedral(*[xyz[idx[n]] for n in quartet])) > 90.0


def entry_files(entry_dir):
    resi = os.path.basename(entry_dir).split('__')[0]
    init_pdb = os.path.join(entry_dir, f'{resi}-init.pdb')
    init_psf = os.path.join(entry_dir, f'{resi}-init.psf')
    confs = sorted(p for p in glob.glob(os.path.join(entry_dir, f'{resi}-[0-9]*.pdb')))
    return init_pdb, init_psf, confs


def flipped_bonds(entry_dir):
    """``[(conformer, quartet)]`` whose double bond is on the other side of 90 deg from the init."""
    init_pdb, init_psf, confs = entry_files(entry_dir)
    names, masses, bonds = read_psf(init_psf)
    pnames, xyz = read_pdb(init_pdb)
    assert pnames == names, f'{entry_dir}: init pdb and psf disagree on atom order'
    first = read_pdb(confs[0])[1] if confs else None
    quartets = double_bond_quartets(names, masses, bonds, xyz, *([first] if confs else []))
    ref = {q: is_trans(q, names, xyz) for q in quartets}
    bad = []
    for c in confs:
        cn, cx = read_pdb(c)
        for q in quartets:
            if is_trans(q, cn, cx) != ref[q]:
                bad.append((os.path.basename(c), q))
    return quartets, bad


# --- the detector, proven able to see what it is for ----------------------------------------

class TestTheDetector:

    def test_oleoyl_is_cis_in_the_ic_built_reference(self, collection):
        init_pdb, init_psf, _ = entry_files(os.path.join(collection, 'POPC'))
        names, masses, bonds = read_psf(init_psf)
        _, xyz = read_pdb(init_pdb)
        q = double_bond_quartets(names, masses, bonds, xyz)
        assert [x[1:3] for x in q] == [('C29', 'C210')]
        assert not is_trans(q[0], names, xyz)

    def test_sphingosine_is_trans_in_the_ic_built_reference(self, collection):
        init_pdb, init_psf, _ = entry_files(os.path.join(collection, 'PSM'))
        names, masses, bonds = read_psf(init_psf)
        _, xyz = read_pdb(init_pdb)
        q = double_bond_quartets(names, masses, bonds, xyz)
        assert ('C4S', 'C5S') in [x[1:3] for x in q]
        assert is_trans(next(x for x in q if x[1:3] == ('C4S', 'C5S')), names, xyz)

    def test_dioleoyl_has_two_double_bonds(self, collection):
        init_pdb, init_psf, _ = entry_files(os.path.join(collection, 'DOPC'))
        names, masses, bonds = read_psf(init_psf)
        _, xyz = read_pdb(init_pdb)
        assert len(double_bond_quartets(names, masses, bonds, xyz)) == 2

    def test_a_bond_turned_180_degrees_reads_as_flipped(self, collection):
        """Rotate POPC's C10-side subtree about C29=C210 by 180 degrees: cis must read trans."""
        import networkx as nx
        init_pdb, init_psf, _ = entry_files(os.path.join(collection, 'POPC'))
        names, masses, bonds = read_psf(init_psf)
        _, xyz = read_pdb(init_pdb)
        q = double_bond_quartets(names, masses, bonds, xyz)[0]
        a, b = names.index('C29'), names.index('C210')
        g = nx.Graph(bonds)
        g.remove_edge(a, b)
        tip = sorted(nx.node_connected_component(g, b))
        axis = (xyz[b] - xyz[a]) / np.linalg.norm(xyz[b] - xyz[a])
        R = 2 * np.outer(axis, axis) - np.eye(3)          # 180-degree rotation about `axis`
        turned = xyz.copy()
        turned[tip] = (xyz[tip] - xyz[b]) @ R.T + xyz[b]
        assert not is_trans(q, names, xyz)
        assert is_trans(q, names, turned)


# --- the sampler must not pivot a double bond -----------------------------------------------

def _popc_mol(collection, **kw):
    init_pdb, init_psf, _ = entry_files(os.path.join(collection, 'POPC'))
    names, masses, bonds = read_psf(init_psf)
    _, xyz = read_pdb(init_pdb)
    heads = [names.index(n) for n in ('N', 'O21', 'O31')]
    tails = [names.index(n) for n in ('C218', 'C316')]
    mol = build_lipid_mc(xyz, None, masses, bonds, heads, tails, np.full(len(names), 0.5),
                         cylinder_radius=float('inf'), **kw)
    return mol, names, xyz


class TestTheSamplerKeepsDoubleBonds:

    def test_the_double_bond_is_not_a_rotatable_bond(self, collection):
        mol, names, _ = _popc_mol(collection)
        pivots = {frozenset((names[r.a], names[r.b])) for r in mol.rotatable}
        assert frozenset(('C29', 'C210')) not in pivots
        # its neighbors, the allylic single bonds, still are -- the chain is not frozen
        assert frozenset(('C28', 'C29')) in pivots
        assert frozenset(('C210', 'C211')) in pivots

    def test_sampled_conformers_keep_the_cis_bond(self, collection):
        mol, names, xyz = _popc_mol(collection)
        q = ('C28', 'C29', 'C210', 'C211')
        assert not is_trans(q, names, xyz)
        samples = run_mc(mol, nsamples=20, n_equil=2000, n_decorr=200, max_angle=np.pi, seed=11)
        assert not any(is_trans(q, names, s) for s in samples)
        idx = [names.index(n) for n in q]
        ref = dihedral(*xyz[idx])
        assert all(abs(dihedral(*s[idx]) - ref) < 1e-6 for s in samples)


# --- the shipped collection -----------------------------------------------------------------

@pytest.mark.xfail(strict=True, reason=(
    'the shipped lipid.tgz (2026-10-08, 037cdfe5) was generated by the sampler before it '
    'stopped pivoting double bonds: 178 entries carry a flipped one. '
    'Regenerate and reship the collection, then delete this marker -- strict, so a clean '
    'collection fails here until it is.'))
def test_every_conformer_keeps_its_reference_isomers(collection):
    report, n_unsat = [], 0
    for entry in sorted(os.listdir(collection)):
        d = os.path.join(collection, entry)
        init_pdb, init_psf, confs = entry_files(d)
        if not (os.path.isfile(init_pdb) and os.path.isfile(init_psf) and confs):
            continue
        quartets, bad = flipped_bonds(d)
        if quartets:
            n_unsat += 1
        if bad:
            report.append(f'{entry}: {len(bad)} flipped, e.g. {bad[0][0]} {bad[0][1][1]}={bad[0][1][2]}')
    assert n_unsat > 50, f'only {n_unsat} unsaturated entries found; is the detector running?'
    assert not report, f'{len(report)} entries with flipped double bonds:\n' + '\n'.join(report)
