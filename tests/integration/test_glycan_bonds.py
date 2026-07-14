# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Regression tests for intra-glycan bond connectivity across chain relabeling.

Background
----------
When pestifer relabels a glycan into its own psfgen segment it stores the
segment's residues in BFS (root-outward) order and renumbers them in that
order, so the LINK-derived ``patch`` commands reference BFS resids.  A bug in
the generic-segment coordinate writer emitted ``$sel set resid [list ...]`` with
the resid list in BFS-residue order, but VMD's ``atomselect`` returns atoms in
ascending-serial order and applies the list in *that* order.  For a branched
glycan the two orderings diverge past the first branch point, so the coordinates
ended up numbered differently from the topology and the branch-point bonds were
wired to the wrong residues -- producing intra-glycan "bonds" spanning >10 A.

These tests build a handful of small glycoproteins from the RCSB (miniproteins
with model N-glycans, all glycosylated at ASN 297) and assert that every
inter-residue bond within a glycan segment is a physically real covalent bond.
Each structure carries a *branched* glycan, so the BFS-vs-serial reordering path
is genuinely exercised; 4b7i/4byh additionally carry chloride ions, so a clean
build also guards the ``CL -> CLA`` ion-placement fix.
"""

import math
import re
import textwrap
from pathlib import Path

import pytest

from pestifer.core.config import Config
from pestifer.core.controller import Controller

# Small glycoproteins with branched N-glycans (and, for 4b7i/4byh, chloride ions).
GLYCAN_PDB_IDS = ['2wah', '4b7i', '4byh']

# Real sugar-sugar glycosidic bonds are ~1.4-1.5 A; the relabeling bug produced
# branch-point "bonds" of 10-16 A.  2.0 A cleanly separates correct from broken.
MAX_INTRA_GLYCAN_BOND = 2.0

_GLYCAN_SEG = re.compile(r'.+G\d\d$')  # pestifer glycan segnames, e.g. 'AG01', 'DG10'


def _build(pdbid: str, workdir: Path) -> dict:
    """Run a minimal fetch + psfgen build for ``pdbid`` in ``workdir``; return the task report."""
    cfg = workdir / f'{pdbid}.yaml'
    cfg.write_text(textwrap.dedent(f"""
        title: {pdbid} intra-glycan bond regression
        tasks:
          - fetch:
              sourceID: {pdbid}
              source_format: pdb
          - psfgen:
              source:
                biological_assembly: 1
    """).lstrip())
    config = Config(userfile=str(cfg)).configure_new()
    controller = Controller().configure(config)
    return controller.do_tasks()


def _parse_psf_pdb(psf: Path, pdb: Path):
    """Return (atoms, bonds, coords).

    atoms: {atomid -> (segname, resid, resname, name)}
    bonds: [(atomid, atomid), ...]
    coords: {serial -> (x, y, z)}  (pestifer writes PSF and PDB in the same atom order)
    """
    lines = psf.read_text().splitlines()
    ai = next(k for k, l in enumerate(lines) if '!NATOM' in l)
    natom = int(lines[ai].split()[0])
    atoms = {}
    for l in lines[ai + 1:ai + 1 + natom]:
        p = l.split()
        atoms[int(p[0])] = (p[1], p[2], p[3], p[4])
    bi = next(k for k, l in enumerate(lines) if '!NBOND' in l)
    nbond = int(lines[bi].split()[0])
    ints: list[int] = []
    k = bi + 1
    while len(ints) < nbond * 2:
        ints += [int(x) for x in lines[k].split()]
        k += 1
    bonds = [(ints[2 * b], ints[2 * b + 1]) for b in range(nbond)]

    coords = {}
    idx = 0
    for l in pdb.read_text().splitlines():
        if l.startswith(('ATOM', 'HETATM')):
            idx += 1
            coords[idx] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return atoms, bonds, coords


def _inter_residue_glycan_bonds(atoms, bonds):
    """Yield (a, b, seg_a) for each bond that joins two different residues of one glycan segment."""
    for a, b in bonds:
        sa, sb = atoms[a], atoms[b]
        if _GLYCAN_SEG.match(sa[0]) and sa[0] == sb[0] and sa[1] != sb[1]:
            yield a, b, sa[0]


@pytest.mark.slow
@pytest.mark.parametrize('pdbid', GLYCAN_PDB_IDS)
def test_intraglycan_bond_geometry(pdbid, tmp_path, monkeypatch):
    """Every inter-residue bond inside a glycan segment must be a real covalent bond."""
    monkeypatch.chdir(tmp_path)
    report = _build(pdbid, tmp_path)

    # A clean build (all tasks result == 0) is itself the guard for the CL -> CLA
    # ion-placement fix: an unplaced chloride trips psfgen's origin-atom check and
    # fails the psfgen task.
    failed = {r['taskname']: r['result'] for r in report.values() if r['result'] != 0}
    assert not failed, f'{pdbid}: build task(s) failed: {failed}'

    psf = tmp_path / 'my_system.psf'
    pdb = tmp_path / 'my_system.pdb'
    assert psf.exists() and pdb.exists(), f'{pdbid}: expected build outputs not found'
    atoms, bonds, coords = _parse_psf_pdb(psf, pdb)

    long_bonds = []
    n_inter = 0
    children: dict[tuple, set] = {}
    for a, b, seg in _inter_residue_glycan_bonds(atoms, bonds):
        n_inter += 1
        d = math.dist(coords[a], coords[b])
        if d > MAX_INTRA_GLYCAN_BOND:
            sa, sb = atoms[a], atoms[b]
            long_bonds.append(f'{seg}:{sa[3]}{sa[1]}--{sb[3]}{sb[1]} ({d:.1f} A)')
        # track parent->children to confirm a branch is present (parent donates the
        # non-anomeric atom; the child contributes its anomeric C1/C2)
        sa, sb = atoms[a], atoms[b]
        if sa[3] in ('C1', 'C2'):
            parent, child = (sb[0], sb[1]), (sa[0], sa[1])
        elif sb[3] in ('C1', 'C2'):
            parent, child = (sa[0], sa[1]), (sb[0], sb[1])
        else:
            continue
        children.setdefault(parent, set()).add(child)

    assert n_inter > 0, f'{pdbid}: no inter-residue glycan bonds found -- structure/parse problem'
    # ensure this structure genuinely exercises the branched-glycan reordering path
    assert any(len(kids) >= 2 for kids in children.values()), (
        f'{pdbid}: no branched glycan detected; this case would not exercise the '
        f'BFS-vs-serial reordering bug')
    assert not long_bonds, (
        f'{pdbid}: {len(long_bonds)} intra-glycan bond(s) span implausible distances '
        f'(> {MAX_INTRA_GLYCAN_BOND} A), indicating branch connectivity wired to the '
        f'wrong residues: {long_bonds}')
