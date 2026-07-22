# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Regression test for biological-assembly (BIOMT) image generation.

Assembly image coordinates are now authored in pure numpy
(``PsfgenScripter._author_subsegment_pdb`` applies the image ``tmat`` via
``Transform.apply``) instead of VMD's ``$sel move {tmat}; writepdb``.  A correct
image is a *rigid-body-displaced copy* of the asymmetric-unit chain: superposing
the image onto the AU chain must give ~0 RMSD (it is a rigid copy), while the
image must sit far from the AU chain in place (the transform actually moved it).

This builds 4zmj (HIV-1 Env trimer) as biological assembly 1 -- a 3-fold
assembly, so the gp120 protein chain and its N-glycans each appear as one AU
copy plus two symmetry images -- and checks both invariants on the protein.
"""

import numpy as np
import pytest

from pestifer.core.config import Config
from pestifer.core.controller import Controller
from pestifer.util.coord import kabsch

pytestmark = pytest.mark.needs_tools


def _ca_by_chain(pdb):
    """{chainID -> {resid -> (x,y,z)}} for CA atoms of a built pdb."""
    out = {}
    for l in open(pdb):
        if l.startswith('ATOM') and l[12:16].strip() == 'CA':
            chain = l[21]
            resid = int(l[22:26])
            out.setdefault(chain, {})[resid] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return out


@pytest.mark.slow
def test_biomt_images_are_displaced_rigid_copies(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = tmp_path / 'biomt.yaml'
    cfg.write_text(
        "title: BIOMT image regression\n"
        "tasks:\n"
        "  - fetch:\n"
        "      sourceID: 4zmj\n"
        "  - psfgen:\n"
        "      source:\n"
        "        biological_assembly: 1\n"
    )
    config = Config(userfile=str(cfg)).configure_new()
    report = Controller().configure(config).do_tasks()
    failed = {r['taskname']: r['result'] for r in report.values() if r['result'] != 0}
    assert not failed, f'build task(s) failed: {failed}'

    ca = _ca_by_chain(tmp_path / 'my_system.pdb')
    # 4zmj assembly 1 has three protein copies of each subunit; each chain here is a full copy
    big = sorted((c for c in ca if len(ca[c]) > 50), key=lambda c: -len(ca[c]))
    assert len(big) >= 2, f'expected >=2 large protein copies, got chains {big}'

    au_chain = big[0]
    n_images_checked = 0
    for img_chain in big[1:]:
        shared = sorted(set(ca[au_chain]) & set(ca[img_chain]))
        if len(shared) < 50:
            continue                      # a different subunit, not an image of au_chain
        P = np.array([ca[au_chain][r] for r in shared])   # AU copy
        Q = np.array([ca[img_chain][r] for r in shared])  # image copy
        # (a) the transform actually displaced the image away from the AU chain
        direct = float(np.sqrt(np.mean(np.sum((P - Q) ** 2, axis=1))))
        # (b) the image is nonetheless a rigid copy: optimal superposition collapses it
        _, _, fit_rmsd = kabsch(Q, P)
        assert fit_rmsd < 0.05, (
            f'chain {img_chain} is not a rigid copy of {au_chain}: fit RMSD {fit_rmsd:.3f} A')
        if direct > 5.0:                  # a genuine (non-identity) symmetry image
            n_images_checked += 1
    assert n_images_checked >= 1, 'no displaced (non-identity) symmetry image found to verify'
