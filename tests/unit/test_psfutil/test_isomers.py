# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Double-bond isomers held through md stages and checked after them (``md: double_bonds``).

Found 2026-11-24 on example 17: with the sampler and the ring resolver both kept off C=C bonds, the
first minimize of a crowded grid placement still turned oleoyl bonds trans.  The oracle here is
CHARMM's own topology -- the IC a residue is built from -- read through a really provisioned
force-field content, not a table written into this file.  Single molecules come from the shipped
lipid collection's psfgen-built ``-init`` entries, so nothing here needs VMD or NAMD.
"""
import logging
import os
import shutil
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace as NS

import numpy as np
import pytest

from pestifer.core.errors import PestiferBuildError
from pestifer.psfutil.isomers import DoubleBondChecker, isomer_of, isomer_reference
from pestifer.tasks.mdtask import MDTask

from tests.unit.test_charmmff.test_lipid_collection_isomers import _collection_tarball


@pytest.fixture(scope='module')
def collection():
    tmp = tempfile.mkdtemp(prefix='isomers-')
    with tarfile.open(_collection_tarball()) as tf:
        members = [m for m in tf.getmembers()
                   if any(f'lipid/{r}/{r}-init.' in m.name for r in ('POPC', 'PSM', 'CHL1'))]
        tf.extractall(tmp, members=members, filter='data')
    yield Path(tmp) / 'lipid'
    shutil.rmtree(tmp, ignore_errors=True)


@pytest.fixture(scope='module')
def charmm():
    from pestifer.core.resourcemanager import ResourceManager
    cc = ResourceManager().charmmff_content
    cc.provision(resnames=['POPC', 'PSM', 'CHL1'])
    return cc


def files(collection, resi):
    return str(collection / resi / f'{resi}-init.psf'), str(collection / resi / f'{resi}-init.pdb')


def write_pdb(src, dst, xyz):
    i = 0
    with open(src) as f, open(dst, 'w') as g:
        for ln in f:
            if ln.startswith(('ATOM', 'HETATM')):
                x, y, z = xyz[i]
                ln = f'{ln[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{ln[54:]}'
                i += 1
            g.write(ln)


def turned(psf, xyz, a, b, deg=180.0):
    """``xyz`` with everything bonded to ``b``'s side of the a-b bond turned ``deg`` about it."""
    from pestifer.psfutil.isomers import _read_psf
    _, bonds = _read_psf(psf)
    nbr = defaultdict(list)
    for i, j in bonds:
        nbr[i].append(j)
        nbr[j].append(i)
    side, todo = {b}, [b]
    while todo:
        x = todo.pop()
        for y in nbr[x]:
            if y != a and y not in side:
                side.add(y)
                todo.append(y)
    assert a not in side, 'the bond must be a bridge'
    out = xyz.copy()
    axis = (xyz[b] - xyz[a]) / np.linalg.norm(xyz[b] - xyz[a])
    rows = sorted(side)
    p = xyz[rows] - xyz[b]
    t = np.radians(deg)
    out[rows] = (xyz[b] + p * np.cos(t) + np.cross(axis, p) * np.sin(t)
                 + np.outer(p @ axis, axis) * (1 - np.cos(t)))
    return out


def ref(cc):
    return lambda rn: isomer_reference(cc, rn)


class TestTheReferenceIsCharmms:

    def test_oleoyl_is_cis_and_sphingosine_trans(self, charmm):
        table = isomer_reference(charmm, {'POPC', 'PSM'})
        assert isomer_of(table['POPC'][frozenset(('C29', 'C210'))][4]) == 'cis'
        assert isomer_of(table['PSM'][frozenset(('C4S', 'C5S'))][4]) == 'trans'

    def test_unprovisioned_content_still_gives_the_reference_and_stays_unprovisioned(self):
        """A continuation build provisions nothing.  The first end-to-end run of this check read its
        reference only from provisioned residues, found none, and reported 114 double bonds -- two
        of them trans -- as all right."""
        from pestifer.core.resourcemanager import ResourceManager
        cc = ResourceManager().charmmff_content
        assert not cc.provisioned
        table = isomer_reference(cc, {'SOPE', 'PSM'})
        assert isomer_of(table['SOPE'][frozenset(('C29', 'C210'))][4]) == 'cis'
        assert isomer_of(table['PSM'][frozenset(('C4S', 'C5S'))][4]) == 'trans'
        assert not cc.provisioned and len(cc.residues) == 0, 'reading must not provision'

    def test_a_name_nothing_defines_has_no_reference(self, charmm):
        assert isomer_reference(charmm, {'NOTARESIDUE'}) == {}
        assert isomer_reference(None, {'POPC'}) == {}


class TestTheChecker:

    def test_popc_has_one_double_bond_and_it_is_right(self, collection, charmm):
        psf, pdb = files(collection, 'POPC')
        ch = DoubleBondChecker(psf, reference=ref(charmm))
        assert [t[4] for t in ch.terms] == [('C28', 'C29', 'C210', 'C211')]
        assert ch.nreferenced == 1
        assert ch.check_coords(ch.load_coords(pdb)) == []

    def test_psm_trans_is_right(self, collection, charmm):
        psf, pdb = files(collection, 'PSM')
        ch = DoubleBondChecker(psf, reference=ref(charmm))
        assert len(ch.terms) == 1 and ch.terms[0][5] == 'trans'
        assert ch.check_coords(ch.load_coords(pdb)) == []

    def test_a_ring_double_bond_is_not_a_term(self, collection, charmm):
        """Cholesterol's C5=C6 is in ring B: it satisfies the coordination rule and has no isomer."""
        psf, _ = files(collection, 'CHL1')
        assert DoubleBondChecker(psf, reference=ref(charmm)).terms == []

    def test_a_turned_oleoyl_is_a_defect_and_is_not_held(self, collection, charmm):
        psf, pdb = files(collection, 'POPC')
        ch = DoubleBondChecker(psf, reference=ref(charmm))
        rows = ch.terms[0][0]
        xyz = turned(psf, ch.load_coords(pdb), rows[1], rows[2])
        wrong = ch.check_coords(xyz)
        assert len(wrong) == 1
        assert 'C29=C210 is trans' in wrong[0].describe() and 'CHARMM defines it cis' in wrong[0].describe()
        assert ch.restraint_terms(xyz) == [], 'holding a wrong bond would cement it'

    def test_a_right_bond_is_held_at_its_reference(self, collection, charmm):
        psf, pdb = files(collection, 'POPC')
        ch = DoubleBondChecker(psf, reference=ref(charmm))
        held = ch.restraint_terms(ch.load_coords(pdb))
        assert [(h[0], h[1]) for h in held] == [(ch.terms[0][0], 0.0)]

    def test_without_a_reference_it_is_held_as_it_is_and_not_checked(self, collection):
        psf, pdb = files(collection, 'POPC')
        ch = DoubleBondChecker(psf, reference=None)
        xyz = turned(psf, ch.load_coords(pdb), *ch.terms[0][0][1:3])
        assert ch.nreferenced == 0
        assert ch.check_coords(xyz) == []
        assert [h[1] for h in ch.restraint_terms(xyz)] == [180.0]


class _TaskCase:
    """A bare MDTask in a temporary directory with its state pointed at one molecule."""

    def make(self, tmp_path, monkeypatch, collection, charmm, resi='POPC', pdb=None):
        monkeypatch.chdir(tmp_path)
        psf, init = files(collection, resi)
        t = MDTask.__new__(MDTask)
        t.basename, t.taskname = '01-06-000_md-minimize', 'md'
        t.pipeline = NS()
        t.resource_manager = NS(charmmff_content=charmm)
        t.registered = []
        t.register = lambda data, key, artifact_type=None: t.registered.append((data, key))
        state = NS(psf=NS(path=psf, name=psf), pdb=NS(path=pdb or init, name=pdb or init))
        t.get_current_artifact = lambda key: state
        return t, state, psf, init


class TestTheTask(_TaskCase):

    def test_the_hold_is_on_by_default_for_minimize_and_nvt_only(self):
        assert MDTask._double_bond_hold_default('minimize')
        assert MDTask._double_bond_hold_default('nvt')
        assert not MDTask._double_bond_hold_default('npt')
        assert not MDTask._double_bond_hold_default('npgt'), 'a profile must not carry restraints'

    def test_the_file_is_merged_after_a_prior_one(self, tmp_path, monkeypatch, collection, charmm):
        t, state, *_ = self.make(tmp_path, monkeypatch, collection, charmm)
        Path('prior.extrabonds').write_text('dihedral 0 1 2 3 100.0 -33.0\n')
        name = t.make_double_bond_restraints_file(state, 100.0, prior='prior.extrabonds')
        lines = Path(name).read_text().splitlines()
        assert lines[0] == 'dihedral 0 1 2 3 100.0 -33.0'
        assert len(lines) == 2 and lines[1].endswith(' 100.0 0.0')
        assert t.registered == [('01-06-000_md-minimize-isomers', 'isomer_extrabonds')]

    def test_a_wrong_isomer_stops_the_build(self, tmp_path, monkeypatch, collection, charmm):
        t, state, psf, init = self.make(tmp_path, monkeypatch, collection, charmm)
        ch = DoubleBondChecker(psf, reference=ref(charmm))
        write_pdb(init, str(tmp_path / 'turned.pdb'),
                  turned(psf, ch.load_coords(init), *ch.terms[0][0][1:3]))
        state.pdb = NS(path=str(tmp_path / 'turned.pdb'), name=str(tmp_path / 'turned.pdb'))
        with pytest.raises(PestiferBuildError, match='C29=C210 is trans'):
            t.check_double_bonds()

    def test_warn_reports_instead(self, tmp_path, monkeypatch, collection, charmm, caplog):
        t, state, psf, init = self.make(tmp_path, monkeypatch, collection, charmm)
        ch = DoubleBondChecker(psf, reference=ref(charmm))
        write_pdb(init, 'turned.pdb', turned(psf, ch.load_coords(init), *ch.terms[0][0][1:3]))
        state.pdb = NS(path='turned.pdb', name='turned.pdb')
        with caplog.at_level(logging.WARNING, logger='pestifer.tasks.basetask'):
            t.check_double_bonds(stop=False)
        assert any('not the isomer CHARMM defines' in r.getMessage() for r in caplog.records)

    def test_a_right_structure_passes_and_says_so(self, tmp_path, monkeypatch, collection, charmm, caplog):
        t, *_ = self.make(tmp_path, monkeypatch, collection, charmm)
        with caplog.at_level(logging.INFO, logger='pestifer.tasks.basetask'):
            t.check_double_bonds()
        assert any('all 1 C=C double bonds with a CHARMM reference' in r.getMessage()
                   for r in caplog.records)

    @pytest.mark.parametrize('check,called', [(None, [True]), ('stop', [True]), ('warn', [False]),
                                              ('off', []), (False, [])])
    def test_do_runs_the_check_after_a_successful_stage(self, check, called):
        for rc in (0, 1):
            seen = []
            t = MDTask.__new__(MDTask)
            t.specs = {} if check is None else {'double_bonds': {'check': check}}
            t.namdrun = lambda rc=rc: rc
            t.report_stereo_strain = lambda: None
            t.track_stereochemistry = lambda **kw: None
            t.check_double_bonds = lambda stop=True: seen.append(stop)
            t.do()
            assert seen == (called if rc == 0 else [])
