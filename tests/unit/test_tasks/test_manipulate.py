# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Tests for ManipulateTask — 'align' and 'transfer_coords' coormods.

Unit tests (no VMD/NAMD needed):
  - TestWriteAlign:          write_align() TCL generation
  - TestTransferCoordsObj:   TransferCoords object construction and validation
  - TestWriteTransferCoords: write_transfer_coords() TCL generation

Integration tests (require VMD; run with --runslow):
  - TestManipulateAlignIntegration:         align a displaced system back to a reference
  - TestManipulateTransferCoordsIntegration: transfer coordinates with and without pre-alignment
"""
import math
import os
import shutil
import unittest
from pathlib import Path

import pytest

from pestifer.objs.align import Align
from pestifer.objs.rottrans import RotTrans
from pestifer.objs.transfer_coords import TransferCoords
from pestifer.scripters.vmd import VMDScripter

pytestmark = pytest.mark.needs_tools


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

class _MockScripter(VMDScripter):
    """Minimal stub that captures addline() calls without touching the filesystem."""
    def __init__(self, molid_varname='mCM'):
        self.lines = []
        self.molid_varname = molid_varname

    def addline(self, line):
        self.lines.append(line)


def _translate_pdb(src: Path, dst: Path, dx: float, dy: float, dz: float):
    """Write a copy of a PDB with all ATOM/HETATM coordinates shifted by (dx, dy, dz)."""
    with open(src) as fh:
        lines = fh.readlines()
    with open(dst, 'w') as fh:
        for line in lines:
            if line[:6] in ('ATOM  ', 'HETATM'):
                x = float(line[30:38]) + dx
                y = float(line[38:46]) + dy
                z = float(line[46:54]) + dz
                line = f'{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}'
            fh.write(line)


def _perturb_pdb(src: Path, dst: Path, dx: float, dy: float, dz: float,
                 resid_set: set[int]):
    """Write a copy of a PDB with coordinates of the given resids shifted by (dx, dy, dz)."""
    with open(src) as fh:
        lines = fh.readlines()
    with open(dst, 'w') as fh:
        for line in lines:
            if line[:6] in ('ATOM  ', 'HETATM'):
                try:
                    resid = int(line[22:26])
                except ValueError:
                    fh.write(line)
                    continue
                if resid in resid_set:
                    x = float(line[30:38]) + dx
                    y = float(line[38:46]) + dy
                    z = float(line[46:54]) + dz
                    line = f'{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}'
            fh.write(line)


def _pdb_coords(path: Path) -> list[tuple[float, float, float]]:
    """Return (x, y, z) for every ATOM/HETATM record in order."""
    coords = []
    with open(path) as fh:
        for line in fh:
            if line[:6] in ('ATOM  ', 'HETATM'):
                coords.append((
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ))
    return coords


def _rmsd(coords_a, coords_b) -> float:
    assert len(coords_a) == len(coords_b), 'atom count mismatch'
    n = len(coords_a)
    total = sum(
        (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
        for (ax, ay, az), (bx, by, bz) in zip(coords_a, coords_b)
    )
    return math.sqrt(total / n)


# ---------------------------------------------------------------------------
# Unit tests: write_align() TCL generation
# ---------------------------------------------------------------------------

class TestWriteAlign(unittest.TestCase):

    def _lines(self, **kwargs) -> list[str]:
        vm = _MockScripter()
        a = Align(**kwargs)
        vm.write_align(a)
        return vm.lines

    def test_pdb_only_load(self):
        lines = self._lines(ref_pdb='ref.pdb')
        self.assertTrue(any('mol load pdb {ref.pdb}' in l for l in lines))
        self.assertFalse(any('mol load psf' in l for l in lines))

    def test_psf_pdb_load(self):
        lines = self._lines(ref_pdb='ref.pdb', ref_psf='ref.psf')
        self.assertTrue(any('mol load psf {ref.psf} pdb {ref.pdb}' in l for l in lines))

    def test_ref_sel_defaults_to_mobile_sel(self):
        lines = self._lines(ref_pdb='ref.pdb', mobile_sel='backbone')
        fit_line = next(l for l in lines if '_ref_fit' in l and 'atomselect' in l)
        self.assertIn('backbone', fit_line)

    def test_custom_ref_sel(self):
        lines = self._lines(ref_pdb='ref.pdb', mobile_sel='backbone', ref_sel='name CA')
        mob_line = next(l for l in lines if '_mob_fit' in l and 'atomselect' in l)
        ref_line = next(l for l in lines if '_ref_fit' in l and 'atomselect' in l)
        self.assertIn('backbone', mob_line)
        self.assertIn('name CA', ref_line)

    def test_apply_to_default_all(self):
        lines = self._lines(ref_pdb='ref.pdb')
        mover_line = next(l for l in lines if '_mover' in l and 'atomselect' in l)
        self.assertIn('"all"', mover_line)

    def test_custom_apply_to(self):
        lines = self._lines(ref_pdb='ref.pdb', apply_to='protein')
        mover_line = next(l for l in lines if '_mover' in l and 'atomselect' in l)
        self.assertIn('protein', mover_line)

    def test_congruency_check_present(self):
        lines = self._lines(ref_pdb='ref.pdb')
        self.assertTrue(any('_n_mob' in l for l in lines))
        self.assertTrue(any('_n_ref' in l for l in lines))
        self.assertTrue(any('exit 1' in l for l in lines))

    def test_congruency_error_message_includes_selections(self):
        lines = self._lines(ref_pdb='ref.pdb', mobile_sel='backbone', ref_sel='name CA')
        # emitted with the PESTIFER-ERROR marker so the manipulate task's log scan turns a
        # congruency mismatch into a hard error (VMD's exit 1 alone returns process code 0)
        error_line = next((l for l in lines if 'PESTIFER-ERROR' in l and 'align' in l), None)
        self.assertIsNotNone(error_line)
        self.assertIn('backbone', error_line)
        self.assertIn('name CA', error_line)

    def test_ref_mol_deleted_after_fit(self):
        lines = self._lines(ref_pdb='ref.pdb')
        fit_idx = next(i for i, l in enumerate(lines) if 'measure fit' in l)
        delete_idx = next(i for i, l in enumerate(lines) if 'mol delete' in l)
        self.assertGreater(delete_idx, fit_idx)

    def test_move_applied_after_fit(self):
        lines = self._lines(ref_pdb='ref.pdb')
        fit_idx = next(i for i, l in enumerate(lines) if 'measure fit' in l)
        move_idx = next(i for i, l in enumerate(lines) if '$_mover move' in l)
        self.assertGreater(move_idx, fit_idx)

    def test_molid_varname_used(self):
        vm = _MockScripter(molid_varname='myMol')
        vm.write_align(Align(ref_pdb='ref.pdb'))
        self.assertTrue(any('$myMol' in l for l in vm.lines))


# ---------------------------------------------------------------------------
# Unit tests: write_rottrans() TCL generation (transrot mods)
# ---------------------------------------------------------------------------

class TestWriteRotTrans(unittest.TestCase):

    def _lines(self, rt: RotTrans, molid=None) -> list[str]:
        vm = _MockScripter()
        vm.write_rottrans(rt, molid=molid) if molid else vm.write_rottrans(rt)
        return vm.lines

    def test_rotation_emitted_about_com(self):
        lines = self._lines(RotTrans(movetype='ROT', axis='z', angle=-9.654))
        self.assertTrue(any('measure center $mover weight mass' in l for l in lines))
        # rotate about the COM (`trans center`), not about the global origin (`trans origin`)
        self.assertTrue(any('trans center $COM axis z -9.654' in l for l in lines),
                        f'rotation not about COM: {lines}')
        self.assertFalse(any('trans origin' in l for l in lines))

    def test_translation_emitted(self):
        lines = self._lines(RotTrans(movetype='TRANS', x=1.0, y=2.0, z=3.0))
        self.assertTrue(any('moveby [list 1.0 2.0 3.0]' in l for l in lines))

    def test_molid_kwarg_accepted(self):
        # regression: ManipulateTask.coormods calls write_rottrans(obj, molid='mCM');
        # the method must accept the kwarg (it previously did not -> TypeError)
        vm = _MockScripter(molid_varname='mCM')
        vm.write_rottrans(RotTrans(movetype='ROT', axis='y', angle=-168.646), molid='mCM')
        self.assertTrue(any('$mCM' in l for l in vm.lines))

    def test_default_sel_is_all_no_guard(self):
        # the default (whole-system) transform selects "all" and needs no disconnection guard
        lines = self._lines(RotTrans(movetype='ROT', axis='z', angle=90.0))
        mover_line = next(l for l in lines if 'set mover' in l)
        self.assertIn('"all"', mover_line)
        self.assertFalse(any('getbonds' in l for l in lines))
        self.assertFalse(any('not fully disconnected' in l for l in lines))

    def test_custom_sel_drives_atomselect(self):
        lines = self._lines(RotTrans(movetype='TRANS', x=1.0, y=2.0, z=3.0,
                                     sel='segid GLYC'))
        mover_line = next(l for l in lines if 'set mover' in l)
        self.assertIn('"segid GLYC"', mover_line)

    def test_custom_sel_emits_disconnection_guard(self):
        lines = self._lines(RotTrans(movetype='ROT', axis='z', angle=90.0,
                                     sel='fragment 3'))
        self.assertTrue(any('getbonds' in l for l in lines))
        self.assertTrue(any('dict exists $_in_sel' in l for l in lines))
        guard_err = next((l for l in lines if 'not fully disconnected' in l), None)
        self.assertIsNotNone(guard_err)
        self.assertIn('fragment 3', guard_err)
        # guard must hard-error
        self.assertTrue(any('exit 1' in l for l in lines))

    def test_guard_precedes_transform(self):
        lines = self._lines(RotTrans(movetype='ROT', axis='z', angle=90.0, sel='fragment 3'))
        guard_idx = next(i for i, l in enumerate(lines) if 'getbonds' in l)
        move_idx = next(i for i, l in enumerate(lines) if '$mover move' in l)
        self.assertLess(guard_idx, move_idx)

    def test_align_literal_vectors_minimal_rotation(self):
        lines = self._lines(RotTrans(movetype='ALIGN', source=[1.0, 0.0, 0.0],
                                     target=[0.0, 0.0, 1.0]))
        text = '\n'.join(lines)
        # source/target set as literal vectors
        self.assertIn('set _src [list 1.0 0.0 0.0]', text)
        self.assertIn('set _tgt [list 0.0 0.0 1.0]', text)
        # minimal rotation: axis = source x target, angle = atan2(|cross|, dot), about the COM
        self.assertIn('veccross $_src $_tgt', text)
        self.assertIn('atan2($_sin, $_cos)', text)
        self.assertTrue(any('trans center $COM axis [vecnorm $_cross]' in l for l in lines))
        # degenerate antiparallel branch present
        self.assertTrue(any('180.0' in l for l in lines))

    def test_align_selection_pair_source(self):
        lines = self._lines(RotTrans(movetype='ALIGN',
                                     source=['resid 1 and name CA', 'resid 30 and name CA'],
                                     target=[0.0, 0.0, 1.0]))
        text = '\n'.join(lines)
        # source vector built from the difference of two selection centers
        self.assertTrue(any('resid 1 and name CA' in l and 'atomselect' in l for l in lines))
        self.assertTrue(any('resid 30 and name CA' in l and 'atomselect' in l for l in lines))
        self.assertIn('vecsub [measure center $_vb weight mass] [measure center $_va weight mass]', text)
        # temporary selections are cleaned up
        self.assertTrue(any('$_va delete' in l for l in lines))
        self.assertTrue(any('$_vb delete' in l for l in lines))

    def test_align_respects_sel_and_guard(self):
        lines = self._lines(RotTrans(movetype='ALIGN', source=[1.0, 0.0, 0.0],
                                     target=[0.0, 0.0, 1.0], sel='fragment 2'))
        mover_line = next(l for l in lines if 'set mover' in l)
        self.assertIn('"fragment 2"', mover_line)
        self.assertTrue(any('getbonds' in l for l in lines))


class TestRotTransObj(unittest.TestCase):

    def test_sel_optional_defaults_none(self):
        rt = RotTrans(movetype='ROT', axis='z', angle=90.0)
        self.assertIsNone(rt.sel)

    def test_sel_settable_via_dict(self):
        rt = RotTrans(movetype='TRANS', x=1.0, y=0.0, z=0.0, sel='segid QQQ')
        self.assertEqual(rt.sel, 'segid QQQ')

    def test_yaml_header_and_alias(self):
        self.assertEqual(RotTrans._yaml_header, 'transrot')
        self.assertIn('rottrans', RotTrans._yaml_aliases)

    def test_rottrans_alias_maps_to_rottrans_class(self):
        from pestifer.core.objmanager import ObjManager
        self.assertIs(ObjManager._obj_classes_byYAML['rottrans'], RotTrans)
        self.assertIs(ObjManager._obj_classes_byYAML['transrot'], RotTrans)

    def test_alias_ingests_under_canonical_header(self):
        # ingesting a 'rottrans' mod must store it under the canonical 'transrot' key so the
        # manipulate task (which dispatches on 'transrot') still finds it
        from pestifer.core.objmanager import ObjManager
        om = ObjManager()
        om.ingest({'rottrans': ['ROT,z,45.0']})
        coord = om.get('coord', {})
        self.assertIn('transrot', coord)
        self.assertEqual(len(coord['transrot']), 1)
        self.assertIsInstance(coord['transrot'][0], RotTrans)

    def test_align_literal_vectors_ok(self):
        rt = RotTrans(movetype='ALIGN', source=[1, 0, 0], target=[0, 0, 1])
        self.assertEqual(rt.source, [1, 0, 0])
        self.assertEqual(rt.target, [0, 0, 1])

    def test_align_selection_pair_ok(self):
        rt = RotTrans(movetype='ALIGN', source=['resid 1', 'resid 30'], target=[0, 0, 1])
        self.assertEqual(rt.source, ['resid 1', 'resid 30'])

    def test_align_requires_source_and_target(self):
        with self.assertRaises(Exception):
            RotTrans(movetype='ALIGN', source=[0, 0, 1])

    def test_align_rejects_malformed_vector(self):
        # 3 strings is neither a literal 3-vector nor a selection pair
        with self.assertRaises(Exception):
            RotTrans(movetype='ALIGN', source=['a', 'b', 'c'], target=[0, 0, 1])
        # a 2-number source is not a valid literal vector or selection pair
        with self.assertRaises(Exception):
            RotTrans(movetype='ALIGN', source=[0, 1], target=[0, 0, 1])


# ---------------------------------------------------------------------------
# Integration test: full VMD run  (--runslow)
# ---------------------------------------------------------------------------

_FIXTURES = (Path(__file__).parent / 'fixtures' / 'continuation_inputs').resolve()
_PSF = 'my_6pti.psf'
_PDB = 'my_6pti.pdb'
_COOR = 'my_6pti.coor'
_RMSD_THRESHOLD = 0.1   # Å; pure-translation alignment should be essentially exact


class TestManipulateAlignIntegration(unittest.TestCase):

    def setUp(self):
        from pestifer.core.config import Config
        from pestifer.core.controller import Controller
        self.controller = Controller().configure(Config().configure_new(), terminate=False)

    def _setup_subdir(self, name: str, dx: float, dy: float, dz: float):
        """Create a clean working subdirectory, symlink PSF/COOR, write displaced PDB."""
        if os.path.exists(name):
            shutil.rmtree(name)
        os.makedirs(name)
        os.symlink(_FIXTURES / _PSF,  Path(name) / _PSF)
        os.symlink(_FIXTURES / _COOR, Path(name) / _COOR)
        # Reference: original coordinates
        os.symlink(_FIXTURES / _PDB,  Path(name) / 'ref.pdb')
        # Mobile: displaced by (dx, dy, dz)
        _translate_pdb(_FIXTURES / _PDB, Path(name) / _PDB, dx, dy, dz)
        os.chdir(name)

    def _run_align(self):
        task_list = [
            {'continuation': {
                'psf': _PSF,
                'pdb': _PDB,
                'coor': _COOR,
            }},
            {'manipulate': {
                'mods': {
                    'align': [{'ref_pdb': 'ref.pdb'}],
                },
            }},
        ]
        self.controller.reconfigure_tasks(task_list)
        self.controller.do_tasks()
        return self.controller.tasks[1]

    def test_align_pure_translation(self):
        """Aligning a system displaced by a pure translation recovers the original coordinates."""
        self._setup_subdir('__test_align_trans', dx=25.0, dy=-10.0, dz=40.0)
        task = self._run_align()
        from pestifer.core.artifacts import StateArtifacts
        state: StateArtifacts = task.get_current_artifact('state')
        self.assertTrue(state.pdb.exists(), 'output PDB not found after manipulate/align')
        aligned_coords = _pdb_coords(Path(state.pdb.name))
        ref_coords = _pdb_coords(Path('ref.pdb'))
        rmsd = _rmsd(aligned_coords, ref_coords)
        self.assertLess(rmsd, _RMSD_THRESHOLD,
                        f'RMSD {rmsd:.4f} Å after alignment exceeds threshold {_RMSD_THRESHOLD} Å')


class TestManipulateRotTransIntegration(unittest.TestCase):

    def setUp(self):
        from pestifer.core.config import Config
        from pestifer.core.controller import Controller
        self.controller = Controller().configure(Config().configure_new(), terminate=False)

    def test_transrot_actually_rotates(self):
        """A transrot ROT mod must move the structure (regression: coormods checked the wrong
        objtype key and emitted no rotation, leaving coordinates unchanged)."""
        name = '__test_transrot'
        if os.path.exists(name):
            shutil.rmtree(name)
        os.makedirs(name)
        os.symlink(_FIXTURES / _PSF, Path(name) / _PSF)
        os.symlink(_FIXTURES / _COOR, Path(name) / _COOR)
        os.symlink(_FIXTURES / _PDB, Path(name) / _PDB)
        os.chdir(name)
        task_list = [
            {'continuation': {'psf': _PSF, 'pdb': _PDB, 'coor': _COOR}},
            {'manipulate': {'mods': {'transrot': ['ROT,z,90.0', 'ROT,y,45.0']}}},
        ]
        self.controller.reconfigure_tasks(task_list)
        self.controller.do_tasks()
        from pestifer.core.artifacts import StateArtifacts
        state: StateArtifacts = self.controller.tasks[1].get_current_artifact('state')
        self.assertTrue(state.pdb.exists())
        out = _pdb_coords(Path(state.pdb.name))
        orig = _pdb_coords(_FIXTURES / _PDB)
        # the rotation moved the structure substantially (not the ~0 RMSD of a no-op)
        self.assertGreater(_rmsd(out, orig), 1.0)

        def centroid(c):
            n = len(c)
            return (sum(p[0] for p in c) / n, sum(p[1] for p in c) / n, sum(p[2] for p in c) / n)

        def rgyr(c):
            cx, cy, cz = centroid(c)
            return math.sqrt(sum((p[0]-cx)**2 + (p[1]-cy)**2 + (p[2]-cz)**2 for p in c) / len(c))

        # rigid: radius of gyration preserved (true for rotation about any point)
        self.assertAlmostEqual(rgyr(out), rgyr(orig), places=2)
        # and about the center of mass: the centroid barely moves.  The old bug rotated about
        # (and collapsed the COM onto) the global origin, which would fling the centroid ~10 Å
        # away from its original position.
        co, cr = centroid(out), centroid(orig)
        self.assertLess(math.dist(co, cr), 3.0,
                        f'centroid moved {math.dist(co, cr):.1f} Å -- not a COM-centered rotation')

    def test_transrot_connected_sel_hard_errors(self):
        """A transrot whose selection is bonded to the rest of the system must hard-error at
        VMD runtime (the disconnection guard), leaving the manipulate task non-zero."""
        name = '__test_transrot_guard'
        if os.path.exists(name):
            shutil.rmtree(name)
        os.makedirs(name)
        os.symlink(_FIXTURES / _PSF, Path(name) / _PSF)
        os.symlink(_FIXTURES / _COOR, Path(name) / _COOR)
        os.symlink(_FIXTURES / _PDB, Path(name) / _PDB)
        os.chdir(name)
        # 'resid 1 to 5' is covalently bonded to resid 6 -> a bond crosses the selection boundary
        task_list = [
            {'continuation': {'psf': _PSF, 'pdb': _PDB, 'coor': _COOR}},
            {'manipulate': {'mods': {'transrot': [
                {'movetype': 'TRANS', 'x': 10.0, 'y': 0.0, 'z': 0.0, 'sel': 'resid 1 to 5'}]}}},
        ]
        self.controller.reconfigure_tasks(task_list)
        self.controller.do_tasks()
        self.assertNotEqual(self.controller.tasks[1].result, 0,
                            'manipulate should fail when transrot selection is not disconnected')

    def test_transrot_align_vector_onto_z(self):
        """A transrot ALIGN carries a source vector (CA of resid 1 -> CA of resid 30) onto the
        target z-axis; after the move that vector must be parallel to [0,0,1]."""
        name = '__test_transrot_align'
        if os.path.exists(name):
            shutil.rmtree(name)
        os.makedirs(name)
        os.symlink(_FIXTURES / _PSF, Path(name) / _PSF)
        os.symlink(_FIXTURES / _COOR, Path(name) / _COOR)
        os.symlink(_FIXTURES / _PDB, Path(name) / _PDB)
        os.chdir(name)
        task_list = [
            {'continuation': {'psf': _PSF, 'pdb': _PDB, 'coor': _COOR}},
            {'manipulate': {'mods': {'transrot': [
                {'movetype': 'ALIGN',
                 'source': ['resid 1 and name CA', 'resid 30 and name CA'],
                 'target': [0.0, 0.0, 1.0]}]}}},
        ]
        self.controller.reconfigure_tasks(task_list)
        self.controller.do_tasks()
        self.assertEqual(self.controller.tasks[1].result, 0)
        from pestifer.core.artifacts import StateArtifacts
        state: StateArtifacts = self.controller.tasks[1].get_current_artifact('state')
        self.assertTrue(state.pdb.exists())

        def ca(path, resid):
            with open(path) as fh:
                for line in fh:
                    if line[:6] in ('ATOM  ', 'HETATM') and line[12:16].strip() == 'CA':
                        try:
                            if int(line[22:26]) == resid:
                                return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                        except ValueError:
                            pass
            raise AssertionError(f'CA of resid {resid} not found in {path}')

        a = ca(Path(state.pdb.name), 1)
        b = ca(Path(state.pdb.name), 30)
        v = (b[0]-a[0], b[1]-a[1], b[2]-a[2])
        n = math.sqrt(sum(c*c for c in v))
        # the source vector is now (anti)parallel to z: |z-component|/|v| ~ 1, x/y ~ 0
        self.assertAlmostEqual(abs(v[2])/n, 1.0, places=3,
                               msg=f'source vector not aligned to z after ALIGN: {v}')


# ---------------------------------------------------------------------------
# Unit tests: TransferCoords object
# ---------------------------------------------------------------------------

class TestTransferCoordsObj(unittest.TestCase):

    def test_required_fields(self):
        tc = TransferCoords(donor_pdb='d.pdb', donor_sel='chain A', mobile_sel='chain B')
        self.assertEqual(tc.donor_pdb, 'd.pdb')
        self.assertIsNone(tc.donor_psf)
        self.assertEqual(tc.donor_sel, 'chain A')
        self.assertEqual(tc.mobile_sel, 'chain B')
        self.assertIsNone(tc.align_donor_sel)
        self.assertIsNone(tc.align_mobile_sel)
        self.assertFalse(tc.pre_align)

    def test_with_psf(self):
        tc = TransferCoords(donor_pdb='d.pdb', donor_psf='d.psf',
                            donor_sel='all', mobile_sel='all')
        self.assertEqual(tc.donor_psf, 'd.psf')

    def test_with_pre_alignment(self):
        tc = TransferCoords(
            donor_pdb='d.pdb', donor_sel='chain A and resid 10 to 20',
            mobile_sel='chain B and resid 10 to 20',
            align_donor_sel='backbone', align_mobile_sel='chain B and backbone',
        )
        self.assertTrue(tc.pre_align)
        self.assertEqual(tc.align_donor_sel, 'backbone')
        self.assertEqual(tc.align_mobile_sel, 'chain B and backbone')

    def test_unpaired_align_donor_raises(self):
        with self.assertRaises(Exception):
            TransferCoords(donor_pdb='d.pdb', donor_sel='all', mobile_sel='all',
                           align_donor_sel='backbone')

    def test_unpaired_align_mobile_raises(self):
        with self.assertRaises(Exception):
            TransferCoords(donor_pdb='d.pdb', donor_sel='all', mobile_sel='all',
                           align_mobile_sel='backbone')

    def test_yaml_header_and_category(self):
        tc = TransferCoords(donor_pdb='d.pdb', donor_sel='all', mobile_sel='all')
        self.assertEqual(tc._yaml_header, 'transfer_coords')
        self.assertEqual(tc._objcat, 'coord')


# ---------------------------------------------------------------------------
# Unit tests: write_transfer_coords() TCL generation
# ---------------------------------------------------------------------------

class TestWriteTransferCoords(unittest.TestCase):

    def _lines(self, **kwargs) -> list[str]:
        vm = _MockScripter()
        tc = TransferCoords(**kwargs)
        vm.write_transfer_coords(tc)
        return vm.lines

    def test_pdb_only_load(self):
        lines = self._lines(donor_pdb='d.pdb', donor_sel='all', mobile_sel='all')
        self.assertTrue(any('mol load pdb {d.pdb}' in l for l in lines))
        self.assertFalse(any('mol load psf' in l for l in lines))

    def test_psf_pdb_load(self):
        lines = self._lines(donor_pdb='d.pdb', donor_psf='d.psf',
                            donor_sel='all', mobile_sel='all')
        self.assertTrue(any('mol load psf {d.psf} pdb {d.pdb}' in l for l in lines))

    def test_transfer_set_xyz(self):
        lines = self._lines(donor_pdb='d.pdb', donor_sel='chain A', mobile_sel='chain B')
        self.assertTrue(any('set {x y z}' in l and '_donor_sel' in l and '_mob_sel' in l
                            for l in lines))

    def test_transfer_congruency_check(self):
        lines = self._lines(donor_pdb='d.pdb', donor_sel='chain A', mobile_sel='chain B')
        self.assertTrue(any('_n_donor' in l for l in lines))
        self.assertTrue(any('_n_mob' in l for l in lines))
        self.assertTrue(any('exit 1' in l for l in lines))

    def test_transfer_error_message_includes_selections(self):
        lines = self._lines(donor_pdb='d.pdb', donor_sel='chain A', mobile_sel='chain B')
        error_line = next((l for l in lines if 'PESTIFER-ERROR' in l and 'transfer_coords' in l
                           and 'align' not in l), None)
        self.assertIsNotNone(error_line)
        self.assertIn('chain A', error_line)
        self.assertIn('chain B', error_line)

    def test_donor_mol_deleted(self):
        lines = self._lines(donor_pdb='d.pdb', donor_sel='all', mobile_sel='all')
        self.assertTrue(any('mol delete $_donor_mol' in l for l in lines))

    def test_no_alignment_block_when_not_requested(self):
        lines = self._lines(donor_pdb='d.pdb', donor_sel='all', mobile_sel='all')
        self.assertFalse(any('measure fit' in l for l in lines))
        self.assertFalse(any('_aln_' in l for l in lines))

    def test_alignment_block_present_when_requested(self):
        lines = self._lines(
            donor_pdb='d.pdb', donor_sel='resid 10 to 20', mobile_sel='resid 10 to 20',
            align_donor_sel='backbone', align_mobile_sel='protein and backbone',
        )
        self.assertTrue(any('measure fit' in l for l in lines))
        self.assertTrue(any('_aln_donor' in l for l in lines))
        self.assertTrue(any('_aln_mob' in l for l in lines))

    def test_alignment_congruency_check_present(self):
        lines = self._lines(
            donor_pdb='d.pdb', donor_sel='all', mobile_sel='all',
            align_donor_sel='backbone', align_mobile_sel='backbone',
        )
        self.assertTrue(any('_n_aln_donor' in l for l in lines))
        self.assertTrue(any('_n_aln_mob' in l for l in lines))

    def test_alignment_applied_to_all_donor_atoms(self):
        lines = self._lines(
            donor_pdb='d.pdb', donor_sel='all', mobile_sel='all',
            align_donor_sel='backbone', align_mobile_sel='backbone',
        )
        self.assertTrue(any('_aln_all' in l and 'move' in l for l in lines))

    def test_alignment_before_transfer(self):
        lines = self._lines(
            donor_pdb='d.pdb', donor_sel='all', mobile_sel='all',
            align_donor_sel='backbone', align_mobile_sel='backbone',
        )
        fit_idx   = next(i for i, l in enumerate(lines) if 'measure fit' in l)
        xfer_idx  = next(i for i, l in enumerate(lines) if 'set {x y z}' in l)
        self.assertLess(fit_idx, xfer_idx)

    def test_molid_varname_used(self):
        vm = _MockScripter(molid_varname='myMol')
        vm.write_transfer_coords(
            TransferCoords(donor_pdb='d.pdb', donor_sel='all', mobile_sel='all'))
        self.assertTrue(any('$myMol' in l for l in vm.lines))


# ---------------------------------------------------------------------------
# Integration tests: transfer_coords  (--runslow)
# ---------------------------------------------------------------------------

# Resids 1–5 of 6PTI are used as the "perturbed subset" in tests below.
_TRANSFER_RESIDS = {1, 2, 3, 4, 5}
_PERTURB_DX = 8.0   # Å shift applied to the subset in the donor


class TestManipulateTransferCoordsIntegration(unittest.TestCase):

    def setUp(self):
        from pestifer.core.config import Config
        from pestifer.core.controller import Controller
        self.controller = Controller().configure(Config().configure_new(), terminate=False)

    def _setup_subdir(self, name: str):
        if os.path.exists(name):
            shutil.rmtree(name)
        os.makedirs(name)
        os.symlink(_FIXTURES / _PSF,  Path(name) / _PSF)
        os.symlink(_FIXTURES / _COOR, Path(name) / _COOR)
        os.symlink(_FIXTURES / _PDB,  Path(name) / _PDB)
        os.chdir(name)

    def _run_transfer(self, donor_pdb: str, extra_specs: dict = {}):
        specs = {
            'donor_pdb': donor_pdb,
            'donor_sel': f'resid {min(_TRANSFER_RESIDS)} to {max(_TRANSFER_RESIDS)}',
            'mobile_sel': f'resid {min(_TRANSFER_RESIDS)} to {max(_TRANSFER_RESIDS)}',
            **extra_specs,
        }
        task_list = [
            {'continuation': {'psf': _PSF, 'pdb': _PDB, 'coor': _COOR}},
            {'manipulate': {'mods': {'transfer_coords': [specs]}}},
        ]
        self.controller.reconfigure_tasks(task_list)
        self.controller.do_tasks()
        return self.controller.tasks[1]

    def test_transfer_without_prealign(self):
        """Coordinates of the perturbed subset in the donor appear in the pipeline output."""
        self._setup_subdir('__test_transfer_no_align')
        donor = Path('donor_perturbed.pdb')
        _perturb_pdb(_FIXTURES / _PDB, donor, dx=_PERTURB_DX, dy=0.0, dz=0.0,
                     resid_set=_TRANSFER_RESIDS)
        task = self._run_transfer(str(donor))
        from pestifer.core.artifacts import StateArtifacts
        state: StateArtifacts = task.get_current_artifact('state')
        self.assertTrue(state.pdb.exists())
        out_coords  = _pdb_coords(Path(state.pdb.name))
        orig_coords = _pdb_coords(_FIXTURES / _PDB)
        donor_coords = _pdb_coords(donor)
        # Identify ATOM records belonging to the transfer subset by position in file
        # (same order as VMD reads them)
        transfer_mask = []
        with open(_FIXTURES / _PDB) as fh:
            for line in fh:
                if line[:6] in ('ATOM  ', 'HETATM'):
                    try:
                        resid = int(line[22:26])
                    except ValueError:
                        resid = -1
                    transfer_mask.append(resid in _TRANSFER_RESIDS)
        # Transferred atoms should match donor, non-transferred should match original
        for i, (o, d, orig, in_subset) in enumerate(
                zip(out_coords, donor_coords, orig_coords, transfer_mask)):
            if in_subset:
                self.assertAlmostEqual(o[0], d[0], places=2,
                    msg=f'atom {i}: x after transfer {o[0]:.3f} != donor {d[0]:.3f}')
            else:
                self.assertAlmostEqual(o[0], orig[0], places=2,
                    msg=f'atom {i}: x unchanged {o[0]:.3f} != original {orig[0]:.3f}')

    def test_transfer_with_prealign(self):
        """Pre-aligning a displaced donor before transfer recovers near-original coordinates."""
        self._setup_subdir('__test_transfer_with_align')
        # Donor is the same structure translated far away — pre-alignment must bring
        # the subset back into the pipeline frame before transfer.
        donor = Path('donor_displaced.pdb')
        _translate_pdb(_FIXTURES / _PDB, donor, dx=50.0, dy=0.0, dz=0.0)
        task = self._run_transfer(str(donor), extra_specs={
            'align_donor_sel': 'backbone',
            'align_mobile_sel': 'backbone',
        })
        from pestifer.core.artifacts import StateArtifacts
        state: StateArtifacts = task.get_current_artifact('state')
        self.assertTrue(state.pdb.exists())
        out_coords  = _pdb_coords(Path(state.pdb.name))
        orig_coords = _pdb_coords(_FIXTURES / _PDB)
        rmsd = _rmsd(out_coords, orig_coords)
        self.assertLess(rmsd, _RMSD_THRESHOLD,
            f'RMSD {rmsd:.4f} Å after pre-aligned transfer exceeds threshold')

    def test_transfer_count_mismatch_hard_errors(self):
        """A donor_sel/mobile_sel atom-count mismatch must hard-error the manipulate task.
        (Regression: the Tcl congruency check `exit 1`s, but VMD returns process code 0, so the
        failure used to be silent -- no output pdb -- and only surfaced at the next task.)"""
        self._setup_subdir('__test_transfer_mismatch')
        task = self._run_transfer(str(_FIXTURES / _PDB), extra_specs={
            'donor_sel': 'resid 1 to 3', 'mobile_sel': 'resid 1 to 8'})
        self.assertNotEqual(task.result, 0,
                            'transfer_coords count mismatch should hard-error the task')
