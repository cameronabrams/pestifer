import os
import tempfile
import types
import unittest
from unittest import mock

from pestifer.tasks.solvate import SolvateTask
from pestifer.core.errors import PestiferError


def _make_task(repo, generate=True, in_ff=True):
    """A bare SolvateTask with just enough wiring to exercise _solvent_box_args / _solvent_box_entry."""
    t = SolvateTask.__new__(SolvateTask)
    cc = mock.Mock()
    cc.pdbrepository = repo
    cc.generate_missing_coordinates = generate
    cc.__contains__ = mock.Mock(return_value=in_ff)
    cc.charmmff_path = '/x/feb26'
    rm = mock.Mock()
    rm.charmmff_content = cc
    rm._charmmff_config = {'release': 'February2026'}
    t.resource_manager = rm
    return t


def _task_with_specs(specs):
    t = SolvateTask.__new__(SolvateTask)
    t.specs = specs
    return t


class TestIonizationPlan(unittest.TestCase):
    def test_neutralize_default_true(self):
        # with no specs, neutralize defaults True -> ions wanted, for water and non-water
        for solvent, is_water in (('TIP3', True), ('DMSO', False)):
            wants_ions, iw, neutralize = _task_with_specs({})._ionization_plan(solvent)
            self.assertTrue(wants_ions)
            self.assertTrue(neutralize)
            self.assertEqual(iw, is_water)

    def test_neutralize_false_no_salt_skips(self):
        for solvent in ('TIP3', 'DMSO'):
            wants_ions, _iw, neutralize = _task_with_specs({'neutralize': False})._ionization_plan(solvent)
            self.assertFalse(wants_ions)
            self.assertFalse(neutralize)

    def test_salt_forces_ions_even_without_neutralize(self):
        wants_ions, _iw, neutralize = _task_with_specs(
            {'neutralize': False, 'salt_con': 0.15})._ionization_plan('DMSO')
        self.assertTrue(wants_ions)
        self.assertFalse(neutralize)


def _box_entry(key='O1'):
    entry = mock.Mock()
    entry.is_box.return_value = True
    entry.get_box_psf.return_value = 'MEOH-box.psf'
    entry.get_box_pdb.return_value = 'MEOH-box.pdb'
    entry.get_box_edge.return_value = 20.34
    entry.get_key_atom.return_value = key
    return entry


class TestSolventBoxLookup(unittest.TestCase):
    def test_water_returns_no_entry_and_no_args(self):
        t = _make_task(repo=None)
        for w in ('TIP3', 'tip3', 'WATER', 'water', 'WAT', None):
            self.assertIsNone(t._solvent_box_entry(w))
            self.assertEqual(t._solvent_box_args(None, w), '')

    def test_non_water_box_entry_and_args(self):
        repo = mock.MagicMock()
        repo.__contains__.return_value = True
        repo.checkout.return_value = _box_entry()
        t = _make_task(repo)
        entry = t._solvent_box_entry('MEOH')
        args = t._solvent_box_args(entry, 'MEOH')
        self.assertIn('-spsf MEOH-box.psf', args)
        self.assertIn('-spdb MEOH-box.pdb', args)
        self.assertIn('-ws 20.34', args)
        self.assertIn('-ks "name O1"', args)

    def test_missing_solvent_generates_when_enabled(self):
        # default (generate_missing_coordinates True): a miss triggers on-demand generation,
        # then the freshly registered box entry is returned
        repo = mock.MagicMock()
        repo.__contains__.side_effect = [False, True]  # miss, then present after generation
        repo.checkout.return_value = _box_entry()
        t = _make_task(repo, generate=True)
        cc = t.resource_manager.charmmff_content
        with mock.patch('pestifer.charmmff.autocache.ensure_solvent_box',
                        return_value='/cache/feb26/solvent') as m:
            entry = t._solvent_box_entry('NOPE')
        m.assert_called_once_with('NOPE', cc)
        repo.add_resource.assert_called_once_with('/cache/feb26/solvent')
        self.assertTrue(entry.is_box())

    def test_missing_solvent_raises_when_generation_disabled(self):
        repo = mock.MagicMock()
        repo.__contains__.return_value = False
        with mock.patch('pestifer.charmmff.autocache.ensure_solvent_box') as m:
            with self.assertRaises(PestiferError) as ctx:
                _make_task(repo, generate=False)._solvent_box_entry('NOPE')
            m.assert_not_called()
        self.assertIn('generate_missing_coordinates', str(ctx.exception))

    def test_missing_solvent_not_in_forcefield_raises(self):
        # the live CC may not be residue-provisioned, so the FF-membership check happens inside
        # ensure_solvent_box (fresh RM); a not-in-FF solvent surfaces as a ValueError there,
        # which _generate_solvent_box wraps as a PestiferError
        repo = mock.MagicMock()
        repo.__contains__.return_value = False
        with mock.patch('pestifer.charmmff.autocache.ensure_solvent_box',
                        side_effect=ValueError('RESI NOPE is not defined in the CHARMM force field')):
            with self.assertRaises(PestiferError):
                _make_task(repo, generate=True)._solvent_box_entry('NOPE')

    def test_molecule_entry_rejected(self):
        entry = mock.Mock(); entry.is_box.return_value = False
        repo = mock.MagicMock()
        repo.__contains__.return_value = True
        repo.checkout.return_value = entry
        with self.assertRaises(PestiferError):
            _make_task(repo)._solvent_box_entry('SOD')

    def test_key_atom_omitted_when_absent(self):
        t = _make_task(repo=None)
        self.assertNotIn('-ks', t._solvent_box_args(_box_entry(key=None), 'X'))


class TestSolventTopology(unittest.TestCase):
    def test_writes_mass_records_and_resi(self):
        # solvate's isolated psfgen context needs the atom-type MASS records + the RESI block
        A = lambda name, t, m, el: types.SimpleNamespace(name=name, type=t, mass=m, element=el)
        topo = mock.Mock()
        topo.atoms = [A('CB', 'CT3', 12.011, 'C'), A('OG', 'OH1', 15.999, 'O'),
                      A('HG1', 'H', 1.008, 'H'), A('HB1', 'HA3', 1.008, 'H'), A('HB2', 'HA3', 1.008, 'H')]
        topo.to_file = lambda buf: buf.write('RESI MEOH  0.00\nATOM CB CT3 -0.04\n')
        cc = mock.Mock(); cc.get_resi.return_value = topo
        t = SolvateTask.__new__(SolvateTask)
        t.resource_manager = mock.Mock(); t.resource_manager.charmmff_content = cc
        out = os.path.join(tempfile.mkdtemp(), 'MEOH.rtf')
        t._write_solvent_topology('MEOH', out)
        s = open(out).read()
        self.assertIn('MASS', s)
        for typ in ('CT3', 'OH1', 'HA3'):        # each unique atom type gets one MASS record
            self.assertIn(typ, s)
        self.assertEqual(s.count('HA3'), 1)      # deduplicated (HB1/HB2 share HA3)
        self.assertIn('RESI MEOH', s)
        # without this directive, solvate's replica residues get bonds but no angles/dihedrals
        # and the molecules distort under MD
        self.assertIn('AUTOGENERATE ANGLES DIHEDRALS', s)
        self.assertTrue(s.rstrip().endswith('END'))


class TestIonCounts(unittest.TestCase):
    def _task(self):
        t = SolvateTask.__new__(SolvateTask)
        rm = mock.Mock()
        rm.charmmff_content.pdbrepository = None      # fall back to the built-in valence map
        t.resource_manager = rm
        return t

    def test_neutralize_positive_adds_anions(self):
        counts = self._task()._ion_counts(net_charge=3, box_volume_A3=1e5,
                                          cation='SOD', anion='CLA', salt_con=None)
        self.assertEqual(counts, {'CLA': 3})

    def test_neutralize_negative_adds_cations(self):
        counts = self._task()._ion_counts(net_charge=-2, box_volume_A3=1e5,
                                          cation='SOD', anion='CLA', salt_con=None)
        self.assertEqual(counts, {'SOD': 2})

    def test_divalent_cation_halves_count(self):
        counts = self._task()._ion_counts(net_charge=-4, box_volume_A3=1e5,
                                          cation='CAL', anion='CLA', salt_con=None)
        self.assertEqual(counts, {'CAL': 2})       # +2 each -> 2 ions neutralize -4

    def test_salt_plus_neutralize(self):
        # 0.15 M in a (60 A)^3 box: n_pairs = round(0.15 * 2.16e5 * 1e-27 * 6.022e23) ~ 20
        vol = 60.0 ** 3
        counts = self._task()._ion_counts(net_charge=3, box_volume_A3=vol,
                                          cation='SOD', anion='CLA', salt_con=0.15)
        n_pairs = round(0.15 * vol * 1e-27 * 6.02214076e23)
        self.assertEqual(counts['SOD'], n_pairs)
        self.assertEqual(counts['CLA'], n_pairs + 3)   # pairs + 3 neutralizing anions

    def test_neutral_system_no_ions(self):
        self.assertEqual(self._task()._ion_counts(0, 1e5, 'SOD', 'CLA', None), {})

    def test_neutralize_false_adds_salt_only(self):
        vol = 60.0 ** 3
        counts = self._task()._ion_counts(3, vol, 'SOD', 'CLA', 0.15, neutralize=False)
        n_pairs = round(0.15 * vol * 1e-27 * 6.02214076e23)
        self.assertEqual(counts, {'SOD': n_pairs, 'CLA': n_pairs})   # no neutralizing anion


class TestDryMembraneCore(unittest.TestCase):
    """``dry_membrane_core`` deletes water between the phosphate planes.  Pinned: it is written
    between ``solvate`` and ionization (so ions only ever see the water that stays), it is not
    written unless asked, and its count is read back -- a deletion that reports nothing must not
    read as zero deleted."""

    class _VT:
        def __init__(self):
            self.lines = []

        def newscript(self, *a, **k):
            pass

        def addline(self, line):
            self.lines.append(line)

        def writescript(self, *a, **k):
            pass

        def runscript(self, *a, **k):
            return 1      # stop do() right after the script is written

    def _script(self, specs):
        with tempfile.TemporaryDirectory() as d:
            xsc = os.path.join(d, 'cell.xsc')
            with open(xsc, 'w') as f:
                f.write('# NAMD extended system configuration\n'
                        '#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z\n'
                        '0 60 0 0 0 60 0 0 0 90 30 30 0\n')
            t = _task_with_specs(specs)
            t.basename = 'solv'
            t.next_basename = lambda *a, **k: None
            t.register = mock.Mock()
            state = types.SimpleNamespace(psf=types.SimpleNamespace(name='m.psf'),
                                          pdb=types.SimpleNamespace(name='m.pdb'),
                                          xsc=types.SimpleNamespace(name=xsc))
            t.get_current_artifact = lambda k: state
            t._solvent_box_entry = lambda solvent: None
            t._solvent_box_args = lambda entry, solvent: ''
            vt = self._VT()
            t.scripters = {'vmd': vt}
            t.do()
            return vt.lines

    def test_deletion_runs_after_solvate_and_before_ionization(self):
        lines = self._script({'dry_membrane_core': True, 'salt_con': 0.15})
        at = lambda word: next(i for i, ln in enumerate(lines) if ln.lstrip().startswith(word))
        self.assertLess(at('solvate '), at('foreach sr $_dry_sr'))
        self.assertLess(at('foreach sr $_dry_sr'), at('autoionize '))
        self.assertTrue(any('delatom' in ln for ln in lines))

    def test_off_by_default(self):
        lines = self._script({'salt_con': 0.15})
        self.assertFalse(any('delatom' in ln for ln in lines))

    def test_the_count_is_read_back_and_recorded(self):
        t = _task_with_specs({})
        t.record_outcome = mock.Mock()
        with tempfile.TemporaryDirectory() as d:
            log = os.path.join(d, 'solv.log')
            with open(log, 'w') as f:
                f.write('Info) dry_membrane_core: removed 37 water(s) between z -19.80 and 19.95\n')
            t._record_dry_membrane_core(log)
        t.record_outcome.assert_called_once_with(dry_membrane_core_removed=37)

    def test_no_report_is_a_warning_not_a_zero(self):
        t = _task_with_specs({})
        t.record_outcome = mock.Mock()
        with tempfile.TemporaryDirectory() as d:
            log = os.path.join(d, 'solv.log')
            open(log, 'w').close()
            with self.assertLogs('pestifer.tasks.solvate', level='WARNING'):
                t._record_dry_membrane_core(log)
        t.record_outcome.assert_not_called()


if __name__ == '__main__':
    unittest.main()
