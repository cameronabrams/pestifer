# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Nothing ships with a pierced ring in it.

A threaded ring is a topologically wrong molecule: once a bond sits on the wrong side of a
ring, neither minimization nor further dynamics puts it back.  It propagates through however
many microseconds of production follow the package, and every number computed from that run is
computed on the wrong chemistry.  **That is a failure whether or not the integrator complains**,
so this gate is not conditional on anything having crashed.

Measured 2026-10-16: every replica of an example-17 arm carried 2-6 protein-lipid ring
piercings at its post-embed hand-off, all invisible because every shipped example ran
``ring_check`` with ``segtypes: [lipid, glycan]`` -- and nothing checked the final structure at
all.  Three of those builds completed and were packaged.

No VMD, no NAMD: `ring_check` is monkeypatched, because what is under test is the DECISION, not
the geometry engine, which has its own tests.
"""
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from pestifer.core.errors import PestiferBuildError
from pestifer.tasks.terminate import TerminateTask


def _artifact(path):
    return SimpleNamespace(path=path, name=str(path))


def _task(state=SimpleNamespace(psf=_artifact('x.psf'), pdb=_artifact('x.pdb'),
                                xsc=_artifact('x.xsc'))):
    t = TerminateTask.__new__(TerminateTask)
    t.taskname = 'terminate'
    t.get_current_artifact = lambda key: state if key == 'state' else None
    return t


def _piercing(pe_seg='L1', pe_res='CHL1', pe_id='159', pr_seg='C', pr_res='LEU', pr_id='702'):
    return {'piercee': {'segname': pe_seg, 'resname': pe_res, 'resid': pe_id, 'segtype': 'lipid'},
            'piercer': {'segname': pr_seg, 'resname': pr_res, 'resid': pr_id, 'segtype': 'protein'}}


class TestTerminateRefusesToShipAPiercing(unittest.TestCase):

    def test_a_clean_structure_packages(self):
        with patch('pestifer.psfutil.psfring.ring_check', return_value=[]):
            _task()._refuse_to_ship_a_piercing()   # must not raise

    def test_a_piercing_stops_the_package(self):
        with patch('pestifer.psfutil.psfring.ring_check', return_value=[_piercing()]):
            with self.assertRaises(PestiferBuildError) as cm:
                _task()._refuse_to_ship_a_piercing()
        self.assertIn('refusing to package', str(cm.exception))

    def test_the_message_names_both_residues(self):
        """A build that fails at `terminate` has spent its whole runtime; the message has to be
        enough to act on without re-running anything."""
        with patch('pestifer.psfutil.psfring.ring_check', return_value=[_piercing()]):
            with self.assertRaises(PestiferBuildError) as cm:
                _task()._refuse_to_ship_a_piercing()
        msg = str(cm.exception)
        for token in ('CHL1', '159', 'LEU', '702'):
            self.assertIn(token, msg)
        self.assertIn('segtypes: [lipid, glycan, protein]', msg)

    def test_every_piercing_is_listed_not_just_the_first(self):
        many = [_piercing(pe_id=str(i)) for i in (159, 268, 199)]
        with patch('pestifer.psfutil.psfring.ring_check', return_value=many):
            with self.assertRaises(PestiferBuildError) as cm:
                _task()._refuse_to_ship_a_piercing()
        msg = str(cm.exception)
        for i in ('159', '268', '199'):
            self.assertIn(i, msg)
        self.assertIn('3 pierced ring(s)', msg)

    def test_it_checks_protein_rings_too(self):
        """`[lipid, glycan]` is what every example used and is what hid these."""
        with patch('pestifer.psfutil.psfring.ring_check', return_value=[]) as rc:
            _task()._refuse_to_ship_a_piercing()
        self.assertEqual(rc.call_args.kwargs['segtypes'], ['lipid', 'glycan', 'protein'])

    def test_an_unreadable_structure_does_not_lose_the_package(self):
        """This gate must not become a new way for a finished build to be thrown away: it only
        raises for a piercing it actually found."""
        with patch('pestifer.psfutil.psfring.ring_check', side_effect=RuntimeError('bad psf')):
            _task()._refuse_to_ship_a_piercing()   # must not raise

    def test_a_task_with_no_state_is_a_no_op(self):
        for state in (None,
                      SimpleNamespace(psf=None, pdb=_artifact('x.pdb'), xsc=None),
                      SimpleNamespace(psf=_artifact('x.psf'), pdb=None, xsc=None)):
            with self.subTest(state=state):
                _task(state)._refuse_to_ship_a_piercing()

    def test_a_system_with_no_xsc_is_still_checked(self):
        """A vacuum build has no periodic cell; it can still carry a piercing."""
        state = SimpleNamespace(psf=_artifact('x.psf'), pdb=_artifact('x.pdb'), xsc=None)
        with patch('pestifer.psfutil.psfring.ring_check', return_value=[]) as rc:
            _task(state)._refuse_to_ship_a_piercing()
        self.assertIsNone(rc.call_args.kwargs['xsc'])


class TestTheGateIsActuallyWiredIntoDo(unittest.TestCase):
    """Every test above calls `_refuse_to_ship_a_piercing` directly, so all of them stay green
    if the CALL is deleted from `TerminateTask.do`.  That is the defect this repository keeps a
    list of -- a negative control aimed at the helper rather than its call site -- and it was
    caught here by reverting the change in an export and noticing which tests did *not* go red.
    """

    def test_do_calls_the_gate(self):
        import inspect
        src = inspect.getsource(TerminateTask.do)
        self.assertIn('_refuse_to_ship_a_piercing', src,
                      'TerminateTask.do no longer calls the piercing gate')

    def test_it_is_called_before_the_package_is_written(self):
        """Raising after the tarball exists would leave a packaged structure on disk carrying
        the defect the gate refused."""
        import inspect
        src = inspect.getsource(TerminateTask.do)
        gate = src.index('_refuse_to_ship_a_piercing')
        for later in ('make_package', 'copy_state_to_basename'):
            if later in src:
                self.assertLess(gate, src.index(later),
                                f'the gate runs after {later}')
