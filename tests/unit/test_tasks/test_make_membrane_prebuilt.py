# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""``make_membrane_system`` with a ``bilayer: prebuilt:`` block, which skips construction.

A build can supply its own already-equilibrated bilayer instead of having pestifer construct one
-- the restart case: take the membrane a previous run settled and embed into it.  The branch is in
the schema and in ``provision``, and until 2026-09-29 it raised ``AttributeError: 'dict' object has
no attribute 'xsc'`` on the first line that used it, for every config that took it.

It survived because nothing ran it.  ``test_make_membrane_decisions`` covers
``using_prebuilt_bilayer`` as a *routing flag* on a mocked task, so it never reaches ``provision``,
and **no bundled example sets** ``prebuilt``, so no sweep touched it either.  These tests exercise
``provision`` itself, which is the only place the defect lived.
"""
import unittest
from unittest import mock

from pestifer.tasks.make_membrane_system import MakeMembraneSystemTask


class _Recorder:
    """Stands in for the pipeline: records what was registered and hands back the artifact.

    The whole bug was the difference between the dict passed IN and the artifact handed BACK, so
    this fake has to honor that distinction or it cannot reproduce anything.
    """

    def __init__(self):
        self.registered = {}

    def register(self, data, key, requestor, artifact_type=None, **kw):
        artifact = mock.Mock()
        artifact.data = data
        artifact.psf.name = data['psf']
        artifact.pdb.name = data['pdb']
        artifact.xsc.path = data['xsc']
        self.registered[key] = artifact
        return artifact


def _task(prebuilt):
    t = mock.Mock(spec=MakeMembraneSystemTask)
    t.specs = {'bilayer': {'prebuilt': prebuilt}}
    t.provisions = {}
    t.resource_manager = mock.Mock()
    rec = _Recorder()
    t.register = lambda data, key, artifact_type=None, **kw: rec.register(data, key, t, artifact_type, **kw)
    t._recorder = rec
    return t


PREBUILT = {'psf': '/scratch/bilayer.psf', 'pdb': '/scratch/bilayer.pdb', 'xsc': '/scratch/bilayer.xsc'}

#: A plausible equilibrated bilayer cell; only the diagonal is used for the area.
BOX = [[92.0, 0.0, 0.0], [0.0, 88.0, 0.0], [0.0, 0.0, 110.0]]
ORIGIN = [0.0, 0.0, 0.0]


class TestPrebuiltBilayerProvision(unittest.TestCase):

    def _provision(self, prebuilt=PREBUILT):
        t = _task(prebuilt)
        with mock.patch('pestifer.tasks.make_membrane_system.BaseTask.provision'), \
             mock.patch('pestifer.tasks.make_membrane_system._cell_or_raise',
                        return_value=(BOX, ORIGIN)) as cell, \
             mock.patch('pestifer.tasks.make_membrane_system.get_toppar_from_psf',
                        return_value=['toppar_all36_lipid_cholesterol.str']):
            MakeMembraneSystemTask.provision(t, {})
        return t, cell

    def test_a_prebuilt_bilayer_provisions_without_raising(self):
        """The reported failure: AttributeError: 'dict' object has no attribute 'xsc'."""
        t, _ = self._provision()
        self.assertTrue(t.using_prebuilt_bilayer)

    def test_the_cell_is_read_from_the_registered_artifact_not_the_request_dict(self):
        """The distinction the bug turned on.

        `register` is handed a dict of paths and returns an artifact; the cell must be read off
        the artifact.  Asserting only that provision survived would pass if someone reverted to
        the dict and the dict happened to grow an `xsc` attribute, so pin the actual argument.
        """
        t, cell = self._provision()
        cell.assert_called_once()
        passed = cell.call_args[0][0]
        self.assertEqual(passed, PREBUILT['xsc'])
        self.assertIsNot(passed, PREBUILT)

    def test_the_box_and_area_come_from_that_cell(self):
        t, _ = self._provision()
        self.assertEqual(t.quilt.box, BOX)
        self.assertEqual(t.quilt.area, 92.0 * 88.0)

    def test_topologies_are_read_from_the_registered_psf(self):
        t, _ = self._provision()
        self.assertEqual(t.quilt.addl_streamfiles, ['toppar_all36_lipid_cholesterol.str'])

    def test_the_state_is_registered_under_quilt_state(self):
        """Downstream reads it back by that key (`get_current_artifact('quilt_state')`), so the
        key is part of the contract, not an internal detail."""
        t, _ = self._provision()
        self.assertIn('quilt_state', t._recorder.registered)
        self.assertEqual(t._recorder.registered['quilt_state'].data, PREBUILT)

    def test_a_prebuilt_block_without_a_pdb_does_not_take_this_branch(self):
        """The branch is gated on `pdb` being present, so a partial block must fall through to
        normal construction rather than half-initializing a quilt."""
        t = _task({'psf': '/scratch/bilayer.psf'})
        with mock.patch('pestifer.tasks.make_membrane_system.BaseTask.provision'):
            MakeMembraneSystemTask.provision(t, {})
        self.assertFalse(t.using_prebuilt_bilayer)
        t.initialize.assert_called_once()


class TestIgnoredRelaxationProtocolsAreAnnounced(unittest.TestCase):
    """A `quilt:` protocol beside `prebuilt:` never runs, and used to say nothing about it.

    Skipping relaxation is correct -- the bilayer is supplied because it is already equilibrated
    -- but a config that carries a protocol reasonably looks like the protocol will be used, and
    on a restart workflow the difference is a membrane that did or did not get further NPgT time.
    """

    def _provision_with(self, protocols):
        t = _task(PREBUILT)
        t.specs['bilayer']['relaxation_protocols'] = protocols
        with mock.patch('pestifer.tasks.make_membrane_system.BaseTask.provision'), \
             mock.patch('pestifer.tasks.make_membrane_system._cell_or_raise',
                        return_value=(BOX, ORIGIN)), \
             mock.patch('pestifer.tasks.make_membrane_system.get_toppar_from_psf',
                        return_value=[]), \
             mock.patch('pestifer.tasks.make_membrane_system.logger') as log:
            MakeMembraneSystemTask.provision(t, {})
        return [str(c.args[0]) for c in log.warning.call_args_list]

    def test_an_ignored_quilt_protocol_is_warned_about_by_name(self):
        warnings = self._provision_with({'quilt': [{'md': {'nsteps': 1000}}]})
        self.assertEqual(len(warnings), 1)
        self.assertIn('quilt', warnings[0])
        self.assertIn('will NOT run', warnings[0])

    def test_both_protocols_are_named_when_both_are_present(self):
        warnings = self._provision_with({'patch': [{'md': {}}], 'quilt': [{'md': {}}]})
        self.assertIn('patch,quilt', warnings[0])

    def test_no_warning_when_no_protocol_was_requested(self):
        """Guards against warning on every prebuilt build, which would train users to ignore it."""
        self.assertEqual(self._provision_with({}), [])

    def test_no_warning_for_an_empty_protocol_block(self):
        self.assertEqual(self._provision_with({'quilt': []}), [])
