# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A mutation is emitted as a psfgen ``mutate`` command; it never changes the molecule's own
residue list.  So the topology collection that walks that list sees the residue being mutated
AWAY and never the one being mutated TO -- and psfgen stops with ``unknown residue type`` for
any target whose defining stream is not loaded for some other reason.
"""
import unittest
from types import SimpleNamespace
from unittest import mock

from pestifer.tasks.psfgen import PsfgenTask


class _Mut:
    def __init__(self, newresname):
        self.newresname = newresname


def _task(mutations, topfile_of):
    t = PsfgenTask.__new__(PsfgenTask)
    om = {'seq': {'mutations': [_Mut(m) for m in mutations]}}
    t.base_molecule = SimpleNamespace(objmanager=om)
    t.resource_manager = SimpleNamespace(
        charmmff_content=SimpleNamespace(get_topfile_of_resname=lambda r: topfile_of.get(r)))
    return t


class TestMutationTopologies(unittest.TestCase):

    def test_a_mutation_target_pulls_in_its_defining_stream(self):
        t = _task(['SEP'], {'SEP': 'toppar_all36_prot_na_combined.str'})
        self.assertEqual(t.mutation_topologies(), ['toppar_all36_prot_na_combined.str'])

    def test_all_three_phosphorylations_share_one_stream_and_it_is_added_once(self):
        m = {'SEP': 'toppar_all36_prot_na_combined.str',
             'TPO': 'toppar_all36_prot_na_combined.str',
             'PTR': 'toppar_all36_prot_na_combined.str'}
        t = _task(['SEP', 'TPO', 'PTR'], m)
        self.assertEqual(t.mutation_topologies(), ['toppar_all36_prot_na_combined.str'])

    def test_several_targets_in_different_streams(self):
        t = _task(['SEP', 'TYS'], {'SEP': 'a.str', 'TYS': 'b.str'})
        self.assertEqual(sorted(t.mutation_topologies()), ['a.str', 'b.str'])

    def test_no_mutations_asks_for_nothing(self):
        self.assertEqual(_task([], {}).mutation_topologies(), [])

    def test_an_undefined_target_is_left_for_psfgen_to_report(self):
        # mutating to a residue the force field does not define is the user's error; psfgen
        # names it.  Raising here would report it as a topology-resolution failure instead.
        t = _task(['NOPE'], {})
        self.assertEqual(t.mutation_topologies(), [])

    def test_a_standard_target_needs_nothing_extra_beyond_its_own_file(self):
        t = _task(['ALA'], {'ALA': 'top_all36_prot.rtf'})
        self.assertEqual(t.mutation_topologies(), ['top_all36_prot.rtf'])


class TestTheCollectorIsActuallyWired(unittest.TestCase):
    """The tests above call ``mutation_topologies`` directly, so all of them pass if the method
    is written and never called.  This pins the call site."""

    def test_psfgen_asks_for_mutation_topologies(self):
        import inspect
        from pestifer.tasks import psfgen as mod
        src = inspect.getsource(mod.PsfgenTask.psfgen)
        self.assertIn('mutation_topologies()', src,
                      'psfgen() must collect the topologies its mutations need')
        # and they must reach psfgen, not merely be computed
        line = [l for l in src.splitlines() if 'required_topology_files' in l and '=' in l]
        self.assertTrue(line, 'expected a required_topology_files assignment')
        joined = ' '.join(src.splitlines())
        self.assertIn('additional_topologies=required_topology_files', joined)
