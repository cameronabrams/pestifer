import unittest

from pestifer.tasks.taskcollections import TaskList
from pestifer.tasks import TerminateTask
from pestifer.tasks.pipeline_contract import (
    validate_pipeline, TaskContract, SOURCE, STATE, MOLECULE, MD_OUTPUT, SOLVATED)
from pestifer.core.errors import PestiferBuildError


def _build(task_specs, add_terminate=True):
    """Build a real task list from a compact spec (list of dicts) the way the controller
    does, assign indices, and optionally append the default terminate."""
    tl = TaskList.from_yaml(task_specs)
    if add_terminate and (not tl or not isinstance(tl[-1], TerminateTask)):
        tl.append(TerminateTask(specs={}, index=len(tl)))
    for i, t in enumerate(tl):
        t.index = i
    return tl


class TestPipelineContract(unittest.TestCase):

    def _assert_ok(self, specs):
        validate_pipeline(_build(specs))  # raises on failure

    def _assert_error(self, specs, needle):
        with self.assertRaises(PestiferBuildError) as cm:
            validate_pipeline(_build(specs))
        self.assertIn(needle, str(cm.exception))

    # ---- valid pipelines ----
    def test_fetch_psfgen_md_ok(self):
        self._assert_ok([{'fetch': {'sourceID': '6pti'}}, {'psfgen': {}}, {'md': {}}])

    def test_continuation_transform_ok(self):
        self._assert_ok([{'continuation': {'psf': 'x.psf', 'pdb': 'x.pdb'}},
                         {'ring_check': {'segtypes': ['lipid']}}, {'md': {}}])

    def test_membrane_embed_and_mdplot_ok(self):
        self._assert_ok([{'fetch': {'sourceID': '6e8w'}}, {'psfgen': {}},
                         {'make_membrane_system': {'embed': {'xydist': 10}}},
                         {'md': {}}, {'mdplot': {}}])

    # ---- missing prerequisite (category A) ----
    def test_continuation_then_psfgen_is_the_motivating_error(self):
        self._assert_error([{'continuation': {'psf': 'x.psf', 'pdb': 'x.pdb'}}, {'psfgen': {}}],
                           "no `fetch` task precedes it")

    def test_transform_first_needs_state(self):
        self._assert_error([{'md': {}}], "requires an existing molecular system")

    def test_base_molecule_needed_after_continuation(self):
        # ligate needs the in-memory molecule, which continuation does not provide
        self._assert_error([{'continuation': {'psf': 'x.psf', 'pdb': 'x.pdb'}}, {'ligate': {}}],
                           "needs the in-memory molecule")

    def test_mdplot_without_md(self):
        self._assert_error([{'fetch': {'sourceID': '6pti'}}, {'psfgen': {}}, {'mdplot': {}}],
                           "needs molecular-dynamics output")

    def test_dangling_fetch(self):
        self._assert_error([{'fetch': {'sourceID': '6pti'}}], "requires an existing molecular system")

    # ---- terminal rules (category D) ----
    def test_nothing_after_terminate(self):
        self._assert_error([{'fetch': {'sourceID': '6pti'}}, {'psfgen': {}},
                            {'terminate': {}}, {'md': {}}],
                           "scheduled after the terminal task")

    # ---- discarding origin (category B) ----
    def test_second_psfgen_discards_state(self):
        self._assert_error([{'fetch': {'sourceID': '6pti'}}, {'psfgen': {}}, {'psfgen': {}}],
                           "discard")

    # ---- contract details ----
    def test_membrane_contract_is_conditional_on_embed(self):
        from pestifer.tasks.make_membrane_system import MakeMembraneSystemTask as MM
        self.assertIn(STATE, MM.pipeline_contract({'embed': {'xydist': 10}}).requires)
        self.assertEqual(frozenset(), MM.pipeline_contract({}).requires)

    def test_default_contract_is_a_transform(self):
        c = TaskContract()
        self.assertEqual(c.requires, frozenset({STATE}))
        self.assertEqual(c.provides, frozenset({STATE}))


if __name__ == '__main__':
    unittest.main()
