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
        # the minimize between `ring_check` and the dynamics is required (see
        # TestUnminimizedDynamics below); this test is about the STATE hand-off, and the
        # hand-off is unchanged by it
        self._assert_ok([{'continuation': {'psf': 'x.psf', 'pdb': 'x.pdb'}},
                         {'ring_check': {'segtypes': ['lipid']}},
                         {'md': {'ensemble': 'minimize'}}, {'md': {}}])

    def test_membrane_embed_and_mdplot_ok(self):
        self._assert_ok([{'fetch': {'sourceID': '6e8w'}}, {'psfgen': {}},
                         {'make_membrane_system': {'embed': {'z_head_group': 'name CA',
                                                             'z_tail_group': 'name CB'}}},
                         {'md': {}}, {'mdplot': {}}])

    def test_bare_membrane_is_an_origin(self):
        # a protein-free bilayer build needs no preceding state -- it originates one
        self._assert_ok([{'make_membrane_system': {'bilayer': {'packer': 'grid'}}}, {'md': {}}])

    def test_continuation_then_psfgen_is_preserve_mode_ok(self):
        # A psfgen after a STATE-provider (continuation) with no fetched SOURCE is now VALID:
        # psfgen edits the incoming topology in place (readpsf-preserve) rather than rebuilding.
        self._assert_ok([{'continuation': {'psf': 'x.psf', 'pdb': 'x.pdb'}}, {'psfgen': {}}, {'md': {}}])

    # ---- missing prerequisite (category A) ----
    def test_bare_psfgen_needs_source(self):
        # psfgen with neither a fetched SOURCE nor a prior STATE still errors up front.
        self._assert_error([{'psfgen': {}}], "no `fetch` task precedes it")

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

    def test_desolvate_rejected_in_build_pipeline(self):
        # desolvate is a standalone utility; it must not appear in a build tasks list
        self._assert_error([{'continuation': {'psf': 'x.psf', 'pdb': 'x.pdb'}}, {'desolvate': {}}],
                           "standalone utility")

    # ---- contract details ----
    def test_membrane_contract_is_conditional_on_embed(self):
        from pestifer.tasks.make_membrane_system import MakeMembraneSystemTask as MM
        # a genuine embed (protein-positioning directive) requires a protein state to insert
        self.assertIn(STATE, MM.pipeline_contract({'embed': {'z_head_group': 'name CA'}}).requires)
        self.assertIn(STATE, MM.pipeline_contract({'embed': {'no_orient': True}}).requires)
        # no embed at all -> origin
        self.assertEqual(frozenset(), MM.pipeline_contract({}).requires)
        # the schema's default `embed` scaffolding (no positioning directive) is NOT an embed:
        # this is the real-pipeline case that must not spuriously require a protein state
        schema_default_embed = {'no_orient': False, 'q_tolerance': 1.e-4, 'xydist': 15,
                                'zdist': 20, 'margin': 2.4,
                                'orient': {'source': [], 'target': []}, 'z_ref_group': {}}
        self.assertEqual(frozenset(), MM.pipeline_contract({'embed': schema_default_embed}).requires)
        # and it must never be flagged as discarding an incoming protein state
        self.assertFalse(MM.pipeline_contract({'embed': schema_default_embed}).discards_state)

    def test_default_contract_is_a_transform(self):
        c = TaskContract()
        self.assertEqual(c.requires, frozenset({STATE}))
        self.assertEqual(c.provides, frozenset({STATE}))


if __name__ == '__main__':
    unittest.main()


class TestUnminimizedDynamics(unittest.TestCase):
    """Dynamics must not integrate coordinates an earlier task moved and nothing minimized.

    `make_membrane_system._guard_piercings` enforces this for the relaxation protocol it
    generates, and says why: a piercing resolution "frees the other lipid into a strained pose",
    and a threading "RATTLE-fails the first dynamics step".  That guard never saw the user's own
    top-level task list, so example 17 shipped `minimize -> ring_check -> 410 K dynamics` -- the
    same pair with the minimize on the wrong side -- for several releases.
    """

    def _tasks(self, specs):
        return _build(specs)

    def test_ring_check_then_dynamics_is_rejected(self):
        with self.assertRaises(PestiferBuildError) as cm:
            validate_pipeline(self._tasks([
                {'continuation': {'psf': 'x.psf', 'pdb': 'x.pdb'}},
                {'ring_check': {'segtypes': ['lipid']}},
                {'md': {'ensemble': 'NVT'}}]))
        msg = str(cm.exception)
        self.assertIn('ring_check', msg)
        self.assertIn('no minimization in between', msg)
        self.assertIn('ensemble: minimize', msg)   # the message must name the remedy

    def test_a_minimize_in_between_clears_it(self):
        validate_pipeline(self._tasks([
            {'continuation': {'psf': 'x.psf', 'pdb': 'x.pdb'}},
            {'ring_check': {'segtypes': ['lipid']}},
            {'md': {'ensemble': 'minimize'}},
            {'md': {'ensemble': 'NVT'}}]))

    def test_a_minimize_BEFORE_the_mover_does_not_count(self):
        """Example 17's actual defect: it minimized, then ring_checked, then ran dynamics."""
        with self.assertRaises(PestiferBuildError):
            validate_pipeline(self._tasks([
                {'continuation': {'psf': 'x.psf', 'pdb': 'x.pdb'}},
                {'md': {'ensemble': 'minimize'}},
                {'ring_check': {'segtypes': ['lipid']}},
                {'md': {'ensemble': 'NVT'}}]))

    def test_membrane_equilibrate_counts_as_dynamics(self):
        with self.assertRaises(PestiferBuildError):
            validate_pipeline(self._tasks([
                {'continuation': {'psf': 'x.psf', 'pdb': 'x.pdb'}},
                {'ring_check': {'segtypes': ['lipid']}},
                {'membrane_equilibrate': {}}]))

    def test_each_mover_is_reported_once_not_once_per_following_stage(self):
        with self.assertRaises(PestiferBuildError) as cm:
            validate_pipeline(self._tasks([
                {'continuation': {'psf': 'x.psf', 'pdb': 'x.pdb'}},
                {'ring_check': {'segtypes': ['lipid']}},
                {'md': {'ensemble': 'NVT'}},
                {'md': {'ensemble': 'NPT'}}]))
        self.assertEqual(str(cm.exception).count('no minimization in between'), 1)

    def test_every_shipped_example_satisfies_the_rule(self):
        """The rule must not be retro-fitted to a set of pipelines that violate it."""
        from importlib.resources import files as pkg_files
        from pathlib import Path
        import yaml
        examples = Path(str(pkg_files('pestifer.resources'))) / 'examples'
        files = sorted(examples.glob('*/inputs/*.yaml'))
        self.assertGreater(len(files), 25, f'no example configs under {examples}')
        for f in files:
            with self.subTest(example=f.parent.parent.name):
                specs = (yaml.safe_load(f.read_text()) or {}).get('tasks') or []
                validate_pipeline(_build(specs, add_terminate=False))
