# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the ``task-table`` Sphinx directive.

The directive reads a pestifer YAML config and renders the pipeline-summary table that heads each
worked example in the documentation.  It is pure data-to-prose: a per-task summary function, a
grouper that collapses consecutive ``md`` tasks, and a small RST emitter.  None of it needs sphinx
to be running, which is why almost all of it is tested here as plain functions.

Two things make it worth testing at all.  It runs at *docs build* time, so a failure takes the
whole build down rather than one page -- and it degrades silently: an unrecognized task renders as
an em dash, which looks like a deliberate blank rather than a gap.  The last test in this file is
the guard for that: it runs the grouper over every shipped example and asserts that every task
type actually used has a summary.

Note these tests only run where sphinx is installed; the directory's conftest skips them
otherwise.  Sphinx is in the ``test`` extra for exactly that reason.
"""

import glob
import os
import tempfile
import unittest
from types import SimpleNamespace

import yaml

from pestifer.sphinxext.task_table import (
    _equilibrate_details, _fetch_details, _manipulate_details, _md_group_details,
    _membrane_details, _psfgen_details, _single_task_details, _table_rst_lines, _task_rows,
)


class TestFetchDetail(unittest.TestCase):

    def test_an_rcsb_entry_becomes_a_link(self):
        got = _fetch_details({'sourceID': '6pti'})
        self.assertIn('PDB 6PTI', got, 'the identifier should be upper-cased for display')
        self.assertIn('https://www.rcsb.org/structure/6pti', got)

    def test_an_alphafold_model_is_named_as_one(self):
        got = _fetch_details({'source': 'alphafold', 'sourceID': 'P22033'})
        self.assertIn('AlphaFold', got)
        self.assertIn('P22033', got)
        self.assertNotIn('rcsb.org', got, 'an AlphaFold model has no PDB entry to link to')

    def test_empty_specs_do_not_raise(self):
        self.assertIn('?', _fetch_details({}))
        self.assertIn('?', _fetch_details(None))


class TestPsfgenDetail(unittest.TestCase):

    def test_a_plain_build_says_so(self):
        self.assertEqual(_psfgen_details({}), 'standard build')
        self.assertEqual(_psfgen_details(None), 'standard build')

    def test_each_kind_of_modification_is_named(self):
        got = _psfgen_details({'mods': {'mutations': ['A:1A'], 'ssbonds': ['x'],
                                        'ssbondsdelete': ['y'], 'substitutions': ['z']}})
        for expected in ('mutations', 'new disulfides', 'delete disulfides', 'loop substitutions'):
            self.assertIn(expected, got)

    def test_grafts_are_counted(self):
        self.assertIn('3 glycan graft(s)', _psfgen_details({'mods': {'grafts': [1, 2, 3]}}))

    def test_exclusions_live_under_source_not_mods(self):
        self.assertIn('exclusions', _psfgen_details({'source': {'exclude': {'chains': ['B']}}}))

    def test_a_source_given_as_a_list_uses_its_first_entry(self):
        self.assertIn('exclusions',
                      _psfgen_details({'source': [{'exclude': {'chains': ['B']}}]}))

    def test_an_empty_source_list_does_not_raise(self):
        self.assertEqual(_psfgen_details({'source': []}), 'standard build')

    def test_loop_modeling_is_recognized_from_either_place(self):
        self.assertIn('loop modeling', _psfgen_details({'mods': {'loops': ['a']}}))
        self.assertIn('loop modeling', _psfgen_details({'source': {'loops': ['a']}}))


class TestMdGrouping(unittest.TestCase):
    """Consecutive md tasks collapse into one row; this is what keeps a 12-phase warm-up from
    filling the table."""

    def test_a_minimization_is_named_not_counted(self):
        self.assertEqual(_md_group_details([{'ensemble': 'minimize'}]), 'minimize')

    def test_a_single_phase_reports_its_steps(self):
        self.assertEqual(_md_group_details([{'ensemble': 'NVT', 'nsteps': 2000}]),
                         'NVT (2,000 steps)')

    def test_consecutive_phases_of_one_ensemble_are_summed_and_counted(self):
        got = _md_group_details([{'ensemble': 'NPT', 'nsteps': 1000},
                                 {'ensemble': 'NPT', 'nsteps': 3000}])
        self.assertEqual(got, 'NPT (4,000 steps, 2 phases)')

    def test_a_change_of_ensemble_starts_a_new_segment(self):
        got = _md_group_details([{'ensemble': 'minimize'},
                                 {'ensemble': 'NVT', 'nsteps': 1000},
                                 {'ensemble': 'NPT', 'nsteps': 2000}])
        self.assertEqual(got, 'minimize → NVT (1,000 steps) → NPT (2,000 steps)')

    def test_ensembles_are_matched_case_insensitively(self):
        got = _md_group_details([{'ensemble': 'NPT', 'nsteps': 10},
                                 {'ensemble': 'npt', 'nsteps': 20}])
        self.assertIn('2 phases', got)

    def test_defaults_apply_when_a_phase_omits_them(self):
        self.assertEqual(_md_group_details([{}]), 'NVT (2,000 steps)')


class TestEquilibrateDetail(unittest.TestCase):
    """These tasks choose their own duration, so the summary must describe the criterion rather
    than a step count -- and must not present the ceiling as the plan."""

    def test_a_density_equilibration_names_its_observable(self):
        got = _equilibrate_details('density_equilibrate', {})
        self.assertIn('self-terminating', got)
        self.assertIn('density', got)

    def test_a_membrane_equilibration_names_both_observables(self):
        got = _equilibrate_details('membrane_equilibrate', {})
        self.assertIn('density', got)
        self.assertIn('lateral area', got)

    def test_the_two_stage_protocol_is_described_and_is_the_default(self):
        self.assertIn('const-area settle', _equilibrate_details('membrane_equilibrate', {}))
        self.assertNotIn('const-area settle',
                         _equilibrate_details('membrane_equilibrate', {'two_stage': False}))

    def test_max_steps_is_presented_as_a_ceiling(self):
        got = _equilibrate_details('density_equilibrate', {'max_steps': 200000})
        self.assertIn('ceiling 200,000 steps', got)


class TestSingleTaskDetails(unittest.TestCase):

    def test_water_is_the_default_solvent(self):
        self.assertEqual(_single_task_details('solvate', {}), 'water box')

    def test_non_water_uses_the_solvent_name(self):
        self.assertEqual(_single_task_details('solvate', {'solvent': 'CHCL3'}), 'CHCL3 box')

    def test_salt_and_ions_are_reported_together(self):
        got = _single_task_details('solvate',
                                   {'salt_con': 0.154, 'cation': 'SOD', 'anion': 'CLA'})
        self.assertIn('0.154 M SOD/CLA', got)

    def test_salt_without_named_ions_still_reports_the_concentration(self):
        self.assertIn('0.154 M salt', _single_task_details('solvate', {'salt_con': 0.154}))

    def test_ligate_reports_its_steered_md_length(self):
        self.assertIn('8,000 steps',
                      _single_task_details('ligate', {'steer': {'nsteps': 8000}}))
        self.assertEqual(_single_task_details('ligate', {}), 'ligate chain breaks')

    def test_cleave_counts_its_sites(self):
        self.assertIn('2 site(s)', _single_task_details('cleave', {'sites': ['a', 'b']}))

    def test_pdb2pqr_reports_the_ph(self):
        self.assertIn('pH 7.4', _single_task_details('pdb2pqr', {'pH': 7.4}))
        self.assertIn('pdb2pqr', _single_task_details('pdb2pqr', {}))

    def test_validate_counts_its_tests(self):
        self.assertEqual(_single_task_details('validate', {'tests': [1, 2, 3]}), '3 test(s)')

    def test_mdplot_names_its_output_directory(self):
        self.assertIn('mdplots/', _single_task_details('mdplot', {}))
        self.assertIn('figs/', _single_task_details('mdplot', {'output_dir': 'figs'}))

    def test_terminate_reports_its_basenames(self):
        got = _single_task_details('terminate',
                                  {'basename': 'my_sys', 'package': {'basename': 'prod'}})
        self.assertIn('my_sys', got)
        self.assertIn('prod', got)
        self.assertEqual(_single_task_details('terminate', {}), '—')

    def test_manipulate_lists_its_mods(self):
        self.assertIn('``translate``',
                      _manipulate_details({'mods': {'translate': {}}}))
        self.assertEqual(_manipulate_details({}), '—')

    def test_manipulate_is_reachable_through_the_dispatcher(self):
        self.assertIn('``translate``',
                      _single_task_details('manipulate', {'mods': {'translate': {}}}))

    def test_continuation_names_the_psf_it_resumes_from(self):
        self.assertIn('my.psf', _single_task_details('continuation', {'psf': 'my.psf'}))
        self.assertEqual(_single_task_details('continuation', {}), '—')

    def test_the_remaining_fixed_summaries(self):
        self.assertIn('ring-threading', _single_task_details('ring_check', {}))
        self.assertIn('merge', _single_task_details('merge', {}))

    def test_a_membrane_system_lists_its_lipids_deduplicated_and_sorted(self):
        got = _membrane_details({'bilayer': {'composition': {
            'upper_leaflet': [{'name': 'POPC'}, {'name': 'CHL1'}],
            'lower_leaflet': [{'name': 'POPC'}, {'name': 'PSM'}]}}})
        self.assertEqual(got, 'CHL1, POPC, PSM')

    def test_a_membrane_system_without_composition_does_not_raise(self):
        self.assertEqual(_membrane_details({}), '—')

    def test_an_unknown_task_falls_back_to_an_em_dash(self):
        self.assertEqual(_single_task_details('no_such_task', {}), '—')


class TestTaskRows(unittest.TestCase):

    def test_consecutive_md_tasks_become_one_row(self):
        tasks = [{'fetch': {'sourceID': '6pti'}},
                 {'md': {'ensemble': 'minimize'}},
                 {'md': {'ensemble': 'NVT', 'nsteps': 1000}},
                 {'terminate': {'basename': 'x'}}]
        rows = _task_rows(tasks)
        self.assertEqual([name for name, _ in rows], ['fetch', 'md', 'terminate'])

    def test_md_groups_separated_by_another_task_stay_separate(self):
        tasks = [{'md': {'ensemble': 'NVT', 'nsteps': 10}},
                 {'solvate': {}},
                 {'md': {'ensemble': 'NPT', 'nsteps': 20}}]
        self.assertEqual([n for n, _ in _task_rows(tasks)], ['md', 'solvate', 'md'])

    def test_a_task_with_null_specs_does_not_raise(self):
        self.assertEqual([n for n, _ in _task_rows([{'md': None}])], ['md'])


class TestRstEmitter(unittest.TestCase):

    TASKS = [{'fetch': {'sourceID': '6pti'}}, {'solvate': {}}]

    def test_the_table_is_a_well_formed_list_table(self):
        lines = _table_rst_lines(self.TASKS)
        self.assertTrue(lines[0].startswith('.. list-table::'))
        self.assertIn('   :header-rows: 1', lines)
        self.assertIn('   * - Step', lines)

    def test_steps_are_numbered_from_one(self):
        lines = _table_rst_lines(self.TASKS)
        self.assertIn('   * - 1', lines)
        self.assertIn('   * - 2', lines)

    def test_task_names_are_rendered_as_literals(self):
        self.assertIn('     - ``fetch``', _table_rst_lines(self.TASKS))

    def test_an_empty_pipeline_still_yields_a_valid_header(self):
        lines = _table_rst_lines([])
        self.assertTrue(lines[0].startswith('.. list-table::'))


class TestDirective(unittest.TestCase):
    """``TaskTableDirective.run`` needs only four things from sphinx -- the source directory, the
    document name, a reporter and a parser -- so it is driven here with stand-ins rather than a
    full app.  The path resolution and the missing-file branch are the parts worth pinning: the
    latter is what a docs author sees when a relative path is wrong."""

    def _directive(self, tmpdir, arg):
        from pestifer.sphinxext.task_table import TaskTableDirective
        d = TaskTableDirective.__new__(TaskTableDirective)
        d.arguments = [arg]
        d.lineno = 1
        d.content_offset = 0
        # SphinxDirective.env is a read-only property reading through state.document.settings,
        # so the stand-in has to be threaded the same way rather than assigned
        env = SimpleNamespace(srcdir=tmpdir, docname='page')
        d.state = SimpleNamespace(
            document=SimpleNamespace(settings=SimpleNamespace(env=env)),
            nested_parse=self._record_parse)
        d.reporter = SimpleNamespace(error=lambda msg, line=None: ('error', msg))
        self.parsed = None
        return d

    def _record_parse(self, viewlist, offset, node):
        self.parsed = list(viewlist)

    def test_a_config_is_resolved_relative_to_the_rst_file_and_parsed(self):
        with tempfile.TemporaryDirectory() as tmp:
            with open(os.path.join(tmp, 'ex.yaml'), 'w') as f:
                yaml.dump({'tasks': [{'fetch': {'sourceID': '6pti'}}]}, f)
            d = self._directive(tmp, 'ex.yaml')
            d.run()
        self.assertTrue(self.parsed[0].startswith('.. list-table::'))
        self.assertTrue(any('6PTI' in line for line in self.parsed))

    def test_a_missing_config_reports_an_error_naming_the_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = self._directive(tmp, 'nope.yaml')
            result = d.run()
        self.assertEqual(result[0][0], 'error')
        self.assertIn('nope.yaml', result[0][1])
        self.assertIn('not found', result[0][1])

    def test_a_config_without_a_tasks_key_yields_an_empty_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            with open(os.path.join(tmp, 'ex.yaml'), 'w') as f:
                yaml.dump({'title': 'no tasks here'}, f)
            d = self._directive(tmp, 'ex.yaml')
            d.run()
        self.assertTrue(self.parsed[0].startswith('.. list-table::'))


class TestEveryShippedExampleSummarizes(unittest.TestCase):
    """The guard that matters.

    The directive degrades silently: a task it does not recognize renders as an em dash, which in
    a rendered table is indistinguishable from a deliberate blank.  Adding a task type to pestifer
    without adding it here therefore produces documentation that is quietly wrong rather than
    broken.  This walks every shipped example and fails if any task summarizes to nothing.
    """

    @staticmethod
    def _example_configs():
        here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        root = os.path.join(here, '..', 'pestifer', 'resources', 'examples')
        return sorted(glob.glob(os.path.join(root, '*', 'inputs', '*.yaml')))

    def test_the_examples_are_actually_found(self):
        # a path typo here would make the whole guard vacuously pass
        self.assertGreater(len(self._example_configs()), 20)

    def test_every_task_in_every_example_has_a_summary(self):
        blanks = []
        for path in self._example_configs():
            with open(path) as f:
                config = yaml.safe_load(f)
            for name, details in _task_rows((config or {}).get('tasks', []) or []):
                if details == '—':
                    blanks.append(f'{os.path.basename(path)}: {name}')
        self.assertEqual(blanks, [], f'task types with no summary: {blanks}')

    def test_every_example_renders_without_raising(self):
        for path in self._example_configs():
            with open(path) as f:
                config = yaml.safe_load(f)
            lines = _table_rst_lines((config or {}).get('tasks', []) or [])
            self.assertTrue(lines[0].startswith('.. list-table::'), path)


if __name__ == '__main__':
    unittest.main()
