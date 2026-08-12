# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Which recorded MD runs belong in a plot.

`mdplot` used to collect data by walking back through preceding tasks while they were MDTask
instances, stopping at the first that was not.  That conflates "not an MD task" with "invalidated
the earlier data": it is right when the interrupting task is `solvate` (a different system) and
wrong when it is `validate` or `ring_check` (the same system).  Every bundled example was
truncated by it, ex16 and ex17 spuriously.  Selection is now made from the recorded atom count.
"""
import unittest

from pestifer.tasks.mdplot import select_commensurable_runs


def _run(task_index, natoms, name='md'):
    return dict(controller=0, task_index=task_index, taskname=name, basename=f'{task_index:02d}',
                ensemble='npt', natoms=natoms, csv_keys=['energy'], csv_basename=f'{task_index:02d}')


class TestSelectCommensurableRuns(unittest.TestCase):

    def test_solvate_boundary_ends_the_span(self):
        # ex01 shape: minimize on the dry system, solvate, then everything on the solvated one
        entries = [_run(3, 882), _run(5, 14130), _run(6, 14130), _run(7, 14130)]
        keep = select_commensurable_runs(entries, 14130)
        self.assertEqual([e['task_index'] for e in keep], [5, 6, 7])

    def test_non_system_changing_task_does_not_truncate(self):
        # ex17 shape: MD, then ring_check/validate (no atom-count change), then more MD.
        # The old adjacency walk stopped at the intervening task and lost the earlier stage.
        entries = [_run(4, 250000), _run(7, 250000), _run(8, 250000)]
        keep = select_commensurable_runs(entries, 250000)
        self.assertEqual([e['task_index'] for e in keep], [4, 7, 8])

    def test_chunked_task_contributes_many_runs_to_one_stage(self):
        entries = [_run(5, 100)] + [_run(7, 100) for _ in range(20)]
        keep = select_commensurable_runs(entries, 100)
        self.assertEqual(len(keep), 21)
        self.assertEqual(len({(e['controller'], e['task_index']) for e in keep}), 2)

    def test_only_the_trailing_span_is_kept(self):
        # two solvation events: only the runs after the most recent one belong
        entries = [_run(2, 500), _run(4, 9000), _run(6, 20000), _run(7, 20000)]
        keep = select_commensurable_runs(entries, 20000)
        self.assertEqual([e['task_index'] for e in keep], [6, 7])

    def test_unknown_counts_are_not_treated_as_a_boundary(self):
        # an unrecorded count is not evidence of a different system; truncating there would
        # reintroduce the silent-loss failure this replaces
        entries = [_run(3, None), _run(4, 14130)]
        keep = select_commensurable_runs(entries, 14130)
        self.assertEqual(len(keep), 2)

    def test_unknown_current_system_keeps_everything(self):
        entries = [_run(3, 882), _run(5, 14130)]
        self.assertEqual(len(select_commensurable_runs(entries, None)), 2)

    def test_empty_history(self):
        self.assertEqual(select_commensurable_runs([], 100), [])


if __name__ == '__main__':
    unittest.main()


class TestAutoBlockAverage(unittest.TestCase):
    """Automatic smoothing targets only observables whose fluctuations dwarf their drift."""

    def test_instantaneous_pressures_are_smoothed(self):
        from pestifer.tasks.mdplot import auto_block_average
        self.assertGreater(auto_block_average('PRESSURE', 800), 0)
        self.assertGreater(auto_block_average('GPRESSURE', 800), 0)

    def test_well_behaved_observables_are_left_alone(self):
        from pestifer.tasks.mdplot import auto_block_average
        # DENSITY and TEMP read fine raw; PRESSAVG/GPRESSAVG are already NAMD running averages,
        # so smoothing them would average an average
        for col in ('DENSITY', 'TEMP', 'PRESSAVG', 'GPRESSAVG', 'TOTAL', 'VOLUME'):
            with self.subTest(column=col):
                self.assertEqual(auto_block_average(col, 800), 0)

    def test_window_scales_with_series_length(self):
        from pestifer.tasks.mdplot import auto_block_average
        self.assertLess(auto_block_average('PRESSURE', 400), auto_block_average('PRESSURE', 4000))

    def test_short_series_is_not_smoothed_away(self):
        from pestifer.tasks.mdplot import auto_block_average
        # a window comparable to the series would flatten it into a single value
        self.assertEqual(auto_block_average('PRESSURE', 8), 0)


class TestOverlayArgParsing(unittest.TestCase):
    """LABEL=PATTERN overlay arguments."""

    def setUp(self):
        import tempfile, os
        self.dir = tempfile.mkdtemp()
        for n in ('b.log', 'a.log', 'c.log'):
            open(os.path.join(self.dir, n), 'w').close()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_glob_expands_sorted(self):
        # chunked runs must concatenate chronologically, not in filesystem order
        from pestifer.subcommands.mdplot import _parse_overlays
        out = _parse_overlays([f'run one={self.dir}/*.log'])
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]['label'], 'run one')
        self.assertEqual([p.split('/')[-1] for p in out[0]['logs']], ['a.log', 'b.log', 'c.log'])

    def test_comma_separated_list_is_kept_in_order(self):
        from pestifer.subcommands.mdplot import _parse_overlays
        out = _parse_overlays([f'x={self.dir}/c.log,{self.dir}/a.log'])
        self.assertEqual([p.split('/')[-1] for p in out[0]['logs']], ['c.log', 'a.log'])

    def test_malformed_and_empty_are_skipped_not_fatal(self):
        from pestifer.subcommands.mdplot import _parse_overlays
        self.assertEqual(_parse_overlays(['no-equals-sign']), [])
        self.assertEqual(_parse_overlays([f'empty={self.dir}/nothing-*.log']), [])
