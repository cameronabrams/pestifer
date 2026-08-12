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
