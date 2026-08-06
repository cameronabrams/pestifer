# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Task numbering for branch pipelines.

A branch pipeline opens with a ``continuation`` that loads an intermediate state file from another
run (e.g. ``00-02-00_md.coor``).  Its task numbering should continue the parent's lineage: the
continuation adopts the state's encoded task index, and the next task is +1 -- so a branch's output
filenames read as one continuous lineage across runs, rather than restarting at 00.
"""
import unittest

from pestifer.tasks.taskcollections import TaskList
from pestifer.tasks.basetask import parse_basename_task_index


class TestParseBasenameTaskIndex(unittest.TestCase):

    def test_index_encoded_names(self):
        self.assertEqual(parse_basename_task_index('00-02-00_md.coor'), 2)
        self.assertEqual(parse_basename_task_index('01-13-02_solvate.pdb'), 13)
        self.assertEqual(parse_basename_task_index('/a/b/00-05-00_md.pdb'), 5)   # basename only

    def test_non_encoded_names(self):
        for n in ('my_structure.pdb', 'protein.psf', 'foo.coor', '', None):
            self.assertIsNone(parse_basename_task_index(n))


class TestBranchNumbering(unittest.TestCase):

    def _indices(self, tl):
        return [t.index for t in tl]

    def test_fresh_pipeline_starts_at_zero(self):
        tl = TaskList.from_yaml([{'fetch': {'sourceID': '6pti'}}, {'psfgen': {}}, {'md': {}}])
        self.assertEqual(self._indices(tl), [0, 1, 2])

    def test_branch_adopts_encoded_index(self):
        # continuation from the state at task 02 -> continuation IS 02, next task is 03
        tl = TaskList.from_yaml([
            {'continuation': {'psf': '00-01-00_psfgen.psf', 'pdb': '00-02-00_md.pdb'}},
            {'md': {}}, {'solvate': {}}])
        self.assertEqual(self._indices(tl), [2, 3, 4])
        self.assertEqual(tl[0].taskname, 'continuation')

    def test_index_taken_from_coordinate_not_psf(self):
        # psf carries an earlier index (MD does not rewrite it); the coordinate file is authoritative
        tl = TaskList.from_yaml([
            {'continuation': {'psf': '00-01-00_psfgen.psf', 'pdb': '00-06-00_md.pdb'}}, {'md': {}}])
        self.assertEqual(self._indices(tl), [6, 7])

    def test_coor_used_when_no_pdb(self):
        tl = TaskList.from_yaml([
            {'continuation': {'psf': 'x.psf', 'coor': '00-07-00_md.coor'}}, {'md': {}}])
        self.assertEqual(self._indices(tl), [7, 8])

    def test_pdb_preferred_over_coor(self):
        tl = TaskList.from_yaml([
            {'continuation': {'psf': 'x.psf', 'pdb': '00-03-00_md.pdb', 'coor': '00-09-00_md.coor'}},
            {'md': {}}])
        self.assertEqual(self._indices(tl), [3, 4])

    def test_non_encoded_continuation_starts_at_zero(self):
        tl = TaskList.from_yaml([
            {'continuation': {'psf': 'x.psf', 'pdb': 'mystructure.pdb'}}, {'md': {}}])
        self.assertEqual(self._indices(tl), [0, 1])

    def test_continuation_not_first_no_offset(self):
        # offset only applies when the pipeline OPENS with a continuation
        tl = TaskList.from_yaml([{'fetch': {'sourceID': '6pti'}}, {'psfgen': {}}])
        self.assertEqual(self._indices(tl), [0, 1])

    def test_branch_basenames_are_continuous(self):
        tl = TaskList.from_yaml([
            {'continuation': {'psf': 'x.psf', 'pdb': '00-02-00_md.pdb'}},
            {'md': {}}, {'solvate': {}}])
        # controller_index is assigned at provision time; simulate the top-level controller (0)
        for t in tl:
            t.controller_index = 0
            t.next_basename()
        self.assertEqual([t.basename for t in tl],
                         ['00-02-000_continuation', '00-03-000_md', '00-04-000_solvate'])
