import unittest
import pytest
from pestifer.molecule.stateinterval import StateInterval, StateIntervalList

class TestStateInterval(unittest.TestCase):
    def setUp(self):
        self.interval = StateInterval(state='A', bounds=[1, 10], build=False)

    def test_describe(self):
        description = self.interval.describe()
        expected_description = "StateInterval(state=A, bounds=[1, 10], build=False)"
        self.assertEqual(description, expected_description)

    def test_declare_buildable(self):
        self.interval.declare_buildable()
        self.assertTrue(self.interval.build)

    def test_increment_rightbound(self):
        self.interval.increment_rightbound(5)
        self.assertEqual(self.interval.bounds, [1, 15])

    def test_increment_rightbound_invalid(self):
        with self.assertRaises(ValueError):
            self.interval.bounds = ['a', 'b']
            self.interval.increment_rightbound()

class TestStateIntervalList(unittest.TestCase):

    def test_add_interval_internal(self):
        self.interval_list = StateIntervalList()
        self.interval1 = StateInterval(state='A', bounds=[1, 5], build=False)
        self.interval2 = StateInterval(state='B', bounds=[6, 10], build=False)
        self.interval_list.append(self.interval1)
        self.interval_list.append(self.interval2)
        nested_interval = StateInterval(state='A', bounds=[2,3], build=False)
        self.interval_list.add_interval(nested_interval)
        # this should do nothing since the states are the same and the entire interval is already in the list
        self.assertEqual(len(self.interval_list), 2)
        nested_interval = StateInterval(state='B', bounds=[2,3], build=False)
        self.interval_list.add_interval(nested_interval)
        # this should split interval1, and yield a list of four intervals
        self.assertEqual(len(self.interval_list), 4)
        self.assertEqual(self.interval_list[0].bounds, [1, 1])
        self.assertEqual(self.interval_list[0].state, 'A')
        self.assertEqual(self.interval_list[1].bounds, [2, 3])
        self.assertEqual(self.interval_list[1].state, 'B')
        self.assertEqual(self.interval_list[2].bounds, [4, 5])
        self.assertEqual(self.interval_list[2].state, 'A')
        self.assertEqual(self.interval_list[3].bounds, [6, 10])
        self.assertEqual(self.interval_list[3].state, 'B')

    def test_add_interval_across_left(self):
        self.interval_list = StateIntervalList()
        self.interval1 = StateInterval(state='A', bounds=[1, 5], build=False)
        self.interval2 = StateInterval(state='B', bounds=[6, 10], build=False)
        self.interval_list.append(self.interval1)
        self.interval_list.append(self.interval2)
        nested_interval = StateInterval(state='A', bounds=[4,7], build=False)
        self.interval_list.add_interval(nested_interval)
        # this should grow interval1 and shrink interval2, with no change in the number of intervals
        self.assertEqual(len(self.interval_list), 2)
        self.assertEqual(self.interval_list[0].bounds, [1, 7])
        self.assertEqual(self.interval_list[0].state, 'A')
        self.assertEqual(self.interval_list[1].bounds, [8, 10])

    def test_add_interval_across_right(self):
        self.interval_list = StateIntervalList()
        self.interval1 = StateInterval(state='A', bounds=[1, 5], build=False)
        self.interval2 = StateInterval(state='B', bounds=[6, 10], build=False)
        self.interval_list.append(self.interval1)
        self.interval_list.append(self.interval2)
        nested_interval = StateInterval(state='B', bounds=[4,7], build=False)
        self.interval_list.add_interval(nested_interval)
        # this should grow interval2 and shrink interval1, with no change in the number of intervals
        self.assertEqual(len(self.interval_list), 2)
        self.assertEqual(self.interval_list[0].bounds, [1, 3])
        self.assertEqual(self.interval_list[0].state, 'A')
        self.assertEqual(self.interval_list[1].bounds, [4, 10])
        self.assertEqual(self.interval_list[1].state, 'B')

    def test_add_interval_across_inserts(self):
        self.interval_list = StateIntervalList()
        self.interval1 = StateInterval(state='A', bounds=[1, 5], build=False)
        self.interval2 = StateInterval(state='B', bounds=[6, 10], build=False)
        self.interval_list.append(self.interval1)
        self.interval_list.append(self.interval2)
        nested_interval = StateInterval(state='C', bounds=[3,8], build=False)
        self.interval_list.add_interval(nested_interval)
        # this should insert a new interval between the two existing intervals
        self.assertEqual(len(self.interval_list), 3)
        self.assertEqual(self.interval_list[0].bounds, [1, 2])
        self.assertEqual(self.interval_list[0].state, 'A')
        self.assertEqual(self.interval_list[1].bounds, [3, 8])
        self.assertEqual(self.interval_list[1].state, 'C')
        self.assertEqual(self.interval_list[2].bounds, [9, 10])
        self.assertEqual(self.interval_list[2].state, 'B')

    def test_add_interval_no_suitable_member(self):
        self.interval_list = StateIntervalList()
        self.interval1 = StateInterval(state='A', bounds=[1, 5], build=False)
        self.interval2 = StateInterval(state='B', bounds=[6, 10], build=False)
        self.interval_list.append(self.interval1)
        self.interval_list.append(self.interval2)
        nested_interval = StateInterval(state='C', bounds=[11, 15], build=False)
        with self.assertRaises(ValueError):
            self.interval_list.add_interval(nested_interval)