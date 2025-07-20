import unittest
from unittest.mock import MagicMock, patch
from pestifer.core.basetask import BaseTask
from pestifer.core.pipeline import PipelineContext

class TestBaseTask(unittest.TestCase):
    def setUp(self):
        # Mock dependencies
        self.ctx = PipelineContext()
        self.mock_writers = {'vmd': MagicMock()}
        self.mock_prior = MagicMock()
        self.mock_prior.index = 0
        config_specs = {'example_config': 99}
        config = {'user': {'namd': {'harmonic': {'spring_constant': 10}}}}
        controller_specs = {
            'config': config,
            'writers': self.mock_writers,
            'prior': self.mock_prior,
            'controller_index': 0,
            'taskname': 'test_task'
        }
        self.task = BaseTask(self.ctx, config_specs, controller_specs)
    
    def test_basetask_initialization(self):
        assert self.task.specs['example_config'] == 99
        assert self.task.index == 1
        assert self.task.controller_index == 0
        assert self.task.taskname == "test_task"
        assert self.task.ctx == self.ctx
        assert self.task.subtaskcount == 0

    def test_basetask_do_method(self):
        assert self.task.do() == 0

    def test_basetask_str_method(self):
        assert str(self.task) == "0 - 1 - test_task [has prior True]"

    def test_basetask_log_message(self):
        with patch("pestifer.core.basetask.logger.info") as mock_logger:
            self.task.log_message("started", extra_info="test")
            mock_logger.assert_called_once_with("Controller 00 Task 01 'test_task' started  (extra_info: test)")

    def test_basetask_next_basename(self):
        self.task.next_basename("subtask")
        assert self.task.basename == "00-01-00_test_task-subtask"
        assert self.task.subtaskcount == 1

    def test_basetask_override_taskname(self):
        self.task.override_taskname("new_taskname")
        assert self.task.taskname == "new_taskname"