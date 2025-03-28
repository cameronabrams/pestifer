import pytest
from unittest.mock import MagicMock, patch
from pestifer.basetask import BaseTask

@pytest.fixture
def mock_task():
    # Mock dependencies
    mock_writers = {'vmd': MagicMock()}
    mock_prior = MagicMock()
    input_dict = {'index': 0}
    config = {'user': {'namd': {'harmonic': {'spring_constant': 10}}}}
    return BaseTask(input_dict, "test_task", config, mock_writers, mock_prior)

def test_initialization(mock_task):
    assert mock_task.index == 0
    assert mock_task.taskname == "test_task"
    assert mock_task.statevars == {}
    assert mock_task.subtaskcount == 0

def test_do_method(mock_task):
    assert mock_task.do() == 0

def test_str_method(mock_task):
    assert str(mock_task) == "0/0 - test_task [has prior True]"

def test_log_message(mock_task):
    with patch("pestifer.basetask.logger.info") as mock_logger:
        mock_task.log_message("started", extra_info="test")
        mock_logger.assert_called_once_with("Task 00 'test_task' started  (extra_info: test)")

def test_next_basename(mock_task):
    mock_task.next_basename("subtask")
    assert mock_task.basename == "00-00-test_task-subtask"
    assert mock_task.subtaskcount == 1

def test_update_statevars_strict(mock_task):
    with pytest.raises(FileNotFoundError):
        mock_task.update_statevars("key", "nonexistent.file", vtype="file", mode="strict")

def test_update_statevars_permissive(mock_task):
    with patch("os.path.exists", return_value=False):
        mock_task.update_statevars("key", "nonexistent.file", vtype="file", mode="permissive")
        assert "key" not in mock_task.statevars

def test_inherit_state(mock_task):
    mock_task.prior.statevars = {"key": "value"}
    mock_task.inherit_state()
    assert mock_task.statevars == {"key": "value"}

def test_save_state(mock_task):
    with patch("pestifer.basetask.BaseTask.update_statevars") as mock_update:
        mock_task.next_basename("subtask")
        mock_task.save_state(["pdb", "coor"])
        assert mock_update.call_count == 2

def test_write_statefile(mock_task):
    with patch("builtins.open", MagicMock()) as mock_open, patch("yaml.dump") as mock_yaml:
        mock_task.next_basename("subtask")
        mock_task.write_statefile()
        mock_open.assert_called_once_with("00-00-test_task-subtask-state.yaml", "w")
        mock_yaml.assert_called_once_with(mock_task.statevars, mock_open.return_value.__enter__())
