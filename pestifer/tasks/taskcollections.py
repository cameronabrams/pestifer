import logging

logger = logging.getLogger(__name__)

from .basetask import BaseTask
from collections import UserList, UserDict
from . import task_classes

class TaskList(UserList[BaseTask]):
    """ 
    A list of BaseTask objects.
    """
    
    @classmethod
    def from_yaml(cls, task_list: list[dict]):
        """
        Create a TaskList from a YAML list of task specifications.
        """
        data = []
        prior: BaseTask = None
        index: int = 0
        for task_unidict in task_list:
            logger.debug(f'Processing task specification: {task_unidict}')
            assert len(task_unidict) == 1, f"Task dictionary {task_unidict} must have a single key-value pair"
            taskname = list(task_unidict.keys())[0]
            task_specs = task_unidict[taskname]
            if not 'taskname' in task_specs:
                task_specs['taskname'] = taskname
            if not 'index' in task_specs:
                task_specs['index'] = index
            # Ensure the name of the task is among the implemented Tasks
            Cls = task_classes.get(taskname, None)
            if Cls is None:
                raise ValueError(f"Task {taskname} is not implemented.")
            data.append(Cls(specs=task_specs))
            if prior is None:
                prior = data[-1]
            else:
                data[-1].prior = prior
                prior = data[-1]
            index += 1
        return cls(data)

class TaskDict(UserDict[str, dict]):
    """
    A dictionary of task specifications.
    """
    pass