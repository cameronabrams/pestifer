import logging
from collections import UserList, UserDict


from . import task_classes
from .basetask import BaseTask
from ..util.stringthings import my_logger

logger = logging.getLogger(__name__)

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
        for idx, task_unidict in enumerate(task_list):
            logger.debug(f'Processing specification for task {idx:02d}:')
            my_logger(task_unidict, logger.debug)
            assert len(task_unidict) == 1, f"Task dictionary must have a single key-value pair"
            taskname = list(task_unidict.keys())[0]
            task_specs = task_unidict[taskname]
            # Ensure the name of the task is among the implemented Tasks
            Cls = task_classes.get(taskname, None)
            if Cls is None:
                raise ValueError(f"Task {taskname} is not implemented.")
            data.append(Cls(specs=task_specs, index=index))
            if prior is None:
                prior = data[-1]
            else:
                data[-1].prior = prior
                prior = data[-1]
            index += 1
        return cls(data)

class TaskDict(UserDict[str, BaseTask]):
    """
    A dictionary of task specifications.
    """
    pass