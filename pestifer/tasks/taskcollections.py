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
            assert len(task_unidict) == 1, f"Task dictionary {task_unidict} must have a single key-value pair"
            taskname = list(task_unidict.keys())[0]
            task_specs = task_unidict[taskname]
            task_specs['taskname'] = taskname
            # Ensure the name of the task is among the implemented Tasks
            class_name = task_classes.get(taskname, None)
            if not class_name:
                raise ValueError(f"Task {taskname} is not implemented.")
            Cls = task_classes[class_name]
            data.append(Cls(specs=task_specs))
            data[-1].index = index
            if prior is None:
                prior = data[-1]
            else:
                data[-1].prior = prior
                prior = data[-1]
        return cls(data)

class TaskDict(UserDict[str, dict]):
    """
    A dictionary of task specifications.
    """
    pass