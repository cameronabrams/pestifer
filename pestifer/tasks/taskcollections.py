import logging
from collections import UserList, UserDict


from . import task_classes
from .basetask import BaseTask, parse_basename_task_index
from ..core.errors import PestiferBuildError
from ..util.stringthings import my_logger

logger = logging.getLogger(__name__)


def _branch_index_offset(task_list: list[dict]) -> int:
    """Task-index offset for a *branch* pipeline that opens with a ``continuation``.

    When the first task is a ``continuation`` loading an index-encoded state file from another run
    (e.g. ``00-02-00_md.coor``), the branch adopts that task's number so its numbering continues
    from the parent: the continuation takes the state's index and the next task is +1.  The offset is
    that encoded index (so ``index = position + offset``).  Returns 0 for an ordinary from-scratch
    pipeline or a continuation whose state file is not index-encoded (unchanged behavior).

    The index is read from the **coordinate** file (``pdb`` preferred, else ``coor``) -- that is the
    task whose state we are resuming (``psf`` may carry an earlier index since MD does not rewrite it).
    """
    if not task_list:
        return 0
    first_name = list(task_list[0].keys())[0]
    if first_name != 'continuation':
        return 0
    specs = task_list[0][first_name] or {}
    coord = specs.get('pdb') or specs.get('coor')
    encoded = parse_basename_task_index(coord)
    return encoded if encoded is not None else 0


class TaskList(UserList[BaseTask]):
    """
    A list of BaseTask objects.
    """

    @classmethod
    def from_yaml(cls, task_list: list[dict]):
        """
        Create a TaskList from a YAML list of task specifications.

        A branch pipeline (one opening with a ``continuation`` from another run's index-encoded state)
        adopts that state's task number so the branch's output filenames continue the parent's
        lineage numbering (see :func:`_branch_index_offset`); ordinary pipelines start at 0.
        """
        data = []
        prior: BaseTask = None
        offset = _branch_index_offset(task_list)
        for idx, task_unidict in enumerate(task_list):
            logger.debug(f'Processing specification for task {idx:02d}:')
            my_logger(task_unidict, logger.debug)
            assert len(task_unidict) == 1, f"Task dictionary must have a single key-value pair"
            taskname = list(task_unidict.keys())[0]
            task_specs = task_unidict[taskname]
            # Ensure the name of the task is among the implemented Tasks
            Cls = task_classes.get(taskname, None)
            if Cls is None:
                raise PestiferBuildError(f"Task {taskname} is not implemented.")
            data.append(Cls(specs=task_specs, index=idx + offset))
            if prior is None:
                prior = data[-1]
            else:
                data[-1].prior = prior
                prior = data[-1]
        return cls(data)

class TaskDict(UserDict[str, BaseTask]):
    """
    A dictionary of task specifications.
    """
    pass