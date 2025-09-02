"""
Tasks that can be executed by the Pestifer framework.
This package contains various task classes that define the operations to be performed during the execution of a Pestifer run.
Each task is responsible for a specific part of the workflow, such as building molecules, running simulations, and processing results.

Each module defines a specific task class that inherits from the base :class:`BaseTask <pestifer.core.basetask.BaseTask>` class.
"""
from .cleave import CleaveTask
from .continuation import ContinuationTask
from .desolvate import DesolvateTask
from .domainswap import DomainSwapTask
from .fetch import FetchTask
from .ligate import LigateTask
from .make_membrane_system import MakeMembraneSystemTask
from .manipulate import ManipulateTask
from .mdtask import MDTask
from .mdplot import MDPlotTask
from .pdb2pqr import PDB2PQRTask
from .psfgen import PsfgenTask
from .ringcheck import RingCheckTask
from .solvate import SolvateTask
from .terminate import TerminateTask
from .validate import ValidateTask

task_classes: dict[str, type] = {
    'cleave': CleaveTask,
    'continuation': ContinuationTask,
    'desolvate': DesolvateTask,
    'domainswap': DomainSwapTask,
    'fetch': FetchTask,
    'ligate': LigateTask,
    'make_membrane_system': MakeMembraneSystemTask,
    'md': MDTask,
    'manipulate': ManipulateTask,
    'mdplot': MDPlotTask,
    'pdb2pqr': PDB2PQRTask,
    'psfgen': PsfgenTask,
    'ring_check': RingCheckTask,
    'solvate': SolvateTask,
    'terminate': TerminateTask,
    'validate': ValidateTask,
}