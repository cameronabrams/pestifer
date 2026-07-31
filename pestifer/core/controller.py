#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
A controller for the pestifer runtime.  Initialization of a Controller object generates
the configuration, from which the list of tasks is created.  The :meth:`do_tasks()` method executes
the tasks.
"""
import logging
import os

from dataclasses import dataclass

from pestifer.tasks.validate import ValidateTask

from .config import Config
from .errors import PestiferBuildError
from .pipeline import PipelineContext
from ..tasks.taskcollections import TaskList
from ..tasks import TerminateTask
from ..tasks.pipeline_contract import validate_pipeline
from ..util.stringthings import plu, my_logger
from ..util.util import running_under_pytest

logger = logging.getLogger(__name__)

@dataclass
class Controller:
    """ 
    Controller class for managing the execution of tasks in the Pestifer runtime.
    This class initializes with a configuration and a list of user-specified tasks.
    It sets up the necessary file writers and creates a list of tasks to be executed.
    The tasks are executed in the order they are defined, and the results of each task
    are collected in a report. If the last task is not a :class:`pestifer.tasks.terminate.TerminateTask` task,
    a default :class:`pestifer.tasks.terminate.TerminateTask` is added to ensure proper termination.

    A Controller is constructed with no arguments and then set up with :meth:`configure`, whose
    parameters are:

    config : Config
        The configuration object containing user specifications and settings.
    userspecs : dict, optional
        A dictionary of user-specific configurations provided at run-time that will update the base configuration.
    index : int, optional
        An index for the controller instance, useful for identifying multiple controllers in a run.
    terminate : bool, optional
        If True, a default terminate task will be added if the last task is not already a TerminateTask.
        This ensures that the runtime will always have a way to clean up and terminate properly.
        Default is True.
    validate : bool, optional
        If True (default), the assembled task list is statically checked for malformed hand-offs
        (see :func:`pestifer.tasks.pipeline_contract.validate_pipeline`) before execution.  Standalone
        utility subcommands set this False.
    """

    def __init__(self):
        self.index = 0
        self.config = None
        self.tasks = None
        self.pipeline = None
        self.parent = None
        self.restart = False   # --restart: resume from the last cleanly-completed task
        self.fresh = False     # --fresh: ignore any existing manifest
        self.from_task = None  # --from: resume explicitly at this task (index or taskname)
        self.packet = None     # provisioning packet, kept for the resume state-restore

    def configure(self, config: Config, userspecs: dict = {}, index: int = 0, terminate: bool = True,
                  validate: bool = True):
        self.index = index
        self.config = config
        self.config.update_user(userspecs)

        # set up the task list
        self.tasks = TaskList.from_yaml(self.config['user'].get('tasks', []))
        user_supplied_tasks = len(self.tasks) > 0
        if terminate and (len(self.tasks) == 0 or not isinstance(self.tasks[-1], TerminateTask)):
            # If the last task is not a TerminateTask, add one with default specs
            specs = self.config.make_default_specs('tasks','terminate')
            logger.debug('Adding default terminate task')
            self.tasks.append(TerminateTask(specs=specs, index=len(self.tasks)))

        # static pre-execution check of the task hand-offs (only the user-facing top-level
        # pipeline has tasks at configure() time; subcontrollers are populated later and are
        # generated internally, so they are not checked here).  Standalone utility subcommands
        # (mdplot, desolvate) that assemble a single-task controller operating on CLI inputs
        # pass validate=False -- they are not build pipelines.  A controller configured with no
        # user tasks (a placeholder later populated via reconfigure_tasks, or a no-op build)
        # holds only the auto-added terminate and is not a pipeline to validate.
        if validate and user_supplied_tasks:
            validate_pipeline(self.tasks)

        self.pipeline = PipelineContext(controller_index=self.index)
        self.provision_tasks()

        logger.info(f'Running "{self.config["user"]["title"]}"')
        logger.info(f'Controller {self.index:02d} will execute {len(self.tasks)} task{plu(len(self.tasks))}.')
        return self

    @classmethod
    def spawn_subcontroller(cls, progenitor: 'Controller') -> 'Controller':
        logger.debug(f'Spawning subcontroller from {progenitor.index:02d}')
        subcontroller = cls()
        subcontroller.configure(progenitor.config.taskless_subconfig(), index=progenitor.index+1, terminate=False)
        subcontroller.parent = progenitor
        return subcontroller
    
    def reconfigure_tasks(self, tasks: list[dict]):
        new_data = dict(tasks=tasks)
        self.config.update_user(new_data)
        self.tasks = TaskList.from_yaml(self.config['user'].get('tasks', []))
        self.provision_tasks()
        logger.info(f'Reconfigured tasks in "{self.config["user"]["title"]}"')
        logger.info(f'Controller {self.index:02d} will execute {len(self.tasks)} task{plu(len(self.tasks))}.')

    def provision_tasks(self):
        packet = {
            'controller_index': self.index,
            'pipeline': self.pipeline,
            'resource_manager': self.config.RM,
            'scripters': self.config.scripters,
            'processor-type': self.config.namd_type,
            'shell-commands': self.config.shell_commands,
            'progress-flag': self.config.use_terminal_progress,
            'namd_global_config': self.config['user']['namd'],
        }
        self.packet = packet   # kept so a --restart can provision a state-restore continuation
        logger.debug(f'Controller {self.index:02d} provisioning {len(self.tasks)} task{plu(len(self.tasks))} with packet:')
        my_logger(packet, logger.debug)
        for task in self.tasks:
            logger.debug('*'*70)
            logger.debug(f'Task {self.index:02d}:{task.index:02d} ** {task.taskname:>20s} **')
            task.provision(packet)
            if task.specs.get('requires_subcontroller', False):
                logger.debug('*'*70)
                logger.debug(f'Controller {self.index:02d} spawning subcontroller for Task {task.index:02d} ** {task.taskname:>20s} **')
                task.subcontroller = Controller.spawn_subcontroller(self)

    def do_tasks(self) -> dict:
        """
        Execute the tasks in the order they were defined.

        This method iterates through the list of tasks, executing each one in turn.
        It collects the results of each task in a report dictionary, which maps task indices to their
        names, indices, and results. If any task fails (i.e., returns a non-zero result),
        a warning is logged, and the execution of subsequent tasks is aborted.

        Returns
        -------
        dict
            A dictionary containing the results of each task, with task indices as keys.
            Each value is itself a dictionary with the following keys:

            - ``taskname``: The name of the task
            - ``taskindex``: The index of the task
            - ``result``: The task's return code
        """
        task_report = {}
        task_durations = 0
        manifest, resume_from = self._init_run_manifest()   # (None, 0) for subcontrollers
        if resume_from > 0:
            self._restore_state_for_resume(manifest, resume_from)
            self._clean_resumed_task_outputs(resume_from)
        for task in self.tasks:
            if task.index < resume_from:
                logger.info(f"--restart: skipping completed task {task.index:02d} '{task.taskname}'")
                continue
            pre_hash = None
            if manifest is not None:
                from .run_manifest import spec_hash
                pre_hash = spec_hash(task.specs)   # BEFORE execute() -- some tasks mutate their specs
            returned_result = task.execute()
            task_report[task.index] = dict(taskname=task.taskname, taskindex=task.index, result=returned_result)
            task_durations += task.duration
            if task.result != 0:
                logger.warning(f'Task {task.taskname} failed; task.result {task.result} returned result {returned_result} controller is aborted.')
                break
            if manifest is not None:
                manifest.record(task, self.pipeline, task_spec_hash=pre_hash)   # durable resume point
        if manifest is not None and task_report and all(r['result'] == 0 for r in task_report.values()):
            manifest.mark_complete()   # a finished build has nothing to resume
        for task in self.tasks:
            if task.index not in task_report:
                continue
            task_report[task.index]['duration'] = task.duration
            task_report[task.index]['duration_frac'] = task.duration/task_durations if task_durations > 0 else 0
        return task_report

    def _init_run_manifest(self):
        """Return ``(manifest, resume_from)`` for a top-level build, or ``(None, 0)``.

        The top-level controller (index 0) writes ``.pestifer-manifest.json``; subcontroller
        pipelines are opaque and appear as a single entry under their parent task.  With
        ``--restart`` (and not ``--fresh``), an existing manifest is loaded, the resume point
        computed, and recording continues into it; otherwise a fresh manifest is started.  Any
        failure is non-fatal (falls back to a from-scratch build).
        """
        if self.index != 0:
            return None, 0
        from .run_manifest import RunManifest, tasks_fingerprint, MANIFEST_NAME
        path = os.path.join(os.getcwd(), MANIFEST_NAME)
        version = ''
        try:
            from importlib.metadata import version as _pkg_version
            version = _pkg_version('pestifer')
        except Exception:
            pass
        try:
            fingerprint = tasks_fingerprint(self.tasks)
        except Exception as exc:
            logger.warning(f'run manifest: could not initialize ({exc}); resume will be unavailable')
            return None, 0

        from_idx = self._resolve_from_task()          # explicit --from index (or None)
        want_resume = (self.restart or from_idx is not None) and not self.fresh
        if want_resume and os.path.exists(path):
            try:
                manifest = RunManifest.load(path)
                if from_idx is not None:
                    resume_from = from_idx
                    logger.info(f'--from: resuming explicitly at task {resume_from:02d}')
                else:
                    rp = manifest.resume_point(self.tasks)
                    if rp < 0:
                        logger.info('--restart: nothing to resume (no matching manifest or build '
                                    'already complete); building from scratch')
                        resume_from = 0
                    else:
                        resume_from = rp + 1
                        logger.info(f'--restart: tasks 00-{rp:02d} already complete; resuming at '
                                    f'task {resume_from:02d}')
                manifest.data.update(pestifer_version=version, tasks_fingerprint=fingerprint,
                                     complete=False)
                return manifest, resume_from
            except Exception as exc:
                logger.warning(f'--restart: could not read run manifest ({exc}); building from scratch')
        elif from_idx is not None and not os.path.exists(path):
            raise PestiferBuildError(
                f'--from {self.from_task!r}: no run manifest ({MANIFEST_NAME}) in this directory to '
                f'restore state from; run the build here first')

        try:
            return RunManifest(path, version=version, fingerprint=fingerprint), 0
        except Exception as exc:
            logger.warning(f'run manifest: could not initialize ({exc}); resume will be unavailable')
            return None, 0

    def _resolve_from_task(self):
        """Resolve ``--from`` (a task index or taskname) to a task index in this pipeline, or None."""
        if not self.from_task:
            return None
        ft = str(self.from_task)
        if ft.isdigit():
            return int(ft)
        for task in self.tasks:
            if task.taskname == ft:
                return task.index
        raise PestiferBuildError(
            f'--from {self.from_task!r}: no task with that index or name in the pipeline')

    def _clean_resumed_task_outputs(self, resume_from):
        """Delete on-disk output files of every task that will re-run (index >= resume_from), so a
        half-written partial left by the interrupt can't be mistaken for a completed output.  The
        restored STATE (from task ``resume_from - 1``) is a lower index, so it is never touched."""
        import glob
        removed = 0
        for task in self.tasks:
            if task.index < resume_from:
                continue
            for f in glob.glob(f'{self.index:02d}-{task.index:02d}-*'):
                if os.path.isfile(f):
                    try:
                        os.remove(f)
                        removed += 1
                    except OSError:
                        pass
        if removed:
            logger.info(f'--restart: cleared {removed} stale/partial output file(s) from re-run tasks')

    def _restore_state_for_resume(self, manifest, resume_from):
        """Re-establish the pipeline STATE (psf/pdb/coor/xsc/vel + CHARMM streams) at the resume
        boundary from the manifest, by running a `continuation` against the last completed task's
        recorded fileset -- the same primitive a manual branch uses.  Raises on failure (a resume
        with no restored state would fail the first tail task confusingly)."""
        k = resume_from - 1
        state = manifest.state_entry(k)
        if not state:
            logger.info(f'--restart: task {k:02d} recorded no STATE fileset; the resume task is an '
                        f'origin, nothing to restore')
            return
        from ..tasks.continuation import ContinuationTask
        logger.info(f'--restart: restoring STATE from task {k:02d}: {state}')
        ct = ContinuationTask(specs=dict(state), index=k)
        if self.packet is not None:
            ct.provision(self.packet)
        ct.controller_index = self.index
        if ct.execute() != 0:
            raise PestiferBuildError(
                f'--restart: could not restore state from task {k:02d} ({state}); '
                f'the recorded files may be missing or corrupt -- rebuild with --fresh')

    def write_complete_config(self, filename='complete-user.yaml'):
        """ 
        Write the complete user configuration to a YAML file.
        This method dumps the user configuration, including all modifications made during the execution of tasks,
        to a specified YAML file. The default filename is 'complete-user.yaml'.
        
        Parameters
        ----------
        filename : str, optional
            The name of the file to which the user configuration will be written. Default is 'complete-user.yaml'.
        """
        self.config.dump_user(filename=filename)