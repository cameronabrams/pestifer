# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Implementation of the ``build`` subcommand for launching system preparations.
"""

import datetime
import json
import logging
import os
import shutil
import sys
import time
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

from ..core.controller import Controller
from ..core.config import Config
from ..core.errors import PestiferError
from ..cli.subcommand import Subcommand
from ..util.stringthings import __pestifer_version__
from ..util.provenance import log_environment
from ..util.citations import log_citations
from ..core.run_record import build_run_record, write_run_record
from ..util.util import hmsf

#: Extensions that only ever name an *input* to a build.  A config value ending in one of
#: these is a file the run will need to read, so ``--check`` verifies it exists rather than
#: letting the build discover it missing several tasks in.  Output names in pestifer configs
#: are basenames without extensions (``basename: my_system``), so this does not flag them.
_INPUT_EXTENSIONS = ('.pdb', '.cif', '.mmcif', '.psf', '.coor', '.xsc', '.vel',
                     '.str', '.rtf', '.prm', '.mol2', '.sdf', '.dcd')

def _walk_strings(node):
    """Yield every string value in a nested dict/list structure."""
    if isinstance(node, str):
        yield node
    elif isinstance(node, dict):
        for v in node.values():
            yield from _walk_strings(v)
    elif isinstance(node, (list, tuple)):
        for v in node:
            yield from _walk_strings(v)


def _missing_input_files(config):
    """Config values that name an input file which is not on disk.

    Heuristic and advisory: a value is checked only when it ends in an extension that
    pestifer never uses for an output basename, so a hit is a file the build will try to
    read.  Reported as a warning rather than an error because a task may legitimately
    resolve a name against a search path rather than the working directory.
    """
    missing = []
    for s in _walk_strings(config['user'].get('tasks', [])):
        if s.lower().endswith(_INPUT_EXTENSIONS) and not os.path.exists(s):
            missing.append(s)
    return sorted(set(missing))


def _preflight(configname, args, **kwargs):
    """Assemble the config and controller without running anything.

    Everything that makes a build fail in its first seconds -- an unparsable config, a key
    the schema does not define, a missing ``vmd``/``namd3``, an impossible task hand-off --
    is already detected by ``Config.configure_new`` and ``Controller.configure``.  This runs
    exactly those, collects what they found, and stops short of ``do_tasks``.

    Returns the report dict; ``report['ok']`` is False if the build would not start.
    """
    report = {
        'ok': False,
        'config': configname,
        'pestifer_version': __pestifer_version__,
        'title': None,
        'executables': {},
        'namd_type': None,
        'ncpus': None,
        'charmmff': None,
        'tasks': [],
        'warnings': [],
        'errors': [],
    }

    try:
        config = Config(userfile=configname, ncpus_override=args.ncpus,
                        processor_type_override=('gpu' if getattr(args, 'gpu', False) else ''),
                        seed_override=getattr(args, 'seed', None),
                        **kwargs).configure_new()
    except PestiferError as e:
        report['errors'].append(str(e))
        return report
    except Exception as e:
        report['errors'].append(f'{type(e).__name__}: {e}')
        return report

    report['title'] = config['user'].get('title')
    report['executables'] = {name: shutil.which(cmd) for name, cmd in config.shell_commands.items()}
    report['namd_type'] = config.namd_type
    report['ncpus'] = config.ncpus
    try:
        report['charmmff'] = config.RM.charmmff_version_path(
            config['user']['charmmff'].get('release', '')).name
    except Exception as e:                      # a bad release pin should not mask the rest
        report['errors'].append(f'CHARMMFF: {e}')

    try:
        C = Controller().configure(config)
    except PestiferError as e:
        report['errors'].append(str(e))
        return report
    except Exception as e:
        report['errors'].append(f'{type(e).__name__}: {e}')
        return report

    report['tasks'] = [{'index': i, 'taskname': t.taskname} for i, t in enumerate(C.tasks)]

    for f in _missing_input_files(config):
        report['warnings'].append(f'input file not found in the working directory: {f}')

    manifest = os.path.join(os.getcwd(), '.pestifer-manifest.json')
    if os.path.exists(manifest):
        report['warnings'].append(
            'a run manifest from a previous build is present; use --restart to resume it, '
            'or --fresh to ignore it and start over')

    # The diagnostics log is opened by the CLI before this runs, so it is always present and
    # is not evidence that the directory has been used for anything.
    ignore = {os.path.basename(configname)}
    clutter = [e for e in os.listdir('.')
               if not e.startswith('.') and e not in ignore
               and not e.startswith('pestifer_diagnostics.log')
               and not e.endswith(('-diagnostics.log', '-diagnostics.log.bak'))]
    if clutter:
        report['warnings'].append(
            f'working directory is not empty ({len(clutter)} other entries); a build writes '
            f'many files and expects a directory of its own')

    report['ok'] = not report['errors']
    return report


def _print_check_report(report):
    """Human-readable form of a ``--check`` report, on stdout."""
    w = sys.stdout.write
    w(f'pestifer {report["pestifer_version"]} build --check: {report["config"]}\n\n')
    if report['title']:
        w(f'  title       {report["title"]}\n')
    if report['charmmff']:
        w(f'  charmmff    {report["charmmff"]}\n')
    if report['namd_type']:
        w(f'  namd        {report["namd_type"]} mode, {report["ncpus"]} PEs\n')
    for name, path in report['executables'].items():
        w(f'  {name:<11} {path or "NOT FOUND"}\n')

    if report['tasks']:
        w(f'\n  task plan ({len(report["tasks"])} tasks):\n')
        for t in report['tasks']:
            w(f'    {t["index"]:>2}  {t["taskname"]}\n')

    for warning in report['warnings']:
        w(f'\nWARNING: {warning}\n')
    for error in report['errors']:
        w(f'\nERROR: {error}\n')

    w('\ncheck passed; this config would build\n' if report['ok'] else
      '\ncheck FAILED; the build would not start\n')


@dataclass
class RunSubcommand(Subcommand):
    name: str = 'build'
    group: str = 'Build a system'
    aliases: list = field(default_factory=lambda: ['run'])
    long_help: str = 'Prepare a system according to instructions in the config file.'
    short_help: str = 'prepare a system'
    func_returns_type: type = Controller

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('config', type=str, default=None, help='input configuration file in YAML format')
        self.parser.add_argument('--output-dir', type=str, default='./', help='name of output directory relative to CWD (default: %(default)s)')
        self.parser.add_argument('--gpu', default=False, action='store_true', help='force run on GPU')
        self.parser.add_argument('--seed', type=int, default=None, help='base NAMD RNG seed; use a different value per replica (0 = let NAMD draw from the clock)')
        self.parser.add_argument('--ncpus', type=int, default=0, help='number of NAMD processing elements (0 = auto-detect)')
        self.parser.add_argument('--complete-config', default=False, action='store_true', help='write complete config file')
        self.parser.add_argument('--restart', default=False, action='store_true',
                                 help='resume an interrupted build from the last cleanly-completed task '
                                      '(reads the .pestifer-manifest.json in the run directory)')
        self.parser.add_argument('--fresh', default=False, action='store_true',
                                 help='ignore any existing run manifest and build from scratch')
        self.parser.add_argument('--from', dest='from_task', type=str, default=None, metavar='TASK',
                                 help='resume explicitly from this task (index or taskname), overriding '
                                      'auto-detection (implies --restart)')
        self.parser.add_argument('--check', default=False, action='store_true',
                                 help='validate the config, the toolchain and the task pipeline, print '
                                      'the plan, and exit without building (exit 1 if it would not build)')
        self.parser.add_argument('--json', default=False, action='store_true',
                                 help='with --check, write the report to stdout as JSON')
        return self.parser

    def default_log_file(self, args):
        """
        Derive the diagnostics-log name from the input config, so that successive builds
        in one directory get unique logs (e.g. ``foo.yaml`` -> ``foo-diagnostics.log``)
        instead of clobbering a shared ``pestifer_diagnostics.log``.
        """
        config = getattr(args, 'config', None)
        if config:
            stem = os.path.splitext(os.path.basename(config))[0]
            if stem:
                return f'{stem}-diagnostics.log'
        return self.log_file

    @staticmethod
    def func(args, **kwargs):
        """
        Run the ``pestifer build`` command with the specified configuration file.
        """
        if args.output_dir != './':
            exec_dir = os.getcwd()
            if not os.path.exists(args.output_dir):
                os.mkdir(args.output_dir)
            if not os.path.exists(os.path.join(args.output_dir, args.config)):
                shutil.copy(args.config, args.output_dir)
            os.chdir(args.output_dir)
        # Set up the Controller and execute tasks
        begin_time = time.time()
        configname = args.config
        # include date and time in the message below
        logger.info(f'pestifer v. {__pestifer_version__} begins at {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(begin_time))} with config {configname}')

        allowed_extensions = ['.yaml', '.yml', '.y']
        cbase, cext = os.path.splitext(configname)
        if not cext:
            fil = [os.path.exists(f'{cbase}{ext}') for ext in allowed_extensions]
            if any(fil):
                iix = fil.index(True)
                configname = f'{cbase}{allowed_extensions[iix]}'

        if getattr(args, 'check', False):
            report = _preflight(configname, args, **kwargs)
            if getattr(args, 'json', False):
                print(json.dumps(report, indent=2))
            else:
                _print_check_report(report)
            if args.output_dir != './':
                os.chdir(exec_dir)
            return 0 if report['ok'] else 1

        config = Config(userfile=configname, ncpus_override=args.ncpus,
                        processor_type_override=('gpu' if getattr(args, 'gpu', False) else ''),
                        seed_override=getattr(args, 'seed', None),
                        **kwargs).configure_new()
        # Record the whole toolchain before anything runs, so this build's log is self-describing
        # and its results never have to have their environment reconstructed after the fact.
        config.environment_report = log_environment(config)
        C = Controller().configure(config)
        C.restart = getattr(args, 'restart', False)
        C.fresh = getattr(args, 'fresh', False)
        C.from_task = getattr(args, 'from_task', None)
        if args.complete_config:
            C.write_complete_config(f'{cbase}-complete.yaml')
        report = C.do_tasks()
        end_time = time.time()
        elapsed_time_s = datetime.timedelta(seconds=(end_time - begin_time))
        logger.info(f'pestifer ends. Elapsed time {time.strftime("%H:%M:%S", time.gmtime(elapsed_time_s.seconds))}.')
        logger.info(f'Task durations:')
        maxnamelen = max([len(t.taskname) for t in C.tasks]) if len(C.tasks) > 0 else 0
        name_format = f'{{:>{maxnamelen}s}}'
        for task in report.values():
            logger.info(f' - {name_format.format(task["taskname"])}: {hmsf(task["duration"])} ({task["duration_frac"]:>5.1%})')
        # Last thing the user sees: what this particular build owes credit to.  Placed after the
        # run rather than before it so it sits beside the results it applies to -- and only when
        # there are results: an aborted build has nothing to cite, and printing 'please cite'
        # under a stack of task failures reads as though it had succeeded.
        if not C.exit_code:
            cites = log_citations(config)
            # One machine-readable record of what this build actually did, for a Methods draft,
            # a replica comparison, or a summary table.  Written only for a build that finished:
            # a record of an aborted run would describe a system that does not exist.
            write_run_record(build_run_record(
                config, C.tasks,
                environment=getattr(config, 'environment_report', None),
                citations={'entries': [
                    {'subject': c.subject, 'text': c.text, 'doi': c.doi,
                     'reason': c.reason, 'key': c.key}
                    for c in (cites or [])]}))
        if args.output_dir != './':
            os.chdir(exec_dir)
        return C
