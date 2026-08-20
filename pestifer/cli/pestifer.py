# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Defines the command-line interface for pestifer
"""
import argparse as ap
import logging
import os
import shutil
import sys

from ..util.stringthings import banner, __pestifer_version__
from ..subcommands import _subcommands
from ..core.errors import PestiferError

logger = logging.getLogger(__name__)

logging.getLogger("pidibble").setLevel(logging.WARNING)
logging.getLogger("ycleptic").setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("filelock").setLevel(logging.WARNING)

class NiceHelpFormatter(ap.HelpFormatter):
    """
    A custom help formatter that provides a more readable help output.
    It formats the help text with a maximum help position and width.

    It also prints the command list under group headings.  ``argparse`` allows only one
    ``add_subparsers`` call and therefore only one title for the whole list, so the headings are
    injected as each group's first command is formatted.  ``group_of`` maps a command name to its
    heading and is set by :func:`grouped_formatter`; with it empty this behaves exactly as before.
    """

    group_of: dict = {}

    def __init__(self, prog):
        super().__init__(prog, max_help_position=40, width=100)
        self._groups_emitted = set()

    def _format_action(self, action):
        text = super()._format_action(action)
        if isinstance(action, ap._SubParsersAction):
            # The subparsers action's own invocation line is blank -- its metavar is empty so
            # argparse does not print the whole {build,run,...} choice list -- and it heads the
            # block containing every command, so drop it rather than leave a stray blank line
            # between the title and the first group.
            lines = text.split('\n')
            while lines and not lines[0].strip():
                lines.pop(0)
            return '\n'.join(lines)
        group = self.group_of.get(getattr(action, 'dest', None))
        if group and group not in self._groups_emitted:
            first = not self._groups_emitted
            self._groups_emitted.add(group)
            text = f'{"" if first else chr(10)}  {group}\n{text}'
        return text


def grouped_formatter(subcommands):
    """A :class:`NiceHelpFormatter` subclass that knows which group each command belongs to.

    Returned as a class because ``argparse`` instantiates the formatter itself, and baked into a
    subclass rather than set on the shared base so nothing leaks between parsers.
    """
    return type('GroupedHelpFormatter', (NiceHelpFormatter,),
                {'group_of': {s.name: s.group for s in subcommands if s.group}})

def subcommand_exit_code(result) -> int:
    """Exit status implied by whatever a subcommand returned.

    A subcommand that ran tasks reports failure either as an int or via an ``exit_code``
    attribute on the object it returns (the Controller); without honoring that, the process
    exited 0 for an aborted build and any caller keying off the exit status recorded it as a
    success.

    Most subcommands, though, return ``True`` to mean "done" -- and ``bool`` is a subclass of
    ``int``, so treating a bare int as a status turned every one of those into exit 1.  Booleans
    are therefore success, and only a genuine int (or a genuine int ``exit_code``) is a status.
    """
    if result is None or isinstance(result, bool):
        return 0
    if isinstance(result, int):
        return result
    code = getattr(result, 'exit_code', 0)
    if isinstance(code, bool) or not isinstance(code, int):
        return 0
    return code


def cli():
    """
    Command-line interface for pestifer.
    """

    parser = ap.ArgumentParser(formatter_class=grouped_formatter(_subcommands))
    parser.add_argument(
        "--banner",
        default=True,
        action=ap.BooleanOptionalAction,
        help="enable or disable the banner"
    )
    parser.add_argument(
        '--kick-ass',
        default=False,
        action=ap.BooleanOptionalAction,
        help=ap.SUPPRESS)
    parser.add_argument(
        '--log-level',
        type=str,
        default='debug',
        choices=['info', 'debug', 'warning'],
        help='logging level (default: %(default)s)')
    parser.add_argument(
        '--log-file',
        type=str,
        default=None,
        help="log file (default: 'pestifer_diagnostics.log'; build/run derive it from the config name)")
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {__pestifer_version__}'
    )
    subparsers = parser.add_subparsers(
        title="Available commands (use \"pestifer <command> --help\" for help with any command)",
        dest="command",
        metavar="",
        required=False
    )

    for subcommand in _subcommands:
        subcommand.add_subparser(subparsers)

    args = parser.parse_args()
    loglevel_numeric = getattr(logging, args.log_level.upper())
    default_log_file_func = getattr(args, 'default_log_file_func', None)
    log_file = args.log_file or (default_log_file_func(args) if default_log_file_func else None)
    if log_file:
        if os.path.exists(log_file):
            shutil.copyfile(log_file, log_file + '.bak')
        logging.basicConfig(filename=log_file, filemode='w', format='%(asctime)s %(name)s %(message)s', level=loglevel_numeric)
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter('%(levelname)s> %(message)s'))
    logging.getLogger('').addHandler(console)
    # Emit the banner on the same stream the console log handler uses (stderr), not stdout.
    # Under `pestifer ... > run.log 2>&1` the two streams share one file, but a redirected
    # stdout is block-buffered and only flushes at exit, so a banner printed to stdout lands
    # at the *end* of the log -- or mid-file once its buffer fills -- even though it is written
    # first.  stderr stays line-buffered when redirected, so writing there keeps the banner at
    # the top where it belongs.
    banner(lambda line: print(line, file=console.stream, flush=True), args)
    try:
        result = args.func(args)
    except PestiferError as e:
        logger.error(str(e))
        sys.exit(1)
    exit_code = subcommand_exit_code(result)
    if exit_code:
        logger.error(f'pestifer exiting with code {exit_code}')
        sys.exit(exit_code)
