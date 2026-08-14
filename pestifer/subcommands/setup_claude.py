# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The setup-claude subcommand.  Installs pestifer's Claude Code skill into
~/.claude/skills/pestifer/ so that an agent driving pestifer knows how to use it.
"""
import argparse as ap
import shutil
from dataclasses import dataclass
from importlib.resources import files as pkg_files
from pathlib import Path

from . import Subcommand
from ..core.errors import PestiferError

_SKILL_NAME = 'pestifer'


@dataclass
class SetupClaudeSubcommand(Subcommand):
    name: str = 'setup-claude'
    short_help: str = "install pestifer's Claude Code skill so an agent can drive pestifer"
    long_help: str = (
        "Copies pestifer's bundled skill to ~/.claude/skills/pestifer/SKILL.md, which teaches "
        "Claude Code how to use pestifer: check configs before building, start from a worked "
        "example, run builds in the background, and recover from a failed run.  Re-run after "
        "upgrading pestifer to pick up the current version of the skill."
    )

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        # Resolve through pestifer.resources (a real package) rather than treating the
        # claude/ subdirectory as an importable one -- it carries no __init__.py, matching
        # every other resource subdirectory.
        source = pkg_files('pestifer.resources').joinpath('claude', 'SKILL.md')
        if not source.is_file():
            raise PestiferError(f'Bundled skill not found at {source}')

        skill_dir = Path(args.skills_dir).expanduser() / _SKILL_NAME
        target = skill_dir / 'SKILL.md'

        if target.exists() and not args.force:
            print(f'{target} already exists; re-run with --force to overwrite it')
            return True

        skill_dir.mkdir(parents=True, exist_ok=True)
        with source.open('rb') as fh:
            with open(target, 'wb') as out:
                shutil.copyfileobj(fh, out)
        print(f'Written: {target}')
        print('Claude Code picks the skill up on its next session in any directory.')
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('--skills-dir', type=str, default='~/.claude/skills',
                                 help='skills directory to install into (default: %(default)s); '
                                      'use ./.claude/skills to scope the skill to one project')
        self.parser.add_argument('--force', default=False, action='store_true',
                                 help='overwrite an existing installed skill')
        return self.parser
