# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Guards on the repository's own shape -- things that make the suite behave differently on a fresh
clone than on a machine that has been running it for a while.  Pure filesystem and ``git``; no
VMD, no NAMD, so this runs on CI, which is the only place it can catch anything.
"""
import subprocess
import unittest
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
_UNIT = Path(__file__).resolve().parent


def _tracked(path: Path) -> list[str]:
    """Files git tracks under ``path`` (empty list if it tracks none)."""
    out = subprocess.run(['git', 'ls-files', '--', str(path)],
                         cwd=_REPO, capture_output=True, text=True)
    return [ln for ln in out.stdout.splitlines() if ln.strip()]


class TestModuleWorkingDirectoriesAreTracked(unittest.TestCase):
    """``conftest.py``'s ``change_test_dir`` runs a module inside ``<module_name>/`` **if that
    directory exists**, and otherwise one level up, in the module's own directory.  Every relative
    path a module writes is written for one of those two, so a directory that exists on a
    developer machine and not in a clone silently shifts the whole module by one level.

    That has now happened three times here, twice on 2026-10-14
    (``test_tasks/test_make_membrane_system/`` and ``test_core/test_tcl/``), and it is invisible
    to ``git status``: the directory holds only generated output, so nothing in it is tracked and
    nothing reports it missing.  A single tracked file -- a ``.gitkeep`` is enough -- makes the
    directory part of the repository and the two cases identical.

    So: a module-named directory must either be absent, or carry at least one tracked file.
    """

    def test_every_module_working_directory_has_a_tracked_file(self):
        offenders = []
        for module in sorted(_UNIT.rglob('test_*.py')):
            workdir = module.with_suffix('')
            if workdir.is_dir() and not _tracked(workdir):
                offenders.append(str(workdir.relative_to(_REPO)))
        self.assertEqual(offenders, [],
                         'module working directories with no tracked file (add a .gitkeep): '
                         + ', '.join(offenders))
