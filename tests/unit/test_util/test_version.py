# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Version reporting must follow the source tree, not install metadata.

An editable install freezes the version in its distribution metadata at install time, so
after a release bump `pestifer --version` reported the old number until someone reinstalled
-- and CacheableObject.APP_VERSION keyed caches to it.
"""
import re
import tomllib
import unittest
from pathlib import Path

from pestifer.util.stringthings import __pestifer_version__, _resolve_version


class TestVersionResolution(unittest.TestCase):

    def _pyproject_version(self):
        root = Path(__file__).resolve().parents[3]
        pyproject = root / 'pyproject.toml'
        if not pyproject.exists():
            self.skipTest('not running from a source tree')
        with open(pyproject, 'rb') as f:
            return tomllib.load(f)['project']['version']

    def test_reports_source_tree_version(self):
        self.assertEqual(_resolve_version(), self._pyproject_version())

    def test_module_constant_matches(self):
        self.assertEqual(__pestifer_version__, self._pyproject_version())

    def test_cli_and_library_agree(self):
        from pestifer.cli.pestifer import __pestifer_version__ as cli_version
        self.assertEqual(cli_version, __pestifer_version__)


class TestSingleVersionSource(unittest.TestCase):
    """No module may look up pestifer's own version for itself.

    The first fix routed the CLI through the resolver but left two other call sites doing
    their own ``importlib.metadata.version('pestifer')``: the startup INFO banner
    (``core/config.py``) printed 3.15.1 beside a 3.16.1 splash, and the run manifest
    (``core/controller.py``) stamped that stale number into every build's provenance record.
    The per-site tests above could not see them, because the defect is a *second* source of
    truth rather than a wrong one -- so this guard is source-level.
    """

    #: the resolver itself is the one legitimate caller
    ALLOWED = {Path('util/stringthings.py')}

    def test_no_module_resolves_its_own_version(self):
        import pestifer
        root = Path(pestifer.__file__).parent
        pattern = re.compile(r"""version\(\s*['"]pestifer['"]\s*\)""")
        offenders = []
        for src in sorted(root.rglob('*.py')):
            rel = src.relative_to(root)
            # _archive holds retired code that is never imported
            if rel in self.ALLOWED or rel.parts[:2] == ('resources', '_archive'):
                continue
            for n, line in enumerate(src.read_text().splitlines(), 1):
                if pattern.search(line):
                    offenders.append(f'{rel}:{n}: {line.strip()}')
        self.assertEqual(
            offenders, [],
            'these must import __pestifer_version__ from pestifer.util.stringthings '
            'instead of resolving the version themselves:\n  ' + '\n  '.join(offenders))


if __name__ == '__main__':
    unittest.main()
