# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Version reporting must follow the source tree, not install metadata.

An editable install freezes the version in its distribution metadata at install time, so
after a release bump `pestifer --version` reported the old number until someone reinstalled
-- and CacheableObject.APP_VERSION keyed caches to it.
"""
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


if __name__ == '__main__':
    unittest.main()
