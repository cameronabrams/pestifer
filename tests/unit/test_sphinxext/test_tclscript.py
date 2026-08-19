# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the ``tclscript`` Sphinx directive.

The directive documents a Tcl script that ships in the repository: it locates the file relative to
the repo root, lifts the ``##``-prefixed header comment out of it, parses that as RST, and appends
a link to the file on GitHub.

Two properties make it worth pinning.  It resolves paths against the **repo root computed from
this package's own location**, which is what lets an rst file reference a script by repo-relative
path -- and is easy to break by moving the module.  And it fails *quietly by design*: a missing
script emits a build warning and renders nothing, deliberately, because the previous behavior
rendered ``ERROR: File not found: <builder's absolute path>`` into the published page and a
deleted script stayed referenced long after it was gone.  A test that did not assert the warning
would not notice if that warning stopped being emitted, leaving a silently-empty section.

The repo root is patched rather than written to, so these tests never touch the working tree.
"""

import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

from docutils import nodes

from pestifer.sphinxext import tclscript
from pestifer.sphinxext.tclscript import TclScriptDirective


SCRIPT = """\
## My Script
##
## Does a useful thing.

# an ordinary comment, not part of the header
puts "hello"
## this trailing block is past the first code line and must not be picked up
"""


class _DirectiveCase(unittest.TestCase):
    """Runs the directive against a throwaway repo root."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = self._tmp.name
        self.addCleanup(self._tmp.cleanup)
        self.dependencies = []
        self.parsed = None
        self.warnings = []

    def write_script(self, relpath, text=SCRIPT):
        full = os.path.join(self.root, relpath)
        os.makedirs(os.path.dirname(full), exist_ok=True)
        with open(full, 'w') as f:
            f.write(text)
        return full

    def run_directive(self, relpath):
        d = TclScriptDirective.__new__(TclScriptDirective)
        d.arguments = [relpath]
        d.lineno = 7
        env = SimpleNamespace(docname='page',
                              note_dependency=self.dependencies.append)
        d.state = SimpleNamespace(
            document=SimpleNamespace(settings=SimpleNamespace(env=env)))

        def fake_parse(state, viewlist, node):
            self.parsed = list(viewlist)

        # The repo root is derived from this module's __file__ at call time, so patching the
        # module attribute redirects it without touching the real tree.
        fake_file = os.path.join(self.root, 'pestifer', 'sphinxext', 'tclscript.py')
        with mock.patch.object(tclscript, '__file__', fake_file), \
             mock.patch.object(tclscript, 'nested_parse_with_titles', fake_parse), \
             mock.patch.object(tclscript.logger, 'warning',
                               lambda msg, **kw: self.warnings.append(msg)):
            return d.run()


class TestMissingScript(_DirectiveCase):

    def test_nothing_is_rendered(self):
        self.assertEqual(self.run_directive('scripts/gone.tcl'), [],
                         'a missing script must render nothing, not an error message')

    def test_a_build_warning_is_emitted_naming_the_path(self):
        self.run_directive('scripts/gone.tcl')
        self.assertEqual(len(self.warnings), 1)
        self.assertIn('gone.tcl', self.warnings[0])
        self.assertIn('not found', self.warnings[0])

    def test_a_missing_script_is_not_registered_as_a_dependency(self):
        self.run_directive('scripts/gone.tcl')
        self.assertEqual(self.dependencies, [])


class TestHeaderExtraction(_DirectiveCase):

    def header(self, text):
        self.write_script('scripts/s.tcl', text)
        self.run_directive('scripts/s.tcl')
        return self.parsed

    def test_double_hash_lines_become_the_header_text(self):
        got = self.header(SCRIPT)
        self.assertIn('My Script', got)
        self.assertIn('Does a useful thing.', got)

    def test_the_comment_markers_are_stripped(self):
        for line in self.header(SCRIPT):
            self.assertFalse(line.startswith('#'), f'unstripped marker in {line!r}')

    def test_blank_and_single_hash_lines_become_blank(self):
        # a single-# line is an ordinary code comment, not prose: it must not appear as text
        got = self.header(SCRIPT)
        self.assertNotIn('an ordinary comment, not part of the header', ' '.join(got))
        self.assertIn('', got)

    def test_the_header_ends_at_the_first_line_of_code(self):
        got = ' '.join(self.header(SCRIPT))
        self.assertNotIn('trailing block', got,
                         'comments after the first code line are not part of the header')

    def test_a_script_with_no_header_is_still_rendered(self):
        self.write_script('scripts/bare.tcl', 'puts "hi"\n')
        result = self.run_directive('scripts/bare.tcl')
        self.assertIsNone(self.parsed, 'nothing to parse when there is no header')
        self.assertEqual(len(result), 1, 'the section should still be produced')


class TestRenderedSection(_DirectiveCase):

    def section(self, relpath='scripts/my_script.tcl'):
        self.write_script(relpath)
        return self.run_directive(relpath)[0]

    def test_the_section_is_titled_and_anchored_by_basename(self):
        sec = self.section()
        self.assertIsInstance(sec, nodes.section)
        self.assertEqual(sec['ids'], ['my_script.tcl'])
        self.assertEqual(sec.children[0].astext(), 'my_script.tcl')

    def test_a_source_link_is_appended(self):
        refs = [n for n in self.section().findall(nodes.reference)]
        self.assertEqual(len(refs), 1)
        self.assertEqual(refs[0].astext(), '[source]')

    def test_the_link_points_at_the_file_on_github(self):
        uri = [n for n in self.section().findall(nodes.reference)][0]['refuri']
        self.assertEqual(
            uri,
            'https://github.com/cameronabrams/pestifer/blob/main/scripts/my_script.tcl')

    def test_the_link_uses_forward_slashes_regardless_of_platform(self):
        uri = [n for n in self.section(
            'scripts/sub/deep.tcl').findall(nodes.reference)][0]['refuri']
        self.assertIn('scripts/sub/deep.tcl', uri)
        self.assertNotIn('\\', uri)


class TestRebuildDependency(_DirectiveCase):
    """Without this the page goes stale: editing the script would not trigger a rebuild."""

    def test_the_script_is_registered_as_a_dependency(self):
        full = self.write_script('scripts/s.tcl')
        self.run_directive('scripts/s.tcl')
        self.assertEqual(self.dependencies, [full])


class TestPathResolution(_DirectiveCase):

    def test_paths_are_resolved_against_the_repo_root(self):
        self.write_script('scripts/s.tcl')
        self.assertEqual(len(self.run_directive('scripts/s.tcl')), 1)

    def test_a_redundant_path_is_normalized(self):
        self.write_script('scripts/s.tcl')
        self.assertEqual(len(self.run_directive('./scripts/../scripts/s.tcl')), 1)


if __name__ == '__main__':
    unittest.main()
