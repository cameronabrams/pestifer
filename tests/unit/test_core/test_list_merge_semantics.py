# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""No shipped example has a list it wrote silently reordered by schema resolution.

ycleptic merges a user list into the schema default with `special_update`, which dedupes and
**keeps the DEFAULT's ordering**.  So a user who writes ``[B, A]`` for a list whose default
contains ``[A, B]`` gets ``[A, B]`` back: not extended, reordered.

Appending is usually right -- a user's extra `psfgen.aliases` should arrive alongside the
standard set, not instead of it -- so this does not forbid growth.  It forbids the user's own
entries coming back in a different relative order, which nothing would report.

It matters most for `charmmff.standard.str`, where the list IS topology load order: CLAUDE.md
records that a residue's terminal patch is decided by which file was read last, `GLYM`
inheriting whatever the previously loaded file left behind because
`toppar_all36_lipid_prot.str` has its `DEFA` line commented out.  A user reordering stream
files to control that would find the order silently reverted and the symptom only in the built
PSF.

Measured when this was written: 11 lists are appended to, 0 reordered.  Pure config
resolution -- no VMD, no NAMD -- so it runs on CI.
"""
import unittest
from importlib.resources import files as pkg_files
from pathlib import Path

import yaml

from pestifer.core.config import Config

_EXAMPLES = Path(str(pkg_files('pestifer.resources'))) / 'examples'


def _lists(node, path=()):
    if isinstance(node, dict):
        for k, v in node.items():
            yield from _lists(v, path + (k,))
    elif isinstance(node, list):
        yield path, node


def _get(node, path):
    for k in path:
        if not isinstance(node, dict) or k not in node:
            return None
        node = node[k]
    return node


def _is_subsequence(small, big):
    it = iter(big)
    return all(any(x == y for y in it) for x in small)


class TestNoWrittenListIsReordered(unittest.TestCase):

    def test_every_example(self):
        files = sorted(_EXAMPLES.glob('*/inputs/*.yaml'))
        self.assertGreater(len(files), 25, f'no example configs under {_EXAMPLES}')
        for f in files:
            with self.subTest(example=f.parent.parent.name):
                raw = yaml.safe_load(f.read_text()) or {}
                resolved = Config(userdict=yaml.safe_load(f.read_text()),
                                  quiet=True).configure_new()['user']
                for path, written in _lists(raw):
                    if path and path[0] == 'tasks':
                        continue          # the task list is positional, not merged
                    got = _get(resolved, path)
                    if got is None or got == written:
                        continue
                    self.assertTrue(
                        _is_subsequence(written, got),
                        f'{".".join(map(str, path))}: the entries written were reordered\n'
                        f'   wrote    {written}\n   resolved {got}')

    def test_the_check_can_fail(self):
        """A reordered list must be detected -- the whole test is a subsequence check, so a
        bug in it would pass everything silently."""
        self.assertTrue(_is_subsequence(['b', 'c'], ['a', 'b', 'c']))
        self.assertFalse(_is_subsequence(['c', 'b'], ['a', 'b', 'c']),
                         'reordering is not being detected')
