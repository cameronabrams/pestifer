# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for how subcommands are grouped and presented.

Nineteen commands presented as one flat list is a list nobody reads: a user looking for "how do I
look at my density curve" has to scan all of them.  They are therefore grouped by *what you hand
them* -- a config, a structure, a NAMD run's output, or nothing at all -- which is a rule a reader
can apply to an unfamiliar command rather than a list they have to memorize.

What is worth testing is not the grouping's aesthetics but its completeness and its agreement with
the documentation.  A command added without a group, or a group heading that exists in the CLI but
not in the docs, reintroduces exactly the flat list this replaced -- silently, because nothing
fails.  The two lists had already drifted into different arbitrary orders once.
"""

import os
import re
import unittest

from pestifer.subcommands import GROUPS, _subcommands, grouped_subcommands

_HERE = os.path.dirname(os.path.abspath(__file__))
DOCS = os.path.normpath(os.path.join(_HERE, '..', '..', '..', 'docs', 'source', 'usage.rst'))


class TestEverySubcommandIsGrouped(unittest.TestCase):

    def test_no_command_is_left_ungrouped(self):
        """The failure this catches is silent: an ungrouped command still works, it just stops
        being findable."""
        self.assertEqual([s.name for s in _subcommands if not s.group], [])

    def test_every_group_used_is_a_declared_one(self):
        """A typo in a group name would otherwise create a new group of one."""
        self.assertEqual(sorted({s.group for s in _subcommands} - set(GROUPS)), [])

    def test_every_declared_group_has_members(self):
        empty = [g for g in GROUPS if not any(s.group == g for s in _subcommands)]
        self.assertEqual(empty, [])

    def test_the_grouping_partitions_the_commands(self):
        grouped = [s for _g, members in grouped_subcommands() for s in members]
        self.assertEqual(len(grouped), len(_subcommands))
        self.assertEqual({s.name for s in grouped}, {s.name for s in _subcommands})


class TestPresentationOrder(unittest.TestCase):

    def test_groups_are_presented_in_the_declared_order(self):
        self.assertEqual([g for g, _m in grouped_subcommands()],
                         [g for g in GROUPS if any(s.group == g for s in _subcommands)])

    def test_the_registry_itself_is_in_group_order(self):
        """The CLI walks ``_subcommands`` directly, so the list order is the display order."""
        seen = [GROUPS.index(s.group) for s in _subcommands]
        self.assertEqual(seen, sorted(seen))

    def test_building_a_system_comes_first(self):
        self.assertEqual(_subcommands[0].name, 'build')

    def test_the_run_analysis_group_is_where_the_seam_is(self):
        """Four of these five read NAMD output and nothing of pestifer's, so they work on runs
        pestifer never managed; report-methods needs a run-record.json and does not.  Keeping them
        in one named group is what makes that boundary visible without splitting the package."""
        after = dict(grouped_subcommands())['After the run']
        self.assertIn('report-methods', [s.name for s in after])
        self.assertIn('make-namd-restart', [s.name for s in after])


class TestDocsAgreeWithTheCli(unittest.TestCase):
    """The two lists disagreed before -- different orders, neither grouped -- so neither was
    learnable.  These are what stop them drifting apart again."""

    @classmethod
    def setUpClass(cls):
        with open(DOCS) as f:
            cls.usage = f.read()

    def test_the_usage_page_was_actually_found(self):
        # a path typo would make every test below vacuously pass
        self.assertIn('Subcommands', self.usage)

    def test_every_group_heading_appears_in_the_docs(self):
        for g in GROUPS:
            with self.subTest(group=g):
                self.assertIn(g, self.usage)

    def test_the_docs_present_the_groups_in_the_same_order(self):
        positions = [self.usage.index(g) for g in GROUPS]
        self.assertEqual(positions, sorted(positions))

    def test_no_command_page_is_listed_above_the_first_group(self):
        """A page outside any group section renders but is unreachable by the rule."""
        first_group = self.usage.index(GROUPS[0])
        before = re.findall(r'^\s+subs/([a-z0-9-]+)\s*$', self.usage[:first_group], re.M)
        self.assertEqual(before, [])

    def test_every_command_with_a_doc_page_is_listed(self):
        pages = {f[:-4] for f in os.listdir(os.path.join(os.path.dirname(DOCS), 'subs'))
                 if f.endswith('.rst')}
        listed = set(re.findall(r'^\s+subs/([a-z0-9-]+)\s*$', self.usage, re.M))
        self.assertEqual(pages - listed, set(), 'documented but not listed in usage.rst')

    def test_fetch_example_is_accounted_for_despite_having_no_page(self):
        # documented inside build-example; without the note it looks simply missing
        self.assertIn('fetch-example', self.usage)


class TestHelpFormatter(unittest.TestCase):
    """The headings are injected as each group's first command is formatted, because argparse
    allows only one subparsers title for the whole list."""

    def _render(self, subcommands=None):
        import argparse as ap
        from pestifer.cli.pestifer import grouped_formatter
        subcommands = _subcommands if subcommands is None else subcommands
        parser = ap.ArgumentParser(prog='pestifer',
                                   formatter_class=grouped_formatter(subcommands))
        subparsers = parser.add_subparsers(title='Available commands', metavar='')
        for s in subcommands:
            subparsers.add_parser(s.name, aliases=s.aliases, help=s.short_help)
        return parser.format_help()

    def test_each_group_heading_is_printed_once(self):
        text = self._render()
        for g in GROUPS:
            with self.subTest(group=g):
                self.assertEqual(text.count(g), 1)

    def test_the_headings_precede_their_commands(self):
        text = self._render()
        for g, members in grouped_subcommands():
            self.assertLess(text.index(g), text.index(members[0].name),
                            f'{g} heading did not precede {members[0].name}')

    def test_every_command_still_appears(self):
        text = self._render()
        for s in _subcommands:
            with self.subTest(command=s.name):
                self.assertIn(s.name, text)

    def test_no_stray_blank_line_under_the_title(self):
        """argparse formats the subparsers action itself as an empty invocation line."""
        after_title = self._render().split('Available commands', 1)[1].split('\n', 1)[1]
        self.assertTrue(after_title.startswith(f'  {GROUPS[0]}'),
                        f'unexpected content under the title: {after_title[:80]!r}')

    def test_an_untagged_command_set_renders_ungrouped_rather_than_failing(self):
        """The formatter must degrade to plain argparse behavior if nothing carries a group."""
        import argparse as ap
        from pestifer.cli.pestifer import grouped_formatter
        parser = ap.ArgumentParser(prog='pestifer', formatter_class=grouped_formatter([]))
        sub = parser.add_subparsers(title='Available commands', metavar='')
        sub.add_parser('build', help='prepare a system')
        text = parser.format_help()
        self.assertIn('build', text)
        for g in GROUPS:
            self.assertNotIn(g, text)


if __name__ == '__main__':
    unittest.main()
