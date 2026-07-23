import unittest

from pestifer.core.system_inspector import (
    build_findings, _group_runs, Findings, MissingRun, MutationFinding, ExcisionRun,
    interactive_select, make_prompter, ChainIdentity, AssemblyInfo,
)


class TestGroupRuns(unittest.TestCase):
    def test_groups_consecutive(self):
        self.assertEqual(list(_group_runs([1, 2, 3, 7, 8, 20])), [(1, 3), (7, 8), (20, 20)])

    def test_dedups_and_sorts(self):
        self.assertEqual(list(_group_runs([5, 4, 5, 6])), [(4, 6)])

    def test_empty(self):
        self.assertEqual(list(_group_runs([])), [])


class TestClassification(unittest.TestCase):
    def _findings(self, missing, resolved, seqadv=()):
        return build_findings('XXXX', 'pdb', missing, list(seqadv), resolved)

    def test_terminal_vs_interior(self):
        # chain A resolved 10..100; missing 5-7 (before -> N), 40-45 (interior), 105-108 (after -> C)
        missing = [('A', n) for n in (5, 6, 7, 40, 41, 42, 43, 44, 45, 105, 106, 107, 108)]
        f = self._findings(missing, {'A': (10, 100)})
        kinds = {(r.start, r.end): r.kind for r in f.missing_runs}
        self.assertEqual(kinds[(5, 7)], 'N')
        self.assertEqual(kinds[(40, 45)], 'interior')
        self.assertEqual(kinds[(105, 108)], 'C')
        self.assertEqual(f.terminal_tail_chains(), {'n': ['A'], 'c': ['A']})
        self.assertEqual(len(f.interior_gaps()), 1)

    def test_no_resolved_range_is_interior(self):
        f = self._findings([('B', 1), ('B', 2)], {})     # no ATOM records for B
        self.assertEqual(f.missing_runs[0].kind, 'interior')

    def test_seqadv_mutations_and_excisions(self):
        seqadv = [
            ('A', 52, 'ASN', 'THR', 'engineered mutation'),
            ('A', 90, 'CYS', 'ALA', 'conflict'),
            ('A', 200, 'HIS', '', 'expression tag'),
            ('A', 201, 'HIS', '', 'expression tag'),
            ('A', 202, 'HIS', '', 'expression tag'),
            ('A', 5, 'GLY', '', 'microheterogeneity'),   # ignored type
        ]
        f = self._findings([], {'A': (1, 300)}, seqadv)
        self.assertEqual(len(f.mutations), 2)
        self.assertEqual(f.mutations[0].revert_shortcode, 'A:ASN,52,THR')
        self.assertEqual(len(f.excisions), 1)                      # 200-202 grouped into one run
        self.assertEqual(f.excisions[0].delete_shortcode, 'A:200-202')

    def test_missing_run_length_and_str(self):
        r = MissingRun('A', 10, 10, 'interior')
        self.assertEqual(r.length, 1)
        self.assertIn('10 (1 residue)', str(r))


class TestAnnotation(unittest.TestCase):
    def test_annotation_has_correct_nesting_guidance(self):
        f = Findings('XXXX', 'pdb',
                     missing_runs=[MissingRun('A', 1, 3, 'N'), MissingRun('A', 50, 60, 'interior')],
                     mutations=[MutationFinding('A', 52, 'ASN', 'THR', 'engineered mutation')],
                     excisions=[ExcisionRun('A', 200, 202, 'expression tag')])
        text = '\n'.join(f.annotation_lines())
        # every emitted line is a comment
        self.assertTrue(all(l.strip().startswith('#') for l in f.annotation_lines()))
        # sequence goes under source; mods is a sibling of source
        self.assertIn('terminal_tails:', text)
        self.assertIn('n: [A]', text)
        self.assertIn('mods:', text)
        self.assertIn('mutations: [A:ASN,52,THR]', text)
        self.assertIn('deletions: [A:200-202]', text)
        self.assertIn('SIBLING of `source:`', text)

    def test_empty_findings(self):
        f = Findings('XXXX', 'pdb')
        self.assertTrue(f.is_empty())


class TestInteractiveSelect(unittest.TestCase):
    def _findings(self):
        return Findings('X', 'pdb',
                        missing_runs=[MissingRun('A', 1, 3, 'N'), MissingRun('B', 5, 5, 'C'),
                                      MissingRun('A', 50, 60, 'interior')],
                        mutations=[MutationFinding('A', 52, 'ASN', 'THR', 'engineered mutation'),
                                   MutationFinding('A', 90, 'CYS', '', 'conflict')],  # no dbres
                        excisions=[ExcisionRun('A', 200, 202, 'expression tag')])

    def test_accept_all(self):
        sel = interactive_select(self._findings(), ask=lambda q, d=False: True, say=lambda m: None)
        self.assertEqual(sel['sequence']['terminal_tails'], {'n': ['A'], 'c': ['B']})
        self.assertTrue(sel['add_ligate'])
        self.assertEqual(sel['mods']['mutations'], ['A:ASN,52,THR'])   # only the revertible one
        self.assertEqual(sel['mods']['deletions'], ['A:200-202'])

    def test_decline_all(self):
        sel = interactive_select(self._findings(), ask=lambda q, d=False: False, say=lambda m: None)
        self.assertEqual(sel['sequence'], {})
        self.assertEqual(sel['mods'], {})
        self.assertFalse(sel['add_ligate'])

    def test_selective(self):
        # accept only tail and ligate prompts
        ask = lambda q, d=False: ('tail' in q) or ('ligate' in q)
        sel = interactive_select(self._findings(), ask=ask, say=lambda m: None)
        self.assertIn('terminal_tails', sel['sequence'])
        self.assertTrue(sel['add_ligate'])
        self.assertNotIn('mutations', sel['mods'])
        self.assertNotIn('deletions', sel['mods'])

    def test_prompter_default_on_blank(self):
        import builtins
        orig = builtins.input
        builtins.input = lambda prompt='': ''
        try:
            ask = make_prompter()
            self.assertTrue(ask('q?', True))
            self.assertFalse(ask('q?', False))
        finally:
            builtins.input = orig

    def test_prompter_yes(self):
        import builtins
        orig = builtins.input
        builtins.input = lambda prompt='': 'yes'
        try:
            self.assertTrue(make_prompter()('q?', False))
        finally:
            builtins.input = orig


class TestChainAndAssembly(unittest.TestCase):
    def test_chain_describe(self):
        self.assertEqual(ChainIdentity('A', 'protein', 187, ['ALA', 'GLY']).describe(),
                         'protein (187 residues)')
        self.assertEqual(ChainIdentity('W', 'water', 50, ['HOH']).describe(), 'water')
        self.assertTrue(ChainIdentity('A', 'glycan', 3, ['NAG', 'BMA']).describe().startswith('glycan (NAG'))
        self.assertTrue(ChainIdentity('I', 'ion', 1, ['ZN']).describe().startswith('ion (ZN'))
        self.assertIn('nucleic acid', ChainIdentity('T', 'nucleicacid', 12, ['DA']).describe())

    def test_is_empty_accounts_for_multichain(self):
        two = Findings('X', 'pdb', chains=[ChainIdentity('A', 'protein', 10),
                                           ChainIdentity('B', 'protein', 10)])
        self.assertFalse(two.is_empty())
        one = Findings('X', 'pdb', chains=[ChainIdentity('A', 'protein', 10)])
        self.assertTrue(one.is_empty())

    def test_interactive_assembly_and_chain_omit(self):
        f = Findings('X', 'pdb',
                     chains=[ChainIdentity('A', 'protein', 100), ChainIdentity('W', 'water', 20, ['HOH'])],
                     assemblies=[AssemblyInfo(1, 3, ['A', 'W'])])
        ask = lambda q, d=False: ('biological assembly' in q) or ('Omit chain W' in q)
        sel = interactive_select(f, ask=ask, say=lambda m: None)
        self.assertEqual(sel['source']['biological_assembly'], 1)
        self.assertEqual(sel['source']['exclude'], ["chainID in ['W']"])

    def test_interactive_declining_assembly_gives_asymmetric_unit(self):
        f = Findings('X', 'pdb', chains=[ChainIdentity('A', 'protein', 100)],
                     assemblies=[AssemblyInfo(1, 3, ['A'])])
        sel = interactive_select(f, ask=lambda q, d=False: False, say=lambda m: None)
        self.assertEqual(sel['source']['biological_assembly'], 0)

    def test_interactive_no_assemblies_is_asymmetric_unit(self):
        f = Findings('X', 'pdb', chains=[ChainIdentity('A', 'protein', 100)])
        sel = interactive_select(f, ask=lambda q, d=False: True, say=lambda m: None)
        self.assertEqual(sel['source']['biological_assembly'], 0)

    def test_annotation_lists_assemblies_and_chains(self):
        f = Findings('X', 'pdb',
                     chains=[ChainIdentity('A', 'protein', 100), ChainIdentity('B', 'glycan', 3, ['NAG'])],
                     assemblies=[AssemblyInfo(1, 3, ['A', 'B'])])
        text = '\n'.join(f.annotation_lines())
        self.assertIn('Biological assemblies', text)
        self.assertIn('1: 3 copies of chain(s) [A, B]', text)
        self.assertIn('A: protein (100 residues)', text)
        self.assertIn('exclude:', text)


if __name__ == '__main__':
    unittest.main()
