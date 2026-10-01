# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Guards on :mod:`pestifer.tasks.validate`.

The failure these pin down is quiet: a test that emits no check is indistinguishable, in the
VMD log and in the ``pass: N, fail: M`` tally, from a test that ran and passed.  A typo in
``measure:`` therefore used to *shrink the test count* while the build still reported success.
So an unsupported ``measure`` or ``relation`` must be rejected -- at config load by the schema,
and at construction by the class -- and never merely logged.
"""
import unittest

from pestifer.tasks.validate import AttributeTest, ConnectionTest, ResidueTest


class _Recorder:
    """Stands in for a VMDScripter; collects the lines a test writes."""
    def __init__(self):
        self.lines = []

    def addline(self, line):
        self.lines.append(line)

    @property
    def text(self):
        return '\n'.join(self.lines)


class TestResidueTestRejectsUnsupportedSpecs(unittest.TestCase):
    def test_supported_measures_emit_a_check(self):
        for measure in ('atom_count', 'residue_count'):
            r = _Recorder()
            ResidueTest('t', 'protein', measure, 6).write(r)
            self.assertIn('vmdcon "PASS', r.text)
            self.assertIn('vmdcon "FAIL', r.text)

    def test_unsupported_measure_is_rejected_at_construction(self):
        with self.assertRaises(NotImplementedError) as cm:
            ResidueTest('t', 'protein', 'residue_cont', 6)          # plausible typo
        self.assertIn('residue_cont', str(cm.exception))
        self.assertIn('residue_count', str(cm.exception))           # names the alternatives

    def test_unsupported_relation_is_rejected_at_construction(self):
        with self.assertRaises(ValueError):
            ResidueTest('t', 'protein', 'atom_count', 6, relation='=<')

    def test_supported_relations_reach_the_emitted_condition(self):
        for rel in ('==', '!=', '<', '<=', '>', '>='):
            r = _Recorder()
            ResidueTest('t', 'protein', 'atom_count', 6, relation=rel).write(r)
            self.assertIn(f'if {{$count {rel} 6}}', r.text)

    def test_write_never_silently_emits_nothing(self):
        # defense in depth: even an object whose measure was mutated past the constructor guard
        # must raise rather than emit an empty check
        t = ResidueTest('t', 'protein', 'atom_count', 6)
        t.measure = 'nonsense'
        r = _Recorder()
        with self.assertRaises(NotImplementedError):
            t.write(r)


class TestConnectionTestRejectsUnsupportedTypes(unittest.TestCase):
    def test_supported_types_emit_a_check(self):
        for ctype in ('interresidue', 'disulfide', 'glycosylation'):
            r = _Recorder()
            ConnectionTest('t', 'protein', ctype, 2).write(r)
            self.assertIn('vmdcon "PASS', r.text)

    def test_unsupported_type_is_rejected_at_construction(self):
        with self.assertRaises(NotImplementedError):
            ConnectionTest('t', 'protein', 'interesidue', 2)

    def test_write_never_silently_emits_nothing(self):
        # the fallthrough was written `case '_':`, which matches the literal string '_' and so
        # matched nothing at all -- an unsupported type fell out of the match with no check and
        # no message
        t = ConnectionTest('t', 'protein', 'disulfide', 2)
        t.connection_type = 'nonsense'
        r = _Recorder()
        with self.assertRaises(NotImplementedError):
            t.write(r)

    def test_the_literal_underscore_is_not_a_connection_type(self):
        with self.assertRaises(NotImplementedError):
            ConnectionTest('t', 'protein', '_', 2)


class TestAttributeTestStillEmits(unittest.TestCase):
    def test_emits_pass_and_fail_branches(self):
        r = _Recorder()
        AttributeTest('t', 'protein', 'resname', 'CYS', 6).write(r)
        self.assertIn('vmdcon "PASS', r.text)
        self.assertIn('vmdcon "FAIL', r.text)


class TestSchemaEnforcesTheSameSpecsAsTheCode(unittest.TestCase):
    """The schema is the first gate and the class is the second; they have to agree, and the
    schema has to actually be *walked*.

    ``validate`` was declared ``type: list`` while its payload is a mapping.  ycleptic's list
    walker ignores a list-typed element that is not a scalar or dict, so the whole ``validate``
    subtree was skipped: every ``choices:`` under it was inert, and a bad ``measure:`` reached
    the task untouched.  ``yclept check-spec`` passes either way -- even as of 2.4.3, whose new
    structure check catches a *misplaced* attribute but not a *mistyped* one -- so only a walk
    catches this.  :class:`TestCheckSpecCatchesAMisplacedAttribute` below covers the other half.
    """

    @staticmethod
    def _base():
        import yaml
        from importlib.resources import files as pkg_files
        return yaml.safe_load(open(str(pkg_files('pestifer.schema').joinpath('base.yaml'))))

    @staticmethod
    def _node(node, name):
        for a in node.get('attributes', []) or []:
            if a.get('name') == name:
                return a
            found = TestSchemaEnforcesTheSameSpecsAsTheCode._node(a, name)
            if found:
                return found
        return None

    def _residue_test_node(self):
        validate = self._node(self._node(self._base(), 'tasks'), 'validate')
        self.assertIsNotNone(validate)
        # a mapping payload must be declared a dict, or its subtree is never walked
        self.assertEqual(validate['type'], 'dict')
        return self._node(validate, 'residue_test')

    def test_schema_choices_match_the_class_guards(self):
        node = self._residue_test_node()
        self.assertEqual(set(self._node(node, 'measure')['choices']), ResidueTest.measure_supported)
        self.assertEqual(set(self._node(node, 'relation')['choices']), ResidueTest.relation_supported)

    def test_connection_type_choices_match_the_class_guard(self):
        validate = self._node(self._node(self._base(), 'tasks'), 'validate')
        node = self._node(validate, 'connection_test')
        self.assertEqual(set(self._node(node, 'connection_type')['choices']),
                         ConnectionTest.connection_type_supported)

    def _walk(self, spec):
        from ycleptic.walkers import dwalk
        dwalk(self._base(), {'tasks': [{'validate': {'tests': [{'residue_test': spec}]}}]})

    def test_a_good_spec_walks_clean(self):
        self._walk({'name': 't', 'selection': 'protein', 'measure': 'residue_count', 'value': 6})

    def test_a_mistyped_measure_is_rejected_at_config_load(self):
        with self.assertRaises(Exception) as cm:
            self._walk({'name': 't', 'selection': 'protein', 'measure': 'residue_cont', 'value': 6})
        self.assertIn('residue_cont', str(cm.exception))

    def test_a_mistyped_relation_is_rejected_at_config_load(self):
        with self.assertRaises(Exception):
            self._walk({'name': 't', 'selection': 'protein', 'measure': 'atom_count',
                        'value': 6, 'relation': '=<'})


if __name__ == '__main__':
    unittest.main()


class TestCheckSpecCatchesAMisplacedAttribute(unittest.TestCase):
    """``yclept check-spec`` must run on the real schema, and must be able to fail.

    ``check-spec`` used to validate *vocabulary* only -- that keys and type names are recognized
    -- and passed, exit 0, on a ``base.yaml`` in which an attribute indented one level too deep
    had been swallowed into its neighbour's ``default:`` list.  The attribute then declares
    nothing: that is how ``charmmff.custom.prm`` silently ceased to exist here on 2026-09-30,
    taking the conflict-resolution parameter file out of every build with it.  ycleptic 2.4.3
    added the structure check (cameronabrams/ycleptic#7), which is why ``pyproject.toml`` floors
    there.

    Two tests, and the second is the one that earns its place.  Asserting only that the real
    schema is clean would pass against *any* ycleptic, including one with no structure check at
    all -- a floor is satisfied the moment the installed version is new enough, so the version
    number cannot tell you the code is there.  Re-introducing the slip and requiring a complaint
    is what pins that the installed ``check-spec`` actually performs this check.

    Running it here also makes it a gate.  It had been a manual step in
    ``docs/source/contributing.rst`` -- something a contributor was asked to remember -- which is
    not a guard.
    """

    @staticmethod
    def _base():
        import yaml
        from importlib.resources import files as pkg_files
        return yaml.safe_load(open(str(pkg_files('pestifer.schema').joinpath('base.yaml'))))

    def test_the_shipped_schema_is_clean(self):
        from ycleptic.speccheck import check_base_spec
        problems = check_base_spec(self._base())
        self.assertEqual(problems, [], 'pestifer/schema/base.yaml declares something ycleptic '
                                       'ignores:\n  ' + '\n  '.join(problems))

    def test_the_check_can_fail(self):
        """NEGATIVE CONTROL: swallow an attribute the way a real indentation slip does."""
        from ycleptic.speccheck import check_base_spec
        import copy

        base = self._base()

        def find(node, name):
            if isinstance(node, dict):
                if node.get('name') == name:
                    return node
                for v in node.values():
                    hit = find(v, name)
                    if hit:
                        return hit
            elif isinstance(node, list):
                for v in node:
                    hit = find(v, name)
                    if hit:
                        return hit
            return None

        custom = find(base, 'custom')
        self.assertIsNotNone(custom, 'POSITIVE CONTROL: no `custom` node to misplace')
        attrs = custom['attributes']
        prm = next(a for a in attrs if a['name'] == 'prm')
        donor = next(a for a in attrs if a['name'] == 'str')
        # exactly what four spaces too much indentation produces: the attribute leaves the
        # `attributes:` list and lands among its neighbour's default values
        attrs.remove(prm)
        donor = copy.deepcopy(donor)
        donor['default'] = list(donor.get('default') or []) + [prm]
        attrs[attrs.index(next(a for a in attrs if a['name'] == 'str'))] = donor

        problems = check_base_spec(base)
        self.assertTrue(problems, 'the installed ycleptic did not catch a swallowed attribute; '
                                  'check-spec has no structure check, so the >=2.4.3 floor is '
                                  'not doing what pyproject.toml says it does')
        self.assertTrue(any('prm' in p for p in problems),
                        f'complained, but not about the misplaced attribute: {problems}')
