# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
`charmmff.user_pdbcollections` must actually reach the PDB repository.

It did not, for as long as the key has existed.  Three spellings of one feature were in play:

    schema key                  charmmff.user_pdbcollections   (schema/base.yaml)
    the schema's own example    charmmff.pdbcollections
    what the code read          charmmff.pdbrepository         (resourcemanager.py)

So a user following the schema set a key nothing read; a user following the documented example
set a key that is not even schema-valid; and the only spelling that worked was undeclared, so it
could not survive validation.  Nothing failed anywhere -- the collection was silently ignored and
the shipped one used instead, which looks exactly like a correct build.

Found 2026-10-02 while looking for a way to test a regenerated conformer collection WITHOUT
installing it over the shipped one.  That is the feature's real job, so it is worth a guard.
"""
import unittest

import yaml
from importlib.resources import files as pkg_files


def _charmmff_schema_keys():
    base = yaml.safe_load(open(str(pkg_files('pestifer.schema').joinpath('base.yaml'))))
    cff = next(a for a in base['attributes'] if a['name'] == 'charmmff')
    return {a['name'] for a in cff['attributes']}


class TestUserPDBCollectionsIsWiredUp(unittest.TestCase):

    def test_the_key_is_declared_in_the_schema(self):
        self.assertIn('user_pdbcollections', _charmmff_schema_keys())

    def test_the_config_key_reaches_the_charmmff_content(self):
        """The round trip that was broken: config key in, attribute out."""
        from pestifer.core.resourcemanager import ResourceManager
        rm = ResourceManager(charmmff_config={'user_pdbcollections': ['/nonexistent/collection']})
        self.assertEqual(list(rm.charmmff_content.user_pdbcollections), ['/nonexistent/collection'])

    def test_a_key_the_schema_does_not_declare_cannot_be_the_only_way_in(self):
        """NEGATIVE CONTROL for the shape of the original bug.

        Reading a key the schema never declares is the defect, because such a key cannot reach the
        code through a validated config at all.  Whatever spelling resourcemanager honours as its
        PRIMARY source must be a declared one.
        """
        import inspect
        from pestifer.core import resourcemanager
        src = inspect.getsource(resourcemanager)
        line = next(l for l in src.splitlines() if 'user_pdbcollections=' in l and '.get(' in l)
        primary = line.split(".get('", 1)[1].split("'", 1)[0]
        self.assertIn(primary, _charmmff_schema_keys(),
                      f"resourcemanager reads '{primary}' first, which the schema does not declare")

    def test_the_schema_example_uses_the_schema_key(self):
        """The example is what a user copies; it named a third spelling."""
        base = yaml.safe_load(open(str(pkg_files('pestifer.schema').joinpath('base.yaml'))))
        cff = next(a for a in base['attributes'] if a['name'] == 'charmmff')
        node = next(a for a in cff['attributes'] if a['name'] == 'user_pdbcollections')
        example = node.get('docs', {}).get('example', {})
        self.assertIn('user_pdbcollections', example.get('charmmff', {}),
                      f'the documented example names something else: {list(example.get("charmmff", {}))}')


if __name__ == '__main__':
    unittest.main()
