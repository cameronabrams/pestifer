# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""The release defines some terms twice; this pins which value pestifer uses.

CHARMM resolves a duplicated term by load order -- upstream says so itself in
``toppar/toppar_all.history`` at entry ``2026_1_24``: "In all cases the final version of a
parameter is the one used for the calculations."  That makes the answer depend on which file
happened to be read last rather than on any judgment about which value is right.

``toppar_pestifer_conflict_resolutions.prm`` is where that judgment is recorded.  It is listed
in ``charmmff.custom.prm``, which :meth:`NAMDScripter.fetch_standard_charmm_parameters` appends
after the standard set, so an entry in it wins whatever order the conflicting sources were read
in.  These tests pin both halves: the file is still registered, and its value still wins.
"""
import copy
import os
import unittest

import yaml

import pestifer
from pestifer.charmmff.charmmffprm import CharmmParamFile

SCHEMA = os.path.join(os.path.dirname(pestifer.__file__), 'schema', 'base.yaml')
CUSTOM_DIR = os.path.join(os.path.dirname(pestifer.__file__), 'resources', 'charmmff', 'custom')
OVERRIDE = 'toppar_pestifer_conflict_resolutions.prm'

#: The one substantive conflict in the shipped default set: a nitro and a primary amine on one
#: aromatic ring.  par_all36_cgenff.prm fits it at 1.25; toppar_all36_carb_imlab.str appends 3.1
#: derived by analogy from the DIAMINE quartet, under a `!DNAP` label meant for one ligand.
CONFLICTED = ('NG2O1', 'CG2R61', 'CG2R61', 'NG2S3')
CHOSEN_KCHI = 1.25          # CGenFF's fitted value; chosen by CFA 2026-09-28
IMLAB_KCHI = 3.1            # what a release-only build got, by load order alone


def _custom_prm_defaults():
    """The ``charmmff.custom.prm`` default list, walked out of the schema."""
    schema = yaml.safe_load(open(SCHEMA))

    def find(node, name):
        if isinstance(node, dict):
            if node.get('name') == name:
                return node
            for v in node.values():
                r = find(v, name)
                if r:
                    return r
        elif isinstance(node, list):
            for v in node:
                r = find(v, name)
                if r:
                    return r
        return None

    custom = find(schema, 'custom')
    assert custom is not None, 'POSITIVE CONTROL: no `custom` node in the schema'
    for a in custom['attributes']:
        assert isinstance(a, dict), f'schema malformed near custom: {a!r}'
        if a['name'] == 'prm':
            return a.get('default') or []
    raise AssertionError('POSITIVE CONTROL: `custom` has no `prm` attribute')


class TestConflictResolutionsAreLoaded(unittest.TestCase):

    def test_the_override_file_is_registered_as_a_custom_prm(self):
        """Registration is the whole mechanism -- the file wins only because it is loaded last.

        An override that ships but is not listed changes nothing, and nothing else would fail.
        """
        self.assertIn(OVERRIDE, _custom_prm_defaults())

    def test_the_override_file_ships(self):
        self.assertTrue(os.path.exists(os.path.join(CUSTOM_DIR, OVERRIDE)))

    def test_it_is_listed_after_the_degenerate_torsion_fills(self):
        """Ordering among custom files is not load-bearing today, but the fills file states a
        different contract (degenerate torsions only), so the two must stay distinct entries."""
        defaults = _custom_prm_defaults()
        self.assertIn('toppar_pestifer_dihedral_fills.prm', defaults)
        self.assertNotEqual(OVERRIDE, 'toppar_pestifer_dihedral_fills.prm')


class TestChosenValueWins(unittest.TestCase):
    """The override must decide the term regardless of the order the release files were read."""

    def setUp(self):
        self.override = CharmmParamFile.from_file(os.path.join(CUSTOM_DIR, OVERRIDE))

    def _kchi(self, *files):
        combined = CharmmParamFile()
        for f in files:
            combined.merge(f)
        by_key = {CharmmParamFile._dihedral_key(d): d for d in combined.dihedrals}
        key = (CONFLICTED, 2)
        self.assertIn(key, by_key, 'POSITIVE CONTROL: the quartet is absent, so order proves nothing')
        return by_key[key].Kchi

    def _competing_definition(self):
        """A stand-in for the release's other definition of the same quartet.

        deepcopy, not `list(...)`: a shallow copy shares the dihedral OBJECTS, so setting Kchi
        on it also rewrites the override's own value and both sides end up equal -- which reads
        as the override losing.  The first version of this fixture did exactly that.
        """
        other = CharmmParamFile()
        other.dihedrals = copy.deepcopy(list(self.override.dihedrals))
        for d in other.dihedrals:
            d.Kchi = IMLAB_KCHI
        return other

    def test_the_override_carries_the_chosen_value(self):
        self.assertEqual(self._kchi(self.override), CHOSEN_KCHI)

    def test_it_beats_a_conflicting_definition_read_before_it(self):
        imlab_like = self._competing_definition()
        self.assertEqual(self._kchi(imlab_like, self.override), CHOSEN_KCHI)

    def test_the_fixture_can_actually_fail(self):
        """Guards the test above: with the override read FIRST, the other value must win.

        If it did not, `_kchi` would be insensitive to order and the test proves nothing.
        """
        imlab_like = self._competing_definition()
        self.assertEqual(self._kchi(self.override, imlab_like), IMLAB_KCHI)
