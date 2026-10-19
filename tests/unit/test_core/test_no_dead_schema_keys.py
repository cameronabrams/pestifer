# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Every leaf key in ``base.yaml`` is named somewhere in the code that must consume it.

A config key the code ignores is worse than one it rejects, and **nothing else pestifer
has will find one**: ``yclept check-spec`` verifies that declarations are well-formed,
not that anything reads them; every example still parses, because no example uses the
key; and the config reference happily generates a page describing behavior that does not
exist.  Three have been found by hand --

* ``tasks.make_membrane_system.bilayer.z_dist``, accepted, validated, then dropped by
  ``bilayer_embed.tcl``'s arg loop without a word;
* ``charmmff.user_custom.{rtf,str,prm}``, a whole documented extension point read by no
  scripter, so naming your own stream there validated, built, and loaded nothing;
* ``...constraints.consfile`` in seven places, where the code reads ``consref`` -- the
  documented name did nothing and the working name could not be set at all.

-- and each was found only because somebody happened to look.  This is that look, run
every time.

**What it can and cannot prove.**  It is a *name* search, so it proves only that the key
is mentioned; it cannot prove the mention is a read.  It is also blind to a key assembled
at run time, which is why :data:`DYNAMIC_KEYS` exists.  Adding to that allowlist is
cheap, so keep it honest: each entry names the code that builds the key, and
``test_every_allowlisted_key_really_is_dynamic`` requires that code to still be there.
"""
import pathlib
import unittest

import yaml

_ROOT = pathlib.Path(__file__).resolve().parents[3]
_SPEC = _ROOT / 'pestifer' / 'schema' / 'base.yaml'
_CODE = _ROOT / 'pestifer'

#: Leaf keys no literal search can see, each with the code that composes the name.
DYNAMIC_KEYS = {
    # MembraneEquilibrateTask._params_for does specs.get(f'{prefix}_{key}') with
    # prefix='area', so area_<key> never appears as a literal.
    'tasks.membrane_equilibrate.area_precision_p':          'tasks/membrane_equilibrate.py:_params_for',
    'tasks.membrane_equilibrate.area_drift_conf':           'tasks/membrane_equilibrate.py:_params_for',
    'tasks.membrane_equilibrate.area_burn_in':              'tasks/membrane_equilibrate.py:_params_for',
    'tasks.membrane_equilibrate.area_window_frac':          'tasks/membrane_equilibrate.py:_params_for',
    'tasks.membrane_equilibrate.area_autocorr_reliability': 'tasks/membrane_equilibrate.py:_params_for',
}

#: The f-string that builds the allowlisted names.  If this stops being how they are
#: read, the allowlist is stale and the keys need re-checking by hand.
_DYNAMIC_EVIDENCE = (pathlib.Path('tasks/membrane_equilibrate.py'),
                     "specs.get(f'{prefix}_{key}')")


def _leaf_paths(node, path=()):
    if isinstance(node, dict) and 'name' in node:
        path = path + (node['name'],)
        kids = node.get('attributes') or []
        if not kids:
            yield path
        for k in kids:
            yield from _leaf_paths(k, path)
    elif isinstance(node, list):
        for k in node:
            yield from _leaf_paths(k, path)
    elif isinstance(node, dict):
        for v in node.values():
            yield from _leaf_paths(v, path)


class TestNoDeadSchemaKeys(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.leaves = sorted(set(_leaf_paths(yaml.safe_load(_SPEC.read_text()))))
        cls.source = '\n'.join(
            f.read_text(errors='replace')
            for pat in ('**/*.py', '**/*.tcl')
            for f in _CODE.glob(pat)
            if 'schema' not in f.parts)

    def test_the_corpus_was_actually_read(self):
        """Positive control.  An empty corpus would make every key look dead; a corpus
        that somehow contained everything would make every key look live."""
        self.assertGreater(len(self.leaves), 300, 'the schema did not parse')
        self.assertIn('fileset.resolve', self.source, 'the source corpus is not pestifer')
        self.assertNotIn('zzz_not_a_real_key_zzz', self.source,
                         'the corpus matches anything -- this test cannot fail')

    def test_every_key_is_named_in_the_code(self):
        dead = ['.'.join(p) for p in self.leaves
                if p[-1] not in self.source and '.'.join(p) not in DYNAMIC_KEYS]
        self.assertEqual(dead, [], 
            'these config keys are declared in base.yaml and named nowhere in '
            'pestifer/*.py or pestifer/*.tcl, so setting one has no effect and nothing '
            'reports that:\n  ' + '\n  '.join(dead) +
            '\nEither wire the key up, delete it, or -- if the code composes the name at '
            'run time -- add it to DYNAMIC_KEYS with a pointer to the code that does.')

    def test_every_allowlisted_key_really_is_dynamic(self):
        """An allowlist entry that outlives its reason turns this test into a rubber
        stamp, so require the composing code to still exist."""
        rel, needle = _DYNAMIC_EVIDENCE
        self.assertIn(needle, (_CODE / rel).read_text(),
                      f'{rel} no longer composes spec names with {needle!r}; the '
                      f'DYNAMIC_KEYS allowlist is stale')
        for key in DYNAMIC_KEYS:
            self.assertIn(key, {'.'.join(p) for p in self.leaves},
                          f'{key} is allowlisted but is no longer a schema key')
