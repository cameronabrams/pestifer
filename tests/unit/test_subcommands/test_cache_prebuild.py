# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""``pestifer cache prebuild``: generate one lipid conformer set ahead of the build that needs it.

A build generates a missing ``(lipid, phase)`` conformer set itself, which is usually what you
want.  This command exists for when that is inconvenient: generating inside a batch job spends
allocation time on a single-molecule vacuum run, and it was the only workaround while such a run
could not be launched under SLURM at all (reported 2026-09-21, fixed in the scripter).

What these tests pin is that the prebuilt entry is the entry the build would have made -- same
phase-qualified name, same sampler.  An entry built with a different sampler would be a cache hit
that quietly changes the science.
"""
import argparse
import unittest
from unittest import mock

from pestifer.subcommands.cache import CacheSubcommand


def _args(**kw):
    kw.setdefault('action', 'prebuild')
    kw.setdefault('resname', 'PSM')
    kw.setdefault('phase', 'Ld')
    kw.setdefault('charmmff_release', '')
    return argparse.Namespace(**kw)


class _CC(dict):
    """Stands in for CHARMMFFContent: membership is 'is this RESI defined'."""

    def __init__(self, resnames=('PSM',)):
        super().__init__()
        self.release_str = 'feb26'
        self._resnames = set(resnames)

    def __contains__(self, item):
        return item in self._resnames

    def provision(self):
        pass


class TestCachePrebuild(unittest.TestCase):

    def _run(self, args):
        calls = {}

        def _ensure(resname, CC, *, phase=None, sampler=None, **kw):
            calls.update(resname=resname, phase=phase, sampler=sampler)
            from pathlib import Path
            return Path('/cache/lipid')

        with mock.patch('pestifer.charmmff.autocache.ensure_lipid_conformer', _ensure), \
             mock.patch('pestifer.core.resourcemanager.ResourceManager') as RM:
            RM.return_value.charmmff_content = _CC()
            CacheSubcommand.func(args)
        return calls

    def test_a_phased_set_uses_the_mc_sampler_like_a_build_does(self):
        # Bilayer picks 'mc' for Ld/Lo; a prebuilt entry from the legacy 'md' sampler would be
        # extended rods where the build expects a melted ensemble
        self.assertEqual(self._run(_args(phase='Lo')),
                         {'resname': 'PSM', 'phase': 'Lo', 'sampler': 'mc'})

    def test_the_default_phase_is_the_fluid_one(self):
        self.assertEqual(self._run(_args())['phase'], 'Ld')

    def test_a_missing_resname_is_refused_with_the_command_to_type(self):
        with self.assertRaises(ValueError) as e:
            self._run(_args(resname=''))
        self.assertIn('--resname', str(e.exception))

    def test_a_resname_absent_from_the_force_field_is_refused(self):
        with mock.patch('pestifer.core.resourcemanager.ResourceManager') as RM:
            RM.return_value.charmmff_content = _CC(resnames=('POPC',))
            with self.assertRaises(ValueError) as e:
                CacheSubcommand.func(_args(resname='NOPE'))
        self.assertIn('not defined in the CHARMM force field', str(e.exception))
