# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Which collection wins when two of them supply the same residue.

``checkout`` walks ``registration_order`` in REVERSE, so the last collection registered wins.
The on-demand generation cache was auto-registered last, which made it outrank both the shipped
repository and an explicit ``user_pdbcollections`` override -- the exact inverse of what
:mod:`pestifer.charmmff.autocache` documents it to be ("a residue that ... has **no** entry in
the shipped PDB repository").

The cost was not theoretical.  A ``POPC__Lo`` generated in August 2026 -- before the shipped
collection carried readable ``__Lo`` entries -- went on overriding the shipped conformers under
3.25.1 on every machine that had ever built with a declared leaflet phase, so the conformer work
in 3.25.0 and 3.25.1 did not reach those machines at all.  Nothing said so: the one line marking
a shadowed collection claimed the opposite ("already registered; will not add again") while the
collection was in fact added, under a suffixed key, at HIGHER precedence.

Measured on panacea 2026-10-05, building ex17's patchA from the 3.25.1 release with and without
a cached ``POPC__Lo`` present: 42.77 A vs 38.03 A of built phosphate-to-phosphate thickness, from
one cached entry, with no diagnostic distinguishing the two runs.
"""
import unittest

from unittest.mock import patch

from pestifer.charmmff.pdbrepository import PDBCollection, PDBRepository


class _Coll(PDBCollection):
    """A PDBCollection with its info dict injected, bypassing tarball/directory discovery."""

    def __init__(self, info, streamID='lipid', tag='<test>'):
        self.info = info
        self.streamID = streamID
        self.path_or_tarball = tag
        self.registration_place = 0

    def checkout(self, name):
        return self.path_or_tarball if name in self.info else None


def _repo():
    r = PDBRepository.__new__(PDBRepository)
    r.collections = {}
    r.registration_order = []
    return r


class TestPrecedenceOrder(unittest.TestCase):

    def test_a_fallback_collection_loses_to_one_already_registered(self):
        r = _repo()
        r.add_collection(_Coll({'POPC__Lo': {}}, tag='shipped'), collection_key='lipid')
        r.add_collection(_Coll({'POPC__Lo': {}}, tag='cache'), collection_key='lipid',
                         fallback=True)
        self.assertEqual(r.checkout('POPC__Lo'), 'shipped')

    def test_a_fallback_collection_loses_to_one_registered_after_it(self):
        """Order of the calls must not decide it -- only the flag.

        The cache is registered after the shipped repository in one path and could be registered
        before it in another; a fix that depended on call order would hold in one and not the
        other.
        """
        r = _repo()
        r.add_collection(_Coll({'POPC__Lo': {}}, tag='cache'), collection_key='lipid',
                         fallback=True)
        r.add_collection(_Coll({'POPC__Lo': {}}, tag='shipped'), collection_key='lipid')
        self.assertEqual(r.checkout('POPC__Lo'), 'shipped')

    def test_without_the_flag_the_later_registration_still_wins(self):
        """POSITIVE CONTROL for the two above.

        An explicit ``user_pdbcollections`` override is registered after the shipped repository
        and MUST keep beating it -- that is the documented way to test a collection without
        installing it.  If this went the other way, the tests above would pass for the wrong
        reason (nothing ever overriding anything).
        """
        r = _repo()
        r.add_collection(_Coll({'POPC__Lo': {}}, tag='shipped'), collection_key='lipid')
        r.add_collection(_Coll({'POPC__Lo': {}}, tag='user'), collection_key='lipid')
        self.assertEqual(r.checkout('POPC__Lo'), 'user')

    def test_all_three_ranks_at_once(self):
        """explicit user override > shipped > generation cache, in one repository."""
        r = _repo()
        r.add_collection(_Coll({'POPC__Lo': {}}, tag='shipped'), collection_key='lipid')
        r.add_collection(_Coll({'POPC__Lo': {}}, tag='user'), collection_key='lipid')
        r.add_collection(_Coll({'POPC__Lo': {}}, tag='cache'), collection_key='lipid',
                         fallback=True)
        self.assertEqual(r.checkout('POPC__Lo'), 'user')
        del r.collections[r.registration_order.pop()]          # drop the user override
        self.assertEqual(r.checkout('POPC__Lo'), 'shipped')

    def test_a_cached_residue_the_release_lacks_is_still_found(self):
        """What must NOT change: the cache's whole purpose is filling gaps."""
        r = _repo()
        r.add_collection(_Coll({'POPC__Lo': {}}, tag='shipped'), collection_key='lipid')
        r.add_collection(_Coll({'WEIRD': {}}, tag='cache'), collection_key='lipid',
                         fallback=True)
        self.assertEqual(r.checkout('WEIRD'), 'cache')
        self.assertIn('WEIRD', r)

    def test_registration_place_is_renumbered_after_a_front_insert(self):
        r = _repo()
        r.add_collection(_Coll({'A': {}}, tag='shipped'), collection_key='lipid')
        r.add_collection(_Coll({'B': {}}, tag='cache'), collection_key='lipid', fallback=True)
        places = [r.collections[k].registration_place for k in r.registration_order]
        self.assertEqual(places, [1, 2], 'display ranks went stale after the insert')


class TestStaleCacheIsReported(unittest.TestCase):

    def test_overlap_with_a_higher_rank_is_listed(self):
        r = _repo()
        r.add_collection(_Coll({'POPC__Lo': {}, 'PSM__Lo': {}}, tag='shipped'),
                         collection_key='lipid')
        r.add_collection(_Coll({'POPC__Lo': {}, 'WEIRD': {}}, tag='cache'),
                         collection_key='lipid', fallback=True)
        self.assertEqual(r.shadowed_by_higher_precedence('lipid_1'), ['POPC__Lo'])

    def test_a_gap_filling_cache_reports_nothing(self):
        r = _repo()
        r.add_collection(_Coll({'PSM__Lo': {}}, tag='shipped'), collection_key='lipid')
        r.add_collection(_Coll({'WEIRD': {}}, tag='cache'), collection_key='lipid',
                         fallback=True)
        self.assertEqual(r.shadowed_by_higher_precedence('lipid_1'), [])

    def test_the_shipped_collection_is_not_reported_as_shadowed_by_the_cache(self):
        """Direction matters: the question is what the CACHE supplies in vain, not the reverse."""
        r = _repo()
        r.add_collection(_Coll({'POPC__Lo': {}}, tag='shipped'), collection_key='lipid')
        r.add_collection(_Coll({'POPC__Lo': {}}, tag='cache'), collection_key='lipid',
                         fallback=True)
        self.assertEqual(r.shadowed_by_higher_precedence('lipid'), [])


class TestTheCallSiteUsesIt(unittest.TestCase):
    """The helper being right is not evidence the caller calls it right.

    `provision_pdbrepository` is where the cache is auto-registered, and a fix that stopped at
    `add_collection` would leave the defect exactly where it was.  Reverting `fallback=True` at
    that call site must turn this red.
    """

    def _provision(self, cache_dirs):
        from pestifer.charmmff.charmmffcontent import CHARMMFFContent
        c = CHARMMFFContent.__new__(CHARMMFFContent)
        c.charmmff_path = '/nonexistent/feb26'
        c.user_pdbcollections = []
        calls = []

        def fake_add_resource(path, *args, fallback=False, **kwargs):
            calls.append((path, fallback))

        repo = _repo()
        repo.add_resource = fake_add_resource
        with patch('pestifer.charmmff.charmmffcontent.PDBRepository', return_value=repo), \
             patch('pestifer.charmmff.autocache.cached_collection_dirs',
                   return_value=cache_dirs):
            CHARMMFFContent.provision_pdbrepository(c)
        return calls

    def test_cached_collections_are_registered_as_fallbacks(self):
        calls = self._provision(['/home/u/.pestifer/pdbrepository/feb26/lipid'])
        self.assertEqual(calls, [('/home/u/.pestifer/pdbrepository/feb26/lipid', True)])

    def test_an_explicit_user_collection_is_not_a_fallback(self):
        """The other half of the contract, at the same call site."""
        from pestifer.charmmff.charmmffcontent import CHARMMFFContent
        c = CHARMMFFContent.__new__(CHARMMFFContent)
        c.charmmff_path = '/nonexistent/feb26'
        c.user_pdbcollections = ['/some/user/collection']
        calls = []
        repo = _repo()
        repo.add_resource = lambda path, *a, fallback=False, **k: calls.append((path, fallback))
        with patch('pestifer.charmmff.charmmffcontent.PDBRepository', return_value=repo), \
             patch('pestifer.charmmff.autocache.cached_collection_dirs', return_value=[]):
            CHARMMFFContent.provision_pdbrepository(c)
        self.assertEqual(calls, [('/some/user/collection', False)])


if __name__ == '__main__':
    unittest.main()
