# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Display folding of '<resname>__Lo' conformer-ensemble keys.

A '__Lo' key is the same residue's liquid-ordered conformer ensemble, not a distinct residue, so
the listing shows one entry per residue with the phase as a marker.  The two markers stay distinct
because a collection can hold an Lo ensemble without the base conformers -- the on-demand user
cache is exactly that -- and one marker would imply it provides a residue it does not.
"""
import unittest

from pestifer.charmmff.pdbrepository import PDBCollection


class _Coll(PDBCollection):
    """A PDBCollection with its info dict injected, bypassing tarball/directory discovery."""
    def __init__(self, info):
        self.info = info
        self.streamID = 'lipid'
        self.path_or_tarball = '<test>'
        self.registration_place = 1


class TestLoDisplayFolding(unittest.TestCase):

    def test_base_and_lo_fold_to_one_starred_entry(self):
        c = _Coll({'DMPC': {}, 'DMPC__Lo': {}, 'DOPC': {}})
        self.assertEqual(c._display_entries(), [('DMPC', '*'), ('DOPC', '')])

    def test_lo_without_base_is_marked_differently(self):
        # the on-demand cache holds generated Lo ensembles only
        c = _Coll({'DMPC__Lo': {}, 'PSM__Lo': {}})
        self.assertEqual(c._display_entries(), [('DMPC', '^'), ('PSM', '^')])

    def test_counts_report_residues_not_keys(self):
        c = _Coll({'DMPC': {}, 'DMPC__Lo': {}, 'DOPC': {}})
        head = c.show().splitlines()[0]
        self.assertIn('contains 2 residues', head)
        self.assertIn('1 with an Lo conformer ensemble', head)

    def test_lo_only_collection_says_so(self):
        c = _Coll({'DMPC__Lo': {}})
        head = c.show().splitlines()[0]
        self.assertIn('contains 1 residues', head)
        self.assertIn('Lo ensembles only', head)

    def test_legend_only_appears_when_a_marker_is_used(self):
        self.assertNotIn('liquid-ordered', _Coll({'DOPC': {}}).show())
        self.assertIn('liquid-ordered', _Coll({'DOPC': {}, 'DOPC__Lo': {}}).show())


if __name__ == '__main__':
    unittest.main()
