# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A residue lookup matches the insertion code, because the index it returns into keys on one.

``CoordManipulator._residue_index`` assigns each residue a VMD-style index by
``(segname, resseqnum, insertion)`` and its docstring says "Verified to match VMD".  Both lookups
built on top of it -- ``residue_of_segname`` and its chain-keyed twin ``_residue_of`` -- compared
``resseqnum`` alone.  With residues 100, 100A and 100B in one segment that is ambiguous: it
resolved to whichever atom came first in atom order, and an insertion-coded residue could not be
addressed at all.

Three consumers of one idea; the index got it right and both lookups on top of it dropped the
field.  The same shape CLAUDE.md already records for ``cell_from_xsc`` and the atom-type case
comparisons, so the fix is one shared ``_find_residue`` rather than the same correction written
twice.

The insertion was discarded TWICE on the crot path: ``apply_crot`` called
``self._residue_of(chain, crot.resid1.resseqnum)``, unpacking a ``ResID`` into an int before a
lookup that would have dropped it anyway.  Passing the ``ResID`` whole is what makes the loss
visible at the call site instead of invisible inside it.

NOT a defect, and deliberately untouched: ``coordmanip.py:133``'s ``resid N`` selection term.
VMD's own ``resid`` keyword is the numeric field and ignores insertion codes (``insertion`` is a
separate keyword there), so matching on resseqnum alone is correct in the thing that mirrors VMD.
"""
import inspect
import unittest
from pathlib import Path

from pestifer.molecule.coordmanip import CoordManipulator, CoordManipulateError
from pestifer.objs.resid import ResID

#: PROA/chain A residues 100, 100A, 100B -- the classic antibody-CDR numbering, and the smallest
#: structure that distinguishes the three.  Loaded without a PSF, so pestifer reads the segname
#: from the chain; both lookups are still exercised, since they differ only in the field matched.
FIX = str(Path(__file__).parent / 'fixtures' / 'insertion_codes.pdb')


class TestInsertionCodeLookup(unittest.TestCase):

    def setUp(self):
        self.cm = CoordManipulator(None, FIX)

    def test_the_fixture_really_carries_insertion_codes(self):
        """POSITIVE CONTROL. Without this the rest could pass against a file with none."""
        codes = sorted({(a.resid.insertion or '') for a in self.cm.atoms.data})
        self.assertEqual(codes, ['', 'A', 'B'])
        self.assertEqual(sorted({a.resid.resseqnum for a in self.cm.atoms.data}), [100])
        self.assertEqual(len(set(self.cm._residue_index())), 3,
                         'the index itself must see three residues, or there is nothing to resolve')

    def test_the_three_residues_resolve_to_three_different_indices(self):
        got = [self.cm.residue_of_segname('A', r) for r in ('100', '100A', '100B')]
        self.assertEqual(len(set(got)), 3, f'100/100A/100B collapsed onto {got}')

    def test_each_index_is_the_residue_that_was_asked_for(self):
        """The sharp version: distinct indices could still be the wrong three.

        Checks the resolved index back against the insertion code of the atoms it covers, which
        is positive evidence about the specific claim rather than the absence of a collision.
        """
        ridx = self.cm._residue_index()
        for spec, want in (('100', ''), ('100A', 'A'), ('100B', 'B')):
            r = self.cm.residue_of_segname('A', spec)
            codes = {(a.resid.insertion or '') for i, a in enumerate(self.cm.atoms.data)
                     if ridx[i] == r}
            self.assertEqual(codes, {want}, f'{spec} resolved to a residue holding {codes}')

    def test_the_chain_keyed_twin_behaves_the_same(self):
        """``_residue_of`` had the identical defect and six call sites of its own."""
        got = [self.cm._residue_of('A', r) for r in ('100', '100A', '100B')]
        self.assertEqual(len(set(got)), 3, f'100/100A/100B collapsed onto {got}')
        self.assertEqual(got, [self.cm.residue_of_segname('A', r) for r in ('100', '100A', '100B')])

    def test_a_ResID_and_its_string_are_the_same_query(self):
        for spec in ('100', '100A', '100B'):
            self.assertEqual(self.cm.residue_of_segname('A', ResID(spec)),
                             self.cm.residue_of_segname('A', spec))

    def test_a_bare_int_means_NO_insertion_code(self):
        """Strictly, not "any residue with this resseqnum".

        That is the whole ambiguity: a bare 100 used to be able to return 100A. Callers that
        genuinely hold only a number still work and now get a defined answer.
        """
        self.assertEqual(self.cm.residue_of_segname('A', 100),
                         self.cm.residue_of_segname('A', ResID('100')))

    def test_a_residue_that_is_not_there_raises(self):
        for missing in ('100C', '101', 999):
            with self.assertRaises(CoordManipulateError):
                self.cm.residue_of_segname('A', missing)

    def test_the_old_resseqnum_only_lookup_WOULD_fail_this(self):
        """NEGATIVE CONTROL at the level the fix lives.

        Restores the pre-fix comparison and requires the three residues to collapse. A guard
        that cannot fail is the same bug it is guarding against.
        """
        ridx = self.cm._residue_index()

        def old_lookup(segname, resid):
            want = resid if isinstance(resid, int) else ResID(resid).resseqnum
            for i, a in enumerate(self.cm.atoms.data):
                if a.segname == segname and a.resid.resseqnum == want:
                    return int(ridx[i])
            return None

        got = [old_lookup('A', r) for r in ('100', '100A', '100B')]
        self.assertEqual(len(set(got)), 1,
                         'the pre-fix lookup is supposed to collapse all three onto one index')


class TestTheCrotCallSitePassesTheWholeResID(unittest.TestCase):
    """The field was dropped twice: once by the caller, once by the lookup. Fixing only the
    lookup leaves `crot.resid1.resseqnum` throwing the insertion away before it ever arrives."""

    def setUp(self):
        self.src = inspect.getsource(CoordManipulator.apply_crot)

    def test_no_call_site_unpacks_a_resid_into_an_int(self):
        for k in (1, 2, 3):
            self.assertNotIn(f'crot.resid{k}.resseqnum', self.src)

    def test_the_resid_is_passed_whole(self):
        self.assertIn('self._residue_of(chain, crot.resid1)', self.src)


if __name__ == '__main__':
    unittest.main()
