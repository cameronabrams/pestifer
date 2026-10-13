# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""The embed clash test runs on heavy atoms, because 2.4 A between ALL atoms is not a clash.

``bilayer_embed.tcl`` deletes any lipid with an atom within 2.4 A of the protein, and one
contacting atom removes the whole residue. Both selections were ``atomselect ... "all"``, so
hydrogens counted on both sides — and that makes 2.4 A a *contact* criterion, not an overlap one.
A C-H bond is ~1.1 A, so two carbons at ordinary van der Waals separation present hydrogens well
inside 2.4 A.

Measured here on a real pestifer-built membrane protein system (ex17 prefill, 114,738 atoms):
across 22 perimeter lipids the **median closest approach is 1.4 A larger on heavy atoms than on
all atoms**. So the old 2.4 A all-atom cut corresponds to roughly 3.8 A between heavy atoms —
squarely vdW contact for C...C (3.4-4.0 A), which is packing, not a clash.

Reported externally with three independent reproductions: 36 lipids removed where ~21 overlap
(PS5 with its 711-716 linker), 28 vs ~21 (PS5, TM only), and 33 vs ~25 (M2 — a different
protein, a different composition, a different pestifer version). The excess tracks TM+linker
*perimeter* rather than composition, which is the signature of a contact test.

The constant is unchanged. What changed is which atoms it applies to: 2.4 A between heavy atoms
is a genuine overlap.

There is no unit test that can run VMD here, so this pins the script's source. That is weaker
than exercising it — and it is exactly the gap that let a shipped Tcl script carry this for as
long as it did.
"""
import re
import unittest
from importlib.resources import files as pkg_files

SRC = (pkg_files('pestifer.resources.tcl.scripts') / 'bilayer_embed.tcl').read_text()

#: every clash test in the script
CONTACT_LINES = [ln.strip() for ln in SRC.splitlines()
                 if 'measure contacts' in ln and not ln.strip().startswith('#')]


class TestBothClashTestsUseHeavyAtoms(unittest.TestCase):

    def test_there_are_exactly_two_clash_tests(self):
        """One for lipids, one for water. If a third appears, it needs the same decision."""
        self.assertEqual(len(CONTACT_LINES), 2, CONTACT_LINES)

    def test_neither_uses_an_all_atom_selection(self):
        """NEGATIVE CONTROL: the exact pre-fix lines must not come back."""
        for ln in CONTACT_LINES:
            self.assertNotIn('$bilayer_sel', ln, f'all-atom lipid clash test is back: {ln}')
            self.assertNotIn('$water_sel', ln, f'all-atom water clash test is back: {ln}')

    def test_both_test_heavy_against_heavy(self):
        for ln in CONTACT_LINES:
            sels = re.findall(r'\$(\w+)', ln)
            self.assertEqual(len(sels), 2, ln)
            for sel in sels:
                self.assertIn('heavy', sel,
                              f'{sel} in "{ln}" is not a heavy-atom selection')

    def test_every_heavy_selection_actually_excludes_hydrogens(self):
        """A selection NAMED heavy that is not would pass the test above and fail in the field."""
        for name in set(re.findall(r'set (\w*heavy\w*) \[atomselect', SRC)):
            decl = re.search(rf'set {name} \[atomselect[^\]]*\]', SRC).group(0)
            self.assertIn('noh', decl, f'{name} does not exclude hydrogens: {decl}')

    def test_an_empty_heavy_selection_is_refused_rather_than_silently_deleting_nothing(self):
        """VMD answers a bad selection with zero atoms, not an error.

        Without this guard a typo would make the clash test find nothing and the embed would
        look clean while leaving every overlapping lipid in place -- a silent failure in the
        direction nobody checks.
        """
        self.assertEqual(SRC.count('refusing to embed without a working clash test')
                         + SRC.count('heavy-atom selection is empty for the water clash test'), 2)
        for guard in re.findall(r'if \{\[\$\w*heavy\w* num\] == 0', SRC):
            self.assertTrue(guard)

    def test_the_deletions_are_reported(self):
        """The embed stage used to delete lipids and say nothing, so an over-removal was
        invisible until someone counted residues in the output PSF by hand."""
        self.assertIn('deleted [array size _seen_membrane] lipid(s)', SRC)
        self.assertIn('deleted [array size _seen_water] water(s)', SRC)


class TestTheIndexMappingAssumptionIsRecorded(unittest.TestCase):
    """`measure contacts` returns GLOBAL atom indices, not selection-relative ones.

    The (segname,resid) lookup reads the parent "all" lists by that index, so narrowing the
    contact selections was only safe because of this. Verified against VMD on a real psf/pdb: a
    contact index of 13741 came back from a selection of 4870 atoms. Recorded in the script
    because the next person to touch those selections has to know it.
    """

    def test_the_script_says_so(self):
        self.assertIn('GLOBAL atom indices', SRC)


if __name__ == '__main__':
    unittest.main()
