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

**Superseded 2026-10-13 by a fractional-vdW criterion, because the flat heavy-atom cut was wrong
in the other direction.** Keeping the 2.4 A constant and changing the atoms it applied to
LOOSENED the test by ~1.4 A: the old all-atom cut was effectively a ~3.8 A heavy-atom cut (the
measured offset above), so 2.4 A heavy retained everything in the 2.4-3.4 A band, which is
genuine steric overlap. An observed build then hit 148 constraint failures across protein and
lipid at the first post-embed dynamics (job 26469605 rep-03).

No single distance can be right, which is why both flat values failed. A lipid phosphate
hydrogen bonded to a protein nitrogen sits at ~2.8 A and must be KEPT; two carbons at 2.8 A are
badly overlapped and must go. Scaling the summed van der Waals radii separates them:

    C...C  1.70+1.70 = 3.40  x0.85 = 2.89   a 2.8 A contact is removed  (clash)
    O...N  1.52+1.55 = 3.07  x0.85 = 2.61   a 2.8 A contact is kept     (H-bond)
    O...O  1.52+1.52 = 3.04  x0.85 = 2.58   a 2.8 A contact is kept     (H-bond, water)

0.85 is mid-window and not a taste: removing a 2.8 A C...C clash needs f > 0.824 and keeping a
2.8 A hydrogen bond needs f < 0.912, so either endpoint is one rounding from breaking the
property it was picked for. Verified against known geometry on synthetic pairs: C-C at 2.50
flagged, C-C at 2.90 kept, O-N at 2.80 kept, O-O at 2.30 flagged.

Two things the implementation must NOT do, both verified on a built system: VMD's own ``radius``
is a display value, not van der Waals -- 1.5 for carbon against 1.70, 1.3 for oxygen against
1.52, systematically ~0.2 A small. And the element field is ``X`` for every atom, so it cannot be
read either. Radii key off the atom name's first character, which across a built
membrane-protein system is only C, N, O, P.

There is no unit test that can run VMD here, so this pins the script's source. That is weaker
than exercising it — and it is exactly the gap that let a shipped Tcl script carry this for as
long as it did.
"""
import re
import unittest
from importlib.resources import files as pkg_files

SRC = (pkg_files('pestifer.resources.tcl.scripts') / 'bilayer_embed.tcl').read_text()

#: every clash test in the script.  Matched on `[measure contacts`, with the bracket, so the
#: phrase quoted inside the positive-control error message is not counted as a second search.
CONTACT_LINES = [ln.strip() for ln in SRC.splitlines()
                 if '[measure contacts' in ln and not ln.strip().startswith('#')]


class TestBothClashTestsUseHeavyAtoms(unittest.TestCase):

    def test_there_is_ONE_shared_overlap_test(self):
        """Both sites route through `overlapping_indices`, so the rule cannot drift apart
        between them -- which it did when they were two flat constants."""
        self.assertEqual(len(CONTACT_LINES), 1,
                         f'expected the single search inside overlapping_indices: {CONTACT_LINES}')
        self.assertIn('_OVERLAP_SEARCH', CONTACT_LINES[0])
        self.assertEqual(SRC.count('overlapping_indices $'), 2,
                         'both the lipid and the water site must use the shared test')

    def test_neither_uses_an_all_atom_selection(self):
        """NEGATIVE CONTROL: the original pre-fix lines must not come back."""
        for ln in CONTACT_LINES + re.findall(r'overlapping_indices [^\n]*', SRC):
            self.assertNotIn('$bilayer_sel', ln, f'all-atom lipid clash test is back: {ln}')
            self.assertNotIn('$water_sel', ln, f'all-atom water clash test is back: {ln}')

    def test_both_sites_pass_heavy_atom_selections(self):
        """`overlapping_indices molid_a sel_a molid_b sel_b` -- the SELECTIONS are the 2nd and
        4th arguments.  The signature grew a second molid on 2026-10-16 because the old
        three-argument form assumed one molecule and let the lipid site pass a selection from
        another, whose indices were then read against the wrong coordinates."""
        # `overlapping_indices \$...` with the dollar excludes the `proc` definition line,
        # whose argument list is bare names
        calls = re.findall(r'overlapping_indices (\$[^\n\]]*)', SRC)
        self.assertEqual(len(calls), 2, f'expected two call sites, found {calls}')
        for ln in calls:
            args = re.findall(r'\$(\w+)', ln)
            self.assertEqual(len(args), 4, f'not a four-argument call: "{ln}"')
            for sel in (args[1], args[3]):
                self.assertIn('heavy', sel, f'{sel} in "{ln}" is not a heavy-atom selection')

    def test_each_index_list_is_read_in_its_own_molecule(self):
        """The defect that made the whole test a no-op: `measure contacts` returns two index
        lists, each in its own molecule's space, and both were read from one molecule's arrays.
        A lipid C218 and a protein OT2 3.497 A apart computed as 62.873 A, so every pair failed
        and ZERO of 465 overlapping lipids were deleted."""
        body = re.search(r'proc overlapping_indices.*?\n\}', SRC, re.S).group(0)
        for ax in 'xyz':
            self.assertIn(f'[lindex ${ax}a $i] - [lindex ${ax}b $j]', body,
                          f'{ax}: the piercer coordinate is not read from its own molecule')
        self.assertIn('names_b $j', body, 'the piercer element is not read from its own molecule')

    def test_the_arithmetic_has_a_positive_control(self):
        """Checking that the two selections are non-empty did not catch the mis-indexing: the
        selections were fine, the arithmetic between them was not.  Every pair `measure contacts
        R` returns is within R by construction, so a computed distance above R means we are not
        measuring the pair VMD gave us."""
        body = re.search(r'proc overlapping_indices.*?\n\}', SRC, re.S).group(0)
        self.assertIn('_OVERLAP_SEARCH + 1.0e-3', body, 'no bound on the computed distance')
        self.assertIn('error "bilayer_embed:', body, 'the bound does not stop the build')


class TestTheOverlapCriterion(unittest.TestCase):
    """The thresholds, because a flat constant has now been wrong in BOTH directions."""

    @staticmethod
    def _limits():
        radii = dict(re.findall(r'([A-Z])\s+([0-9.]+)', re.search(
            r'array set _VDW_RADIUS \{([^}]*)\}', SRC).group(1)))
        frac = float(re.search(r'set _OVERLAP_FRACTION ([0-9.]+)', SRC).group(1))
        return {k: float(v) for k, v in radii.items()}, frac

    def test_the_radii_are_van_der_waals_not_vmd_display_values(self):
        """VMD reports 1.5 for carbon and 1.3 for oxygen -- ~0.2 A small. Using those would
        quietly loosen every threshold."""
        r, _ = self._limits()
        self.assertAlmostEqual(r['C'], 1.70, places=2)
        self.assertAlmostEqual(r['O'], 1.52, places=2)
        self.assertAlmostEqual(r['N'], 1.55, places=2)

    def test_a_carbon_carbon_clash_at_2_8_is_removed(self):
        r, f = self._limits()
        self.assertLess(2.8, f * (r['C'] + r['C']))

    def test_a_hydrogen_bond_at_2_8_is_KEPT(self):
        """The property no flat cutoff could deliver, and the reason both flat values failed."""
        r, f = self._limits()
        for a, b in (('O', 'N'), ('O', 'O')):
            self.assertGreater(2.8, f * (r[a] + r[b]),
                               f'a {a}...{b} H-bond at 2.8 A would be deleted')

    def test_the_fraction_is_inside_the_window_that_satisfies_BOTH(self):
        """NEGATIVE CONTROL for drifting the constant: outside 0.824-0.912 one of the two
        properties above silently stops holding."""
        _r, f = self._limits()
        self.assertGreater(f, 0.824)
        self.assertLess(f, 0.912)

    def test_the_search_radius_exceeds_every_threshold(self):
        """A search narrower than the widest limit silently misses clashes."""
        r, f = self._limits()
        search = float(re.search(r'set _OVERLAP_SEARCH ([0-9.]+)', SRC).group(1))
        self.assertGreater(search, f * 2 * max(r.values()))

    def test_an_unknown_element_falls_back_to_carbon(self):
        self.assertIn('set ri 1.70', SRC)
        self.assertIn('set rj 1.70', SRC)

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


class TestTheClashSelectionCoversTheWholeInsertedSystem(unittest.TestCase):
    """`noh`, never `protein and noh`.

    The embed deletes lipids and solvent that overlap the inserted system, and that deletion is
    what makes a piercing impossible: a bond cannot thread a ring without coming ~1.6 A from a
    ring atom (1.77 for a 7-ring), far inside the 2.89 A C...C deletion threshold -- evading it
    needs a ring radius of 2.79 A, about seventeen atoms, and `max_ring_size` is 7.

    That guarantee only covers what the selection covers.  Written `protein and noh` -- which
    reads more naturally and is what one would expect -- glycans would be invisible, and a
    **pierced glycan ring has no remedy at all**: it cannot be rotated reliably and cannot be
    deleted without breaking the glycan tree, so it is an unrecoverable build.  The selection is
    inclusive today by good fortune rather than by test, which is what this fixes.
    """

    def test_neither_site_restricts_the_inserted_system_to_protein(self):
        for ln in re.findall(r'atomselect \$\w+ "([^"]*)"', SRC):
            if 'noh' in ln:
                self.assertNotIn('protein', ln,
                                 f'selection "{ln}" excludes glycans from the clash test')

    def test_the_lipid_site_selects_all_heavy_atoms_of_the_inserted_molecule(self):
        self.assertIn('set pro_heavy [atomselect $protein "noh"]', SRC)

    def test_the_water_site_protects_everything_that_is_not_the_new_water(self):
        """Same rule on the solvent pass: glycans are on the protected side, not the deleted
        one, so a water overlapping a glycan goes and the glycan stays."""
        self.assertRegex(SRC, r'set pro_heavy_w \[atomselect \$\w+ "not \(segname [^)]*\) and noh"\]')
