# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import unittest
from unittest import mock


from pathlib import Path

from pidibble.pdbparse import PDBParser
from pidibble.pdbrecord import PDBRecordDict

from pestifer.molecule.atom import AtomList
from pestifer.molecule.residue import Residue, ResidueList
from pestifer.objs.link import Link, LinkList, ic_reference_closest
from pestifer.objs.resid import ResID
from pestifer.psfutil.psfpatch import PSFLinkPatch

logger = logging.getLogger(__name__)

class TestLink(unittest.TestCase):

    def test_link_creation(self):
        link = Link(
            chainID1="A",
            resid1=ResID(1),
            name1="N",
            chainID2="B",
            resid2=ResID(2),
            name2="CA"
        )
        self.assertIsInstance(link, Link)
        self.assertEqual(link.chainID1, "A")
        self.assertEqual(link.resid1, ResID(1))
        self.assertEqual(link.name1, "N")
        self.assertEqual(link.chainID2, "B")
        self.assertEqual(link.resid2, ResID(2))
        self.assertEqual(link.name2, "CA")
        self.assertEqual(link._yaml_header, 'links')
        self.assertEqual(link._objcat, 'topol')
        # self.assertEqual(link.patchhead,1)
        self.assertEqual(repr(link), "Link(chainID1='A', resid1=ResID(resseqnum=1), name1='N', chainID2='B', resid2=ResID(resseqnum=2), name2='CA')")

    def test_link_from_shortcode(self):
        link = Link("A_1_N-B_2_CA")
        self.assertIsInstance(link, Link)
        self.assertEqual(link.chainID1, "A")
        self.assertEqual(link.resid1, ResID(1))
        self.assertEqual(link.name1, "N")
        self.assertEqual(link.chainID2, "B")
        self.assertEqual(link.resid2, ResID(2))
        self.assertEqual(link.name2, "CA")

    def test_link_from_psflinkpatch(self):
        patch = PSFLinkPatch(['NGLB', 'A:1', 'B:2'])
        link = Link(patch)
        self.assertIsInstance(link, Link)
        self.assertEqual(link.chainID1, 'A')
        self.assertEqual(link.resid1, ResID(1))
        self.assertEqual(link.name1, 'ND2')
        self.assertEqual(link.chainID2, 'B')
        self.assertEqual(link.resid2, ResID(2))
        self.assertEqual(link.name2, 'C1')

class TestLinkList(unittest.TestCase):

    def setUp(self):
        # set the inputs path to PACKAGE/tests/inputs
        self.inputs_path = Path(__file__).parent.parent.parent / 'inputs'

    def test_link_list_from_pdb(self):
        p = PDBParser(filepath=self.inputs_path / "4zmj.pdb").parse().parsed
        self.assertIsInstance(p, PDBRecordDict)
        L = LinkList.from_pdb(p)
        self.assertIsInstance(L, LinkList)
        self.assertGreater(len(L), 0)

    def test_link_list_from_cif(self):
        parsed = PDBParser(filepath=str(self.inputs_path / '4zmj.cif'), input_format='mmCIF').parse().parsed
        self.assertIsInstance(parsed, PDBRecordDict)
        L = LinkList.from_cif(parsed)
        self.assertIsInstance(L, LinkList)
        self.assertGreater(len(L), 0)
            
    def test_link_list_assign_residues(self):
        # make a list of 10 mock residues
        residues = ResidueList([Residue(chainID='A' if i%2==1 else 'B', segname='A' if i%2==1 else 'B', resid=ResID(i), resname='ALA', atoms=AtomList([]),segtype='protein', resolved=True) for i in range(1, 11)])

        # create a LinkList with 6 links
        links = LinkList([
            Link(chainID1='A', resid1=ResID(1), name1='N', chainID2='B', resid2=ResID(2), name2='CA'),
            Link(chainID1='A', resid1=ResID(3), name1='C', chainID2='B', resid2=ResID(4), name2='O'),
            Link(chainID1='A', resid1=ResID(5), name1='CA', chainID2='B', resid2=ResID(6), name2='CB'),
            Link(chainID1='A', resid1=ResID(7), name1='N', chainID2='B', resid2=ResID(8), name2='CA'),
            Link(chainID1='A', resid1=ResID(9), name1='C', chainID2='B', resid2=ResID(10), name2='O'),
            Link(chainID1='A', resid1=ResID(11), name1='CA', chainID2='B', resid2=ResID(12), name2='CB')
        ])

        # assign residues to the links
        ignored_by_ptnr1 = links.assign_objs_to_attr('residue1', residues, chainID='chainID1', resid='resid1')
        ignored_by_ptnr2 = links.assign_objs_to_attr('residue2', residues, chainID='chainID2', resid='resid2')

        self.assertEqual(links[0].residue1, residues[0])
        self.assertEqual(links[0].residue2, residues[1])
        self.assertEqual(links[1].residue1, residues[2])
        self.assertEqual(links[1].residue2, residues[3])
        self.assertEqual(links[2].residue1, residues[4])
        self.assertEqual(links[2].residue2, residues[5])
        self.assertEqual(len(ignored_by_ptnr1), 1)  # last link references residues that are not in the list
        self.assertEqual(len(ignored_by_ptnr2), 0)  # all links have valid residue2 since the previous link was removed by assign_objs_to_attr

class TestLinkOrientation(unittest.TestCase):
    """
    A LINK record may list a glycosidic bond either way round.  The PDB convention and the
    CHARMM ``PRES`` definitions both put the anomeric carbon second (``O4 -> C1``), but other
    producers -- Rosetta output is the case that prompted this -- sometimes write ``C1 -> O4``.
    Pestifer canonicalizes on ingest; these pin that it happens, that it is complete, and that
    it never fires on a link that already works.
    """

    def _link(self, n1, n2, st1, st2, rn1='BGLC', rn2='BGLC'):
        L = Link(chainID1='A', resid1=ResID(1), name1=n1,
                 chainID2='B', resid2=ResID(2), name2=n2)
        L.segtype1, L.segtype2, L.resname1, L.resname2 = st1, st2, rn1, rn2
        return L

    def test_paired_fields_are_derived_not_hand_listed(self):
        # A partner field left out of the swap corrupts atom identity silently, so the pairing
        # is derived from the model's fields.  Both naming shapes must be picked up.
        pairs = dict(Link._paired_fields())
        for a, b in [('chainID1', 'chainID2'), ('resid1', 'resid2'), ('name1', 'name2'),
                     ('atom1', 'atom2'), ('residue1', 'residue2'), ('segtype1', 'segtype2'),
                     ('altloc1', 'altloc2'), ('sym1', 'sym2'), ('segname1', 'segname2'),
                     ('resname1', 'resname2'),
                     ('ptnr1_label_asym_id', 'ptnr2_label_asym_id'),
                     ('ptnr1_auth_seq_id', 'ptnr2_auth_seq_id')]:
            self.assertEqual(pairs.get(a), b, f'{a} is not paired with {b}')
        # non-partner fields must never be swapped
        for solo in ('patchname', 'patchhead', 'link_distance', 'empty'):
            self.assertNotIn(solo, pairs)

    def test_reverse_is_involutive(self):
        L = self._link('C1', 'O4', 'glycan', 'glycan')
        before = L.shortcode()
        L.reverse(); L.reverse()
        self.assertEqual(L.shortcode(), before)

    def test_inverted_glycan_glycan_is_reversed(self):
        L = self._link('C1', 'O4', 'glycan', 'glycan')
        self.assertTrue(L.canonicalize_glycan_orientation())
        self.assertEqual(L.name1, 'O4')
        self.assertEqual(L.name2, 'C1')

    def test_inverted_glycan_protein_is_reversed(self):
        L = self._link('C1', 'ND2', 'glycan', 'protein', 'BGLC', 'ASN')
        self.assertTrue(L.canonicalize_glycan_orientation())
        self.assertEqual((L.name1, L.resname1), ('ND2', 'ASN'))
        self.assertEqual((L.name2, L.resname2), ('C1', 'BGLC'))

    def test_inverted_sialic_acid_is_reversed(self):
        # sialic acids are anomeric at C2, not C1
        L = self._link('C2', 'O6', 'glycan', 'glycan', 'ANE5', 'BGAL')
        self.assertTrue(L.canonicalize_glycan_orientation())
        self.assertEqual((L.name1, L.name2), ('O6', 'C2'))

    def test_canonical_links_are_left_alone(self):
        for n1, n2, st1, st2, rn1, rn2 in [('O4', 'C1', 'glycan', 'glycan', 'BGLC', 'BGLC'),
                                           ('ND2', 'C1', 'protein', 'glycan', 'ASN', 'BGLC'),
                                           ('O6', 'C2', 'glycan', 'glycan', 'BGAL', 'ANE5')]:
            L = self._link(n1, n2, st1, st2, rn1, rn2)
            self.assertFalse(L.canonicalize_glycan_orientation(), f'{n1}->{n2} was reversed')
            self.assertEqual((L.name1, L.name2), (n1, n2))

    def test_non_glycan_links_are_left_alone(self):
        for n1, n2, st1, st2, rn1, rn2 in [('ZN', 'NE2', 'ion', 'protein', 'ZN', 'HIS'),
                                           ('NE2', 'FE', 'protein', 'ligand', 'HIS', 'HEM'),
                                           ('C', 'N', 'protein', 'protein', 'ALA', 'GLY')]:
            L = self._link(n1, n2, st1, st2, rn1, rn2)
            self.assertFalse(L.canonicalize_glycan_orientation())

    def test_anomeric_carbon_on_both_sides_is_not_guessed(self):
        # not a glycosidic linkage; there is no basis for a direction, so leave it to fail loudly
        L = self._link('C1', 'C1', 'glycan', 'glycan')
        self.assertFalse(L.canonicalize_glycan_orientation())


class TestLinkOrientationAgainstRealStructure(unittest.TestCase):
    """
    The invariant that matters: the same bonds, written either way round, must produce the same
    patches and the same glycan tree.  4zmj carries 25 glycan links spanning seven patch types.
    """

    def setUp(self):
        self.inputs_path = Path(__file__).parent.parent.parent / 'inputs'

    def _fresh(self):
        p = PDBParser(filepath=str(self.inputs_path / '4zmj.pdb')).parse().parsed
        residues = ResidueList.from_residuegrouped_atomlist(AtomList.from_pdb(p))
        residues.apply_segtypes()
        return LinkList.from_pdb(p), residues

    @staticmethod
    def _as_written_by_an_inverting_producer(L):
        # built by crossing the fields directly, not by calling Link.reverse(), so the test does
        # not depend on the method it is meant to exercise
        return LinkList([Link(chainID1=l.chainID2, resid1=l.resid2, name1=l.name2,
                              chainID2=l.chainID1, resid2=l.resid1, name2=l.name1,
                              resname1=l.resname2, resname2=l.resname1,
                              altloc1=l.altloc2, altloc2=l.altloc1) for l in L])

    @staticmethod
    def _fingerprint(L):
        # patch, emission order, and the parent/child direction link_to established
        return sorted((l.patchname, l.patchhead,
                       (l.residue1.chainID, str(l.residue1.resid)),
                       (l.residue2.chainID, str(l.residue2.resid))) for l in L)

    def test_inverted_links_give_the_same_patches_and_tree(self):
        Lc, rc = self._fresh()
        Lc.assign_residues(rc)
        Li, ri = self._fresh()
        Li = self._as_written_by_an_inverting_producer(Li)
        Li.assign_residues(ri)

        self.assertEqual(len(Lc), 25)
        self.assertEqual(self._fingerprint(Lc), self._fingerprint(Li))
        # and the canonical run really did resolve patches, so the comparison is not
        # two piles of UNFOUND agreeing with each other
        self.assertNotIn('UNFOUND', {l.patchname for l in Lc})
        # These are the patches 4zmj's geometry actually selects.  They changed when the
        # periodicity bug in ic_reference_closest was fixed: the old code misclassified every
        # reference point, so the set pinned here previously was the buggy classification.
        # Within a family the patches are topologically identical -- same dele/ATOM/BOND, same
        # types and charges, differing only in IC seed values -- so the correction changes the
        # geometry unresolved atoms are built from, not the chemistry.
        self.assertEqual({l.patchname for l in Lc},
                         {'NGLA', 'NGLB', '12ab', '13bb', '14aa', '16BT'})


class TestOGlycanPatchSelection(unittest.TestCase):
    """
    Serine and threonine take DIFFERENT O-glycosylation patches.  SGPA/SGPB retype 1CB to CT2 (a
    CH2) and bond through 1OG; TGPA/TGPB retype it to CT1 (a CH) and bond through 1OG1.  Emitting
    the serine patch for a threonine asks NAMD for a CT1 CT2 HA1 angle that does not exist.
    """

    def _patches_offered_for(self, resname1, name1):
        """The patch names set_patchname actually offers the geometry lookup for this link."""
        L = Link(chainID1='A', resid1=ResID(1), name1=name1,
                 chainID2='B', resid2=ResID(2), name2='C1')
        L.resname1, L.resname2 = resname1, 'BGLC'
        L.segtype1, L.segtype2 = 'protein', 'glycan'
        L.residue1, L.residue2 = 'r1', 'r2'
        seen = set()

        def capture(res12, ICmap):
            for entry in ICmap:
                seen.update(entry['mapping'])
            return next(iter(seen))

        with mock.patch('pestifer.objs.link.ic_reference_closest', side_effect=capture):
            L.set_patchname()
        return seen

    def test_threonine_offers_the_threonine_patches(self):
        self.assertEqual(self._patches_offered_for('THR', 'OG1'), {'TGPA', 'TGPB'})

    def test_serine_still_offers_the_serine_patches(self):
        self.assertEqual(self._patches_offered_for('SER', 'OG'), {'SGPA', 'SGPB'})

    def test_threonine_patches_are_addressable_when_read_back_from_a_psf(self):
        # _from_psflinkpatch looks the patch up here; without these entries a PSF carrying a
        # TGPA patch raises KeyError on re-read, so emitting TGPA without adding them would be
        # worse than the bug it fixes
        for p in ('TGPA', 'TGPB'):
            self.assertIn(p, Link._patch_atomnames, f'{p} missing from _patch_atomnames')
            self.assertEqual(Link._patch_atomnames[p], ['OG1', 'C1'])

    def test_serine_and_threonine_use_different_hydroxyl_atoms(self):
        self.assertEqual(Link._patch_atomnames['SGPA'][0], 'OG')
        self.assertEqual(Link._patch_atomnames['TGPA'][0], 'OG1')


class TestICReferenceClosest(unittest.TestCase):
    """
    ``ic_reference_closest`` picks the patch whose reference IC values are nearest the measured
    ones, in periodic dihedral space.  The periodicity correction was two sequential ifs --
    a 180-degree SHIFT, not a wrap -- which sent a perfect match to the maximum per-component
    distance and a near-opposite to nearly zero.
    """

    # every mapping that appears in set_patchname, family by family
    FAMILIES = {
        'NGL': {'NGLA': [168.99], 'NGLB': [-70.91]},
        'SGP': {'SGPA': [45.37], 'SGPB': [19.87]},
        'TGP': {'TGPA': [69.9], 'TGPB': [33.16]},
        '11': {'11aa': [103.46, 103.54], '11ab': [121.75, 51.80], '11bb': [-56.58, -79.64]},
        '12': {'12aa': [-132.81, 47.16], '12ab': [115.32, 86.93],
               '12ba': [-133.78, 168.07], '12bb': [117.14, -168.07]},
        '13': {'13aa': [113.19, 65.46], '13ab': [-141.32, 65.46],
               '13ba': [-131.68, -100.16], '13bb': [-141.32, -130.16]},
        '14': {'14aa': [-86.29, 133.57], '14ab': [72.71, 48.64],
               '14ba': [-86.3, -130.97], '14bb': [81.86, -130.97]},
        '16': {'16AT': [71.24], '16BT': [-63.49]},
    }

    class _Res:
        class _Atoms:
            def get(self, f):
                return type('A', (), {'name': 'x', 'resname': 'R', 'resseqnum': 1})()
        atoms = _Atoms()

    def _pick(self, mapping, measured_degrees):
        import numpy as np
        icmaps = [{'ICatomnames': ['1A', '2B', '2C', '2D'],
                   'mapping': {k: v[i] for k, v in mapping.items()}}
                  for i in range(len(next(iter(mapping.values()))))]
        vals = iter([d * np.pi / 180.0 for d in measured_degrees])
        with mock.patch('pestifer.objs.link.measure_dihedral', side_effect=lambda *a: next(vals)):
            return ic_reference_closest([self._Res(), self._Res()], icmaps)

    def test_every_reference_point_selects_its_own_patch(self):
        """Fed a patch's own reference geometry, the lookup must return that patch.  Before the
        fix this failed for all 23 reference points across all eight families."""
        for family, mapping in self.FAMILIES.items():
            for patch, refs in mapping.items():
                with self.subTest(family=family, patch=patch):
                    self.assertEqual(self._pick(mapping, refs), patch)

    def test_wrap_treats_plus_and_minus_180_as_adjacent(self):
        # 179 and -179 are 2 degrees apart, not 358
        mapping = {'near': [179.0], 'far': [90.0]}
        self.assertEqual(self._pick(mapping, [-179.0]), 'near')


class TestOneToOneLinkICMap(unittest.TestCase):
    """The 1->1 branch had one IC entry keyed ``atomnames`` where every other entry -- and the
    loop that reads them -- uses ``ICatomnames``, so any 1->1 glycosidic link raised KeyError
    before it could be classified."""

    def test_no_ic_entry_uses_the_wrong_key(self):
        import inspect
        src = inspect.getsource(Link.set_patchname)
        self.assertNotIn("{'atomnames'", src.replace(' ', ''))
        self.assertNotIn("{'atomnames':", src)

    def test_one_to_one_link_reaches_the_geometry_lookup(self):
        L = Link(chainID1='A', resid1=ResID(1), name1='O1',
                 chainID2='B', resid2=ResID(2), name2='C1')
        L.resname1, L.resname2 = 'BGLC', 'BGLC'
        L.segtype1, L.segtype2 = 'glycan', 'glycan'
        L.residue1, L.residue2 = 'r1', 'r2'
        with mock.patch('pestifer.objs.link.ic_reference_closest', return_value='11aa') as m:
            L.set_patchname()
        m.assert_called_once()
        self.assertEqual(L.patchname, '11aa')
        # all three IC entries must have been offered, not two
        self.assertEqual(len(m.call_args[0][1]), 3)
