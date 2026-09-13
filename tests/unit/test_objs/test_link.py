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
        # The chemically correct patch for each of 4zmj's 25 links, from what the sugars ARE:
        # NAG is beta-GlcNAc (equatorial at C1), MAN alpha-mannose (axial at C1, axial O2), BMA
        # beta-mannose.  This previously pinned {NGLA, NGLB, 12ab, 13bb, 14aa, 16BT} -- the choice
        # a nearest-reference-DIHEDRAL lookup made, which gave most beta-GlcNAc-Asn links the alpha
        # patch.  A torsion about the glycosidic bond is conformation and cannot see configuration.
        self.assertEqual(sorted((l.resname1, l.name1, l.resname2, l.patchname) for l in Lc), sorted(
            [('ASN', 'ND2', 'NAG', 'NGLB')] * 18 +
            [('NAG', 'O4', 'NAG', '14bb')] * 3 +
            [('NAG', 'O4', 'BMA', '14bb'),
             ('BMA', 'O3', 'MAN', '13ab'),
             ('BMA', 'O6', 'MAN', '16AT'),
             ('MAN', 'O2', 'MAN', '12aa')]))

    def test_identity_wins_over_a_distorted_ring_and_says_so(self):
        # Man A5 has alpha handedness at C1 but sits in a boat-like ring, where "axial" means
        # nothing and the coordinates read equatorial.  The residue's identity decides.
        L, r = self._fresh()
        with self.assertLogs('pestifer.objs.link', level='WARNING') as cm:
            L.assign_residues(r)
        man_man = [l for l in L if l.resname1 == 'MAN' and l.resname2 == 'MAN']
        self.assertEqual([l.patchname for l in man_man], ['12aa'])
        self.assertTrue(any('MAN A5 C1' in m and 'distorted' in m for m in cm.output), cm.output)


class TestGlycanPatchNamesMatchCharmmGeometry(unittest.TestCase):
    """
    CHARMM names each glycosidic patch by axial/equatorial ring geometry.  Every patch pestifer
    chooses between was built by psfgen from its own internal coordinates alone; measuring that
    geometry must give back the same name, for both the coordinate classifier and residue identity.
    The fixture is what makes this checkable without psfgen.
    """

    @classmethod
    def setUpClass(cls):
        import json
        from types import SimpleNamespace
        path = Path(__file__).parent.parent.parent / 'inputs' / 'glycan_patch_reference_geometry.json'
        cls.cases = json.load(open(path))['cases']
        cls.NS = SimpleNamespace

    def _residue(self, atoms, segname, resname):
        class _Atoms(list):
            def get(self, f):
                m = [a for a in self if f(a)]
                return None if not m else (m[0] if len(m) == 1 else m)
        pick = [self.NS(name=a['name'], x=a['xyz'][0], y=a['xyz'][1], z=a['xyz'][2])
                for a in atoms if a['segname'] == segname[0] and a['resid'] == segname[1]]
        return self.NS(resname=resname, chainID=segname[0], resid=self.NS(resid=segname[1]),
                       atoms=_Atoms(pick))

    def _labels(self, patch, case, label_fn):
        from pestifer.objs import link as linkmod
        atoms, r1, r2 = case['atoms'], case['residue1'], case['residue2']
        if case['kind'] == 'prot':
            aa = self._residue(atoms, ('P', '1'), r1)
            sugar = self._residue(atoms, ('G', '1'), r2)
            x = linkmod._atom(aa, linkmod._PROTEIN_GLYCOSYLATION[r1][1])
            return linkmod._PROTEIN_GLYCOSYLATION[r1][0] + label_fn(sugar, '1', x).upper()
        acceptor = self._residue(atoms, ('G', '1'), r1)
        donor = self._residue(atoms, ('G', '2'), r2)
        n = patch[1]
        sub = linkmod._atom(acceptor, f'O{n}')
        d = label_fn(donor, '1', sub)
        if n == '6':
            return '16' + ('AT' if d == 'a' else 'BT')
        if n == '1':   # CHARMM names 1<->1 from residue 1's C1 first
            return f'11{label_fn(acceptor, n, sub)}{d}'
        return f'1{n}{d}{label_fn(acceptor, n, sub)}'

    def test_geometry_reproduces_every_patch_name(self):
        from pestifer.objs import link as linkmod
        geom = lambda res, pos, sub: linkmod._anomeric_label(res, f'C{pos}', sub)
        self.assertEqual(len(self.cases), 23)
        for patch, case in self.cases.items():
            with self.subTest(patch=patch):
                self.assertEqual(self._labels(patch, case, geom), patch)

    def test_residue_identity_reproduces_every_patch_name(self):
        from pestifer.objs import link as linkmod
        for patch, case in self.cases.items():
            with self.subTest(patch=patch):
                ident = lambda res, pos, sub: linkmod._identity_label(res, pos)
                self.assertEqual(self._labels(patch, case, ident), patch)

    def test_set_patchname_itself_reproduces_every_patch_name(self):
        # the wiring, not just the helpers: set_patchname orders the letters itself, and 1<->1 is
        # named residue-1-first where every other link is donor-first
        for patch, case in self.cases.items():
            with self.subTest(patch=patch):
                atoms, r1, r2 = case['atoms'], case['residue1'], case['residue2']
                if case['kind'] == 'prot':
                    res1, res2 = self._residue(atoms, ('P', '1'), r1), self._residue(atoms, ('G', '1'), r2)
                    name1, seg1 = {'ASN': 'ND2', 'SER': 'OG', 'THR': 'OG1'}[r1], 'protein'
                else:
                    res1, res2 = self._residue(atoms, ('G', '1'), r1), self._residue(atoms, ('G', '2'), r2)
                    name1, seg1 = f'O{patch[1]}', 'glycan'
                L = Link(chainID1='A', resid1=ResID(1), name1=name1, chainID2='B', resid2=ResID(2), name2='C1')
                L.resname1, L.resname2, L.segtype1, L.segtype2 = r1, r2, seg1, 'glycan'
                L.residue1, L.residue2 = res1, res2
                L.set_patchname()
                self.assertEqual(L.patchname, patch)

    def test_axial_table_is_consistent_with_its_anomer_prefix(self):
        from pestifer.objs.link import _AXIAL_POSITIONS
        for resname, axial in _AXIAL_POSITIONS.items():
            with self.subTest(residue=resname):
                self.assertEqual(resname.startswith('A'), '1' in axial)
        # textbook stereochemistry, spot-checked: manno O2, galacto O4, allo O3 axial
        self.assertEqual(_AXIAL_POSITIONS['AMAN'], '12')
        self.assertEqual(_AXIAL_POSITIONS['BGAL'], '4')
        self.assertEqual(_AXIAL_POSITIONS['BALL'], '3')
        self.assertEqual(_AXIAL_POSITIONS['BGLCNA'], '')


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
