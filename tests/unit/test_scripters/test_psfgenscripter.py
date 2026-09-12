import os
import unittest
from types import SimpleNamespace

import numpy as np

from pestifer.core.config import Config
from pestifer.molecule.transform import Transform
from pestifer.scripters import PsfgenScripter

class TestPsfgenScripter(unittest.TestCase):

    def setUp(self):
        self.c = Config().configure_new()
        self.p: PsfgenScripter = self.c.get_scripter('psfgen')

    def tearDown(self):
        self.c.RM.charmmff_content.clean_local_charmmff_files()

    def test_header(self):
        self.assertIsInstance(self.p, PsfgenScripter)
        self.p.newscript()
        self.p.B.banner('TESTING')
        self.p.writescript('test-state')
        self.assertTrue(os.path.isfile('pestifer-script.tcl'))
        os.remove(self.p.basename+'.tcl')
        self.p.newscript('testing')
        self.p.B.banner('TESTING')
        self.p.writescript('test-state')
        self.assertTrue(os.path.isfile('testing.tcl'))
        os.remove(self.p.basename+'.tcl')
        self.assertFalse(os.path.exists('testing.tcl'))
    
    def test_fetch_charmff_files(self):
        self.p.fetch_standard_charmm_topologies()
        self.assertTrue(os.path.isfile('top_all36_prot.rtf'))
        self.assertTrue(os.path.isfile('toppar_water_ions.str'))



class TestImageChainMapIsAuthoritative(unittest.TestCase):
    """
    A biological assembly that selects a *subset* of the asymmetric unit must build only the
    chains it selects.

    ``Transform.chainIDmap`` already names exactly those chains, but the stanza writers
    translate segment names with ``chainIDmap.get(seglabel, seglabel)``, so before this guard
    an unmapped segment fell through under its own A.U. name: 8DX0 assembly 1 selects chain A
    alone and pestifer built both protomers, both magnesiums, and both water sets -- silently,
    exit 0.  These pin the filter at the point of use.
    """

    def setUp(self):
        self.c = Config().configure_new()
        self.p: PsfgenScripter = self.c.get_scripter('psfgen')

    def tearDown(self):
        self.c.RM.charmmff_content.clean_local_charmmff_files()

    @staticmethod
    def _transform(chainIDmap):
        return Transform(index=1, tmat=np.identity(4, dtype=float),
                         applies_chainIDs=list(chainIDmap or {}), chainIDmap=chainIDmap)

    def test_includes_reads_the_map_as_authoritative(self):
        t = self._transform({'A': 'A', 'C': 'C', 'E': 'E'})
        self.assertTrue(t.includes('A'))
        self.assertFalse(t.includes('B'))
        # any label is sufficient: a glycan is mapped by structured segname, its chainID is the
        # parent protein chain
        self.assertTrue(t.includes('AG01', 'A'))
        self.assertFalse(t.includes('BG01', 'B'))
        # a None or empty label never matches on its own
        self.assertFalse(t.includes(None, ''))

    def test_empty_map_carries_no_information(self):
        # no assembly active (or the map is not generated yet) -- everything belongs
        for empty in ({}, None):
            self.assertTrue(self._transform(empty).includes('B'))

    def test_write_segments_skips_chains_the_image_omits(self):
        segments = SimpleNamespace(data=[
            SimpleNamespace(segname='A', chainID='A', segtype='protein'),
            SimpleNamespace(segname='B', chainID='B', segtype='protein'),
            SimpleNamespace(segname='AG01', chainID='A', segtype='glycan'),
            SimpleNamespace(segname='BG01', chainID='B', segtype='glycan'),
            SimpleNamespace(segname='C', chainID='C', segtype='ion'),
            SimpleNamespace(segname='D', chainID='D', segtype='ion'),
        ])
        written = []
        self.p.newscript('subset-assembly')
        self.p.write_segment = lambda segment, transform=None: written.append(segment.segname)

        t = self._transform({'A': 'A', 'AG01': 'AG01', 'C': 'C'})
        self.p.write_segments(segments, t, special=set())
        self.assertEqual(written, ['A', 'AG01', 'C'],
                         'segments absent from a non-empty image chain map must not be built')

        # and the unfiltered case is untouched: an empty map builds the whole A.U.
        written.clear()
        self.p.write_segments(segments, self._transform({}), special=set())
        self.assertEqual(written, ['A', 'B', 'AG01', 'BG01', 'C', 'D'])


class _Subsegments(list):
    """The slice of SubsegmentList the polymer stanza touches: list behaviour plus ``.data``."""
    @property
    def data(self):
        return self


class TestBackboneAcylatedFirstResidue(unittest.TestCase):
    """A chain that begins with a residue whose own topology acylates the backbone N (GLYM,
    N-myristoyl-glycine) must be written with ``first none``.

    CHARMM's lipid_prot stream sets no DEFA -- its ``DEFA FIRS NTER LAST CTER`` is commented out
    -- so the residue takes whatever default the previously read topology left in force.  Built
    directly in psfgen with NTER in force and no ``first none``, 1UPH's GLYM came out with an
    ``NH3`` nitrogen carrying HT1/HT2/HT3 *while still bonded to C1*, and psfgen raised nothing.
    pestifer happened to load that stream after a file ending ``DEFAULT FIRST NONE``, so its
    builds were right by load order alone.
    """

    def setUp(self):
        self.c = Config().configure_new()
        self.p: PsfgenScripter = self.c.get_scripter('psfgen')
        self.p.newscript('stanza_test')
        self.p._author_subsegment_pdb = lambda *a, **k: None

    def tearDown(self):
        self.c.RM.charmmff_content.clean_local_charmmff_files()
        for f in ('stanza_test.tcl',):
            if os.path.exists(f):
                os.remove(f)

    def _segment(self, first_resname, mutations=None, patches=None):
        from pestifer.molecule.residue import Residue, ResidueList
        from pestifer.molecule.atom import AtomList
        from pestifer.objs.mutation import MutationList
        from pestifer.objs.patch import PatchList
        from pestifer.objs.resid import ResID
        residues = ResidueList([
            Residue({'resname': rn, 'resid': ResID(rid), 'chainID': 'A', 'segtype': 'protein',
                     'segname': 'A', 'resolved': True, 'atoms': AtomList([])})
            for rn, rid in ((first_resname, 2), ('ALA', 3), ('ARG', 4))])
        sub = SimpleNamespace(state='RESOLVED', bounds=[0, 2], pdb=None, build=False,
                              num_items=lambda: 3, declare_buildable=lambda: None)
        return SimpleNamespace(
            segname='A', segtype='protein', chainID='A', specs={}, residues=residues,
            subsegments=_Subsegments([sub]),
            objmanager={'seq': {'mutations': mutations or MutationList([]),
                                'patches': patches or PatchList([])}})

    def _stanza(self, segment):
        t = Transform(index=1, tmat=np.identity(4, dtype=float), applies_chainIDs=['A'],
                      chainIDmap={'A': 'A'})
        self.p.write_polymer_stanza(segment, t)
        return self.p.B.byte_collector if hasattr(self.p.B, 'byte_collector') else str(self.p.B)

    def test_glym_first_residue_gets_first_none(self):
        text = self._stanza(self._segment('GLYM'))
        seg = text[text.index('segment A {'):text.index('}', text.index('segment A {'))]
        self.assertIn('first none', seg)

    def test_ordinary_first_residue_is_left_to_the_default(self):
        text = self._stanza(self._segment('GLY'))
        seg = text[text.index('segment A {'):text.index('}', text.index('segment A {'))]
        self.assertNotIn('first none', seg)

    def test_a_mutation_to_glym_at_the_first_position_counts(self):
        from pestifer.objs.mutation import Mutation, MutationList
        text = self._stanza(self._segment('GLY', mutations=MutationList([Mutation('A:GLY,2,GLYM')])))
        seg = text[text.index('segment A {'):text.index('}', text.index('segment A {'))]
        self.assertIn('first none', seg)

    def test_a_mutation_elsewhere_does_not(self):
        from pestifer.objs.mutation import Mutation, MutationList
        text = self._stanza(self._segment('GLY', mutations=MutationList([Mutation('A:ALA,3,GLYM')])))
        seg = text[text.index('segment A {'):text.index('}', text.index('segment A {'))]
        self.assertNotIn('first none', seg)
