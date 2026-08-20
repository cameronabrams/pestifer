# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for psfgen's pure helpers and bookkeeping.

``test_psfgen.py`` covers the guards that stop a bad build (unset coordinates, unaliased ions) and
runs real psfgen; ``test_psfgen_preserve.py`` and ``test_psfgen_resegment.py`` cover the two
incoming-PSF modes.  This file takes what is left: the graph and geometry helpers that decide
*what* gets rotated, declashed or carried forward, none of which needs VMD.

These are the parts where an error produces a structure that builds and runs.  A missed rotatable
bond just means a glycan never declashes; a resid parsed as ``None`` silently drops a residue from
a loop.  Nothing downstream complains.
"""

import os
import tempfile
import unittest
from unittest import mock

import networkx as nx
import numpy as np

from pestifer.tasks.psfgen import PsfgenTask, _resid_seqnum


def _task(**attrs):
    t = PsfgenTask.__new__(PsfgenTask)
    t.taskname = 'psfgen'
    t.basename = 'test'
    t.specs = {}
    for k, v in attrs.items():
        setattr(t, k, v)
    return t


class _Atom:
    """Stand-in for a PSF atom, carrying only what the graph helpers consult."""

    def __init__(self, serial, name='C', element=None, is_root=False):
        self.serial = serial
        self.name = name
        self.element = element if element is not None else name[0]
        self.is_root = is_root

    def isH(self):
        return self.element == 'H'

    def is_pep(self, other):
        # a peptide bond joins a backbone carbonyl C to the next residue's amide N
        return {self.name, other.name} == {'C', 'N'}

    def __repr__(self):
        return f'<{self.name}{self.serial}>'


# --- resid coercion --------------------------------------------------------------------------------

class TestResidSeqnum(unittest.TestCase):
    """PDB resids carry an optional insertion code, so a resid is not always an int.  Returning
    None for something unparseable is the documented behavior; the risk is a *wrong* number."""

    def test_a_plain_integer_is_itself(self):
        self.assertEqual(_resid_seqnum(145), 145)

    def test_a_numeric_string_is_converted(self):
        self.assertEqual(_resid_seqnum('145'), 145)

    def test_an_insertion_code_is_stripped(self):
        """'52A' is residue 52, not residue 52A and not a parse failure."""
        self.assertEqual(_resid_seqnum('52A'), 52)

    def test_a_negative_resid_keeps_its_sign(self):
        # negative resids are legal and appear in real entries (e.g. expression tags)
        self.assertEqual(_resid_seqnum('-3'), -3)
        self.assertEqual(_resid_seqnum('-3B'), -3)

    def test_surrounding_whitespace_is_tolerated(self):
        self.assertEqual(_resid_seqnum('  17 '), 17)

    def test_something_unparseable_is_none_rather_than_a_guess(self):
        for bad in (None, '', 'A', 'abc'):
            with self.subTest(bad=bad):
                self.assertIsNone(_resid_seqnum(bad))


# --- rotatable bonds -------------------------------------------------------------------------------

class TestRotatableBonds(unittest.TestCase):
    """Which bonds a pendant group can be rotated about during declash.

    Serials are 1-based in the PSF and 0-based in the coordinate array, so every returned index is
    ``serial - 1``; an off-by-one here rotates the wrong atoms and quietly wrecks the structure.
    """

    @staticmethod
    def _chain(names, root_index=0):
        """A linear graph of atoms with the given names, rooted at one end."""
        atoms = [_Atom(i + 1, name=n) for i, n in enumerate(names)]
        atoms[root_index].is_root = True
        g = nx.Graph()
        g.add_nodes_from(atoms)
        for a, b in zip(atoms, atoms[1:]):
            g.add_edge(a, b)
        return g, atoms

    def test_a_linear_chain_rotates_about_its_interior_bonds(self):
        g, _atoms = self._chain(['C', 'O', 'O', 'C'])
        bonds, movers = PsfgenTask._rotatable_bonds(g)
        self.assertTrue(bonds)
        self.assertEqual(len(bonds), len(movers))

    def test_indices_are_zero_based(self):
        """The bond between serials 2 and 3 must come back as (1, 2).

        The root atom is a legitimate bond endpoint -- rotating about the bond that attaches the
        pendant group is exactly what swings the whole group -- so this checks the conversion
        directly rather than assuming the lowest index belongs to a non-root atom.
        """
        g, _atoms = self._chain(['C', 'O', 'O', 'C'])
        bonds, _movers = PsfgenTask._rotatable_bonds(g)
        self.assertIn((1, 2), bonds)
        self.assertTrue(all(i >= 0 for b in bonds for i in b))

    def test_a_bond_to_hydrogen_is_not_rotatable(self):
        """Rotating a terminal H about its own bond moves nothing."""
        g, _atoms = self._chain(['C', 'O', 'H'])
        bonds, _movers = PsfgenTask._rotatable_bonds(g)
        self.assertEqual(len(bonds), 0)

    def test_a_peptide_bond_is_not_rotatable(self):
        g, _atoms = self._chain(['C', 'N', 'O', 'O'])
        bonds, _movers = PsfgenTask._rotatable_bonds(g)
        endpoint_names = {(b[0], b[1]) for b in bonds}
        self.assertNotIn((0, 1), endpoint_names, 'the C-N peptide bond was treated as rotatable')

    def test_a_bond_rooted_on_both_sides_is_not_rotatable(self):
        """A loop anchored at both ends cannot turn about an interior bond."""
        g, atoms = self._chain(['C', 'O', 'O', 'C'])
        atoms[-1].is_root = True
        bonds, _movers = PsfgenTask._rotatable_bonds(g)
        self.assertEqual(len(bonds), 0)

    def test_the_movers_are_the_unrooted_side(self):
        g, _atoms = self._chain(['C', 'O', 'O', 'C'])
        bonds, movers = PsfgenTask._rotatable_bonds(g)
        i = bonds.index((1, 2))
        self.assertEqual(list(movers[i]), [3], 'only the far side of the bond should move')

    def test_a_bond_whose_far_side_is_a_lone_hydrogen_is_skipped(self):
        # rotating one H about the bond that carries it is a no-op
        g, _atoms = self._chain(['C', 'O', 'H'])
        bonds, _m = PsfgenTask._rotatable_bonds(g)
        self.assertEqual(bonds, [])

    def test_the_graph_is_left_intact(self):
        """Each candidate bond is cut and restored; a leaked removal would corrupt later passes."""
        g, _atoms = self._chain(['C', 'O', 'O', 'C'])
        before = sorted((a.serial, b.serial) for a, b in g.edges())
        PsfgenTask._rotatable_bonds(g)
        after = sorted((a.serial, b.serial) for a, b in g.edges())
        self.assertEqual(before, after)

    def test_ring_bonds_are_never_rotatable(self):
        """A ring bond is not a bridge, so turning it would tear the ring open.

        The bond *attaching* the ring to the root is still rotatable -- that one swings the whole
        ring -- so the assertion is about which bonds appear, not that none do.
        """
        atoms = [_Atom(i + 1, name='C') for i in range(4)]
        atoms[0].is_root = True
        g = nx.Graph()
        g.add_edges_from([(atoms[0], atoms[1]), (atoms[1], atoms[2]),
                          (atoms[2], atoms[3]), (atoms[3], atoms[1])])
        bonds, _m = PsfgenTask._rotatable_bonds(g)
        ring_bonds = {(1, 2), (2, 3), (1, 3)}
        self.assertFalse(ring_bonds & set(bonds), f'a ring bond was rotated: {bonds}')
        self.assertEqual(bonds, [(0, 1)], 'the attachment bond should still be rotatable')


# --- crotation bookkeeping -------------------------------------------------------------------------

class TestCrotationTargetedResids(unittest.TestCase):
    """Residues a user positioned by hand must be left as built.

    Crotations run before the terminal-tail modeler, so a tail a user folded into a helix with an
    ALPHA crotation would otherwise be re-modeled and the user's work silently discarded.  The
    result is expanded over the biological-assembly images, because a crotation written against
    chain A applies to every image of chain A.
    """

    @staticmethod
    def _crot(chain, r1, r2=None):
        c = mock.Mock()
        c.chainID = chain
        c.resid1 = mock.Mock(resseqnum=r1)
        c.resid2 = mock.Mock(resseqnum=r2) if r2 is not None else None
        return c

    def _resids(self, crots=None, irots=None, chainmaps=None):
        objmgr = {'coord': {'crotations': crots or [], 'irotations': irots or []}}
        t = _task(objmanager=objmgr)
        ba = mock.Mock()
        ba.transforms = [mock.Mock(chainIDmap=m) for m in (chainmaps or [{}])]
        t.base_molecule = mock.Mock(active_biological_assembly=ba)
        return t._crotation_targeted_resids()

    def test_no_object_manager_yields_nothing(self):
        self.assertEqual(_task(objmanager=None)._crotation_targeted_resids(), {})

    def test_no_rotations_yields_nothing(self):
        self.assertEqual(self._resids(), {})

    def test_a_single_residue_rotation(self):
        self.assertEqual(self._resids([self._crot('A', 10)]), {'A': {10}})

    def test_a_residue_range_is_expanded_inclusively(self):
        self.assertEqual(self._resids([self._crot('A', 10, 13)]), {'A': {10, 11, 12, 13}})

    def test_a_reversed_range_is_still_expanded(self):
        self.assertEqual(self._resids([self._crot('A', 13, 10)]), {'A': {10, 11, 12, 13}})

    def test_irotations_count_too(self):
        self.assertEqual(self._resids(irots=[self._crot('A', 5)]), {'A': {5}})

    def test_the_result_is_expanded_over_assembly_images(self):
        """A crotation on chain A applies to every image of A, under each image's own label."""
        got = self._resids([self._crot('A', 7)],
                           chainmaps=[{'A': 'A'}, {'A': 'B'}, {'A': 'C'}])
        self.assertEqual(got, {'A': {7}, 'B': {7}, 'C': {7}})

    def test_an_unmapped_chain_keeps_its_own_label(self):
        self.assertEqual(self._resids([self._crot('X', 3)], chainmaps=[{'A': 'B'}]), {'X': {3}})

    def test_a_rotation_naming_no_chain_is_ignored(self):
        c = self._crot(None, 5)
        self.assertEqual(self._resids([c]), {})

    def test_rotations_on_several_chains_are_kept_apart(self):
        got = self._resids([self._crot('A', 1), self._crot('B', 2)])
        self.assertEqual(got, {'A': {1}, 'B': {2}})


# --- incoming stream files -------------------------------------------------------------------------

class TestIncomingStreamFiles(unittest.TestCase):
    """Which CHARMM stream files an imported PSF needs.

    Prefer what a prior continuation or merge already registered; only fall back to reading the
    PSF's own topology remarks, which is guesswork by comparison.
    """

    PSF_HEADER = ('PSF EXT CMAP\n'
                  '\n'
                  '       4 !NTITLE\n'
                  ' REMARKS topology /path/to/top_all36_prot.rtf\n'
                  ' REMARKS topology /path/to/toppar_all36_carb_glycopeptide.str\n'
                  ' REMARKS topology /elsewhere/toppar_water_ions.str\n'
                  ' REMARKS topology /path/to/toppar_all36_carb_glycopeptide.str\n'
                  '\n'
                  '     100 !NATOM\n'
                  ' REMARKS topology /never/reached.str\n')

    def _resolve(self, registered=None, present=()):
        with tempfile.TemporaryDirectory() as d:
            cwd = os.getcwd()
            os.chdir(d)
            try:
                with open('in.psf', 'w') as f:
                    f.write(self.PSF_HEADER)
                for name in present:
                    open(name, 'w').close()
                arts = None
                if registered is not None:
                    arts = [mock.Mock(name=n) for n in registered]
                    for a, n in zip(arts, registered):
                        a.name = n
                t = _task()
                t.get_current_artifact = lambda k: arts
                cc = mock.Mock()
                cc.copy_charmmfile_local = lambda n: open(n, 'w').close()
                t.resource_manager = mock.Mock(charmmff_content=cc)
                return t._incoming_stream_files('in.psf')
            finally:
                os.chdir(cwd)

    def test_registered_artifacts_are_preferred_over_parsing(self):
        got = self._resolve(registered=['a.str', 'b.str'], present=['a.str', 'b.str'])
        self.assertEqual(got, ['a.str', 'b.str'])

    def test_a_registered_file_that_is_gone_is_dropped(self):
        got = self._resolve(registered=['a.str', 'missing.str'], present=['a.str'])
        self.assertEqual(got, ['a.str'])

    def test_stream_files_are_read_from_the_psf_remarks(self):
        got = self._resolve()
        self.assertIn('toppar_all36_carb_glycopeptide.str', got)
        self.assertIn('toppar_water_ions.str', got)

    def test_rtf_topologies_are_not_stream_files(self):
        self.assertNotIn('top_all36_prot.rtf', self._resolve())

    def test_a_repeated_remark_is_listed_once(self):
        got = self._resolve()
        self.assertEqual(got.count('toppar_all36_carb_glycopeptide.str'), 1)

    def test_parsing_stops_at_the_atom_section(self):
        """The remarks block ends at !NATOM; reading on would scan the whole file for nothing."""
        self.assertNotIn('reached.str', self._resolve())

    def test_only_the_basename_is_used(self):
        """The paths in a PSF's remarks are wherever it was built, which is not here."""
        for name in self._resolve():
            self.assertEqual(name, os.path.basename(name))


# --- ring-piercing resolution ------------------------------------------------------------------------

class _PiercingCase(unittest.TestCase):
    """A ring piercing has no atomic clash, so no declash pass can see it: a bond simply threads
    through a ring.  Both resolvers rotate the offending branch out and report whatever they
    cannot clear, non-fatally, so a downstream ``ring_check`` can try harder."""

    @staticmethod
    def _piercing(piercee_res='PRO', piercer_seg='G1'):
        return {'piercee': {'resname': piercee_res, 'segname': 'A', 'resid': 42},
                'piercer': {'resname': 'BGLC', 'segname': piercer_seg, 'resid': 7,
                            'segtype': 'glycan'}}

    def _task_with_state(self):
        t = _task()
        state = mock.Mock()
        state.psf.name = 'in.psf'
        state.pdb.name = 'in.pdb'
        state.xsc = None
        t.get_current_artifact = lambda k: state
        t.next_basename = mock.Mock(side_effect=lambda n=None: setattr(t, 'basename', 'out'))
        t.register = mock.Mock()
        return t


class TestResolveTailPiercings(_PiercingCase):

    def _resolve(self, serials, fixed=None, unresolved=()):
        t = self._task_with_state()
        with mock.patch('pestifer.psfutil.ring_resolve.resolve_model_built_piercings',
                        return_value=(fixed, list(unresolved))) as engine, \
             mock.patch('pestifer.tasks.psfgen.shutil.copy') as copy:
            t.resolve_tail_piercings(serials)
        return t, engine, copy

    def test_no_modeled_tail_atoms_means_nothing_to_check(self):
        t, engine, _c = self._resolve([])
        engine.assert_not_called()
        t.next_basename.assert_not_called()

    def test_a_clean_scan_registers_nothing(self):
        t, engine, copy = self._resolve([1, 2, 3], fixed=None)
        engine.assert_called_once()
        copy.assert_not_called()
        t.register.assert_not_called()

    def test_a_resolved_piercing_replaces_the_coordinates(self):
        t, _e, copy = self._resolve([1, 2, 3], fixed='fixed.pdb')
        copy.assert_called_once_with('fixed.pdb', 'out.pdb')
        t.register.assert_called_once()

    def test_the_tail_atoms_are_passed_as_the_movers(self):
        """A tail has a free end, so its distal portion can swing without constraint -- which is
        only true if the engine is told the tail atoms are what may move."""
        _t, engine, _c = self._resolve([5, 6, 7])
        self.assertEqual(engine.call_args.args[4], [5, 6, 7])
        self.assertEqual(engine.call_args.kwargs['mover_kind'], 'tail')

    def test_an_unresolvable_piercing_warns_and_names_both_residues(self):
        with self.assertLogs('pestifer.tasks.psfgen', level='WARNING') as cm:
            self._resolve([1], unresolved=[self._piercing()])
        out = ''.join(cm.output)
        self.assertIn('PRO', out)
        self.assertIn('A-42', out)
        self.assertIn('ring_check', out, 'the warning should say what to do next')

    def test_an_unresolvable_piercing_is_not_fatal(self):
        t, _e, _c = self._resolve([1], unresolved=[self._piercing()])
        self.assertIsNotNone(t)


class TestResolveGlycanPiercings(_PiercingCase):

    @staticmethod
    def _specs(check=True):
        return {'source': {'sequence': {'glycans': {'declash': {'check_piercings': check}}}}}

    def _resolve(self, nglycans=1, check=True, fixed=None, unresolved=()):
        t = self._task_with_state()
        t.specs = self._specs(check)
        t.declash_counts = {'glycan': nglycans}
        with mock.patch('pestifer.psfutil.ring_resolve.resolve_glycan_piercings',
                        return_value=(fixed, list(unresolved))) as engine, \
             mock.patch('pestifer.tasks.psfgen.shutil.copy') as copy:
            t.resolve_glycan_piercings()
        return t, engine, copy

    def test_a_structure_with_no_glycans_is_not_scanned(self):
        _t, engine, _c = self._resolve(nglycans=0)
        engine.assert_not_called()

    def test_the_scan_can_be_switched_off(self):
        _t, engine, _c = self._resolve(check=False)
        engine.assert_not_called()

    def test_glycans_are_scanned_even_with_no_declash_cycles(self):
        """The clash-count declash cannot see a piercing, so this pass is independent of it."""
        _t, engine, _c = self._resolve(nglycans=3)
        engine.assert_called_once()

    def test_a_resolved_piercing_replaces_the_coordinates(self):
        t, _e, copy = self._resolve(fixed='fixed.pdb')
        copy.assert_called_once_with('fixed.pdb', 'out.pdb')
        t.register.assert_called_once()

    def test_protein_and_glycan_rings_are_both_checked(self):
        """A grafted glycan can spear a sugar ring as readily as a proline."""
        _t, engine, _c = self._resolve()
        self.assertEqual(engine.call_args.kwargs['segtypes'], ['protein', 'glycan'])

    def test_an_unresolvable_piercing_names_the_piercer_segment_type(self):
        with self.assertLogs('pestifer.tasks.psfgen', level='WARNING') as cm:
            self._resolve(unresolved=[self._piercing()])
        out = ''.join(cm.output)
        self.assertIn('BGLC', out)
        self.assertIn('glycan', out)


# --- topology collection ---------------------------------------------------------------------------

class TestResiTopologies(unittest.TestCase):
    """Every residue name in the structure must resolve to a CHARMM topology file.

    Failing loudly here matters: psfgen's own error for an unknown residue arrives much later and
    reads as a segment-building failure rather than a missing parameter set.
    """

    def _collect(self, resnames, topfile_of=None):
        t = _task()
        residues = [mock.Mock(resname=n) for n in resnames]
        t.base_molecule = mock.Mock()
        t.base_molecule.asymmetric_unit.residues.data = residues
        cc = mock.Mock()
        cc.get_topfile_of_resname = (topfile_of if topfile_of is not None
                                     else (lambda n: f'top_{n.lower()}.rtf'))
        t.resource_manager = mock.Mock(charmmff_content=cc)
        return t.resi_topologies()

    def test_each_distinct_residue_contributes_its_topology(self):
        got = self._collect(['ALA', 'TIP3'])
        self.assertCountEqual(got, ['top_ala.rtf', 'top_tip3.rtf'])

    def test_a_shared_topology_is_listed_once(self):
        got = self._collect(['ALA', 'GLY', 'SER'], topfile_of=lambda n: 'top_all36_prot.rtf')
        self.assertEqual(got, ['top_all36_prot.rtf'])

    def test_repeated_residues_do_not_repeat_the_topology(self):
        self.assertEqual(len(self._collect(['ALA'] * 50)), 1)

    def test_an_unresolvable_residue_is_a_build_error_naming_it(self):
        from pestifer.core.errors import PestiferBuildError
        with self.assertRaises(PestiferBuildError) as cm:
            self._collect(['ALA', 'XYZ'], topfile_of=lambda n: None if n == 'XYZ' else 'x.rtf')
        self.assertIn('XYZ', str(cm.exception))


class TestPatchTopologies(unittest.TestCase):
    """Patches and links carry their own topology requirements, and unlike residues a missing one
    is not fatal -- the patch may be defined in a file already loaded."""

    def _collect(self, patches=(), links=(), topfile_of=None):
        t = _task()
        objmanager = {'seq': {'patches': [mock.Mock(patchname=p) for p in patches]},
                      'topol': {'links': [mock.Mock(patchname=l) for l in links]}}
        t.base_molecule = mock.Mock(objmanager=objmanager)
        cc = mock.Mock()
        cc.get_topfile_of_resname = (topfile_of if topfile_of is not None
                                     else (lambda n: f'top_{n.lower()}.str'))
        t.resource_manager = mock.Mock(charmmff_content=cc)
        return t.patch_topologies()

    def test_nothing_to_patch_needs_nothing(self):
        self.assertEqual(self._collect(), [])

    def test_a_patch_contributes_its_topology(self):
        self.assertEqual(self._collect(patches=['DISU']), ['top_disu.str'])

    def test_a_link_contributes_its_topology(self):
        self.assertEqual(self._collect(links=['NGLB']), ['top_nglb.str'])

    def test_patches_and_links_are_merged_without_duplicates(self):
        got = self._collect(patches=['NGLB'], links=['NGLB'])
        self.assertEqual(got, ['top_nglb.str'])

    def test_an_unresolvable_patch_is_skipped_rather_than_fatal(self):
        got = self._collect(patches=['DISU', 'UNKNOWN'],
                            topfile_of=lambda n: 'top_disu.str' if n == 'DISU' else None)
        self.assertEqual(got, ['top_disu.str'])


# --- declash dispatch ------------------------------------------------------------------------------

class TestDeclashSegtypeDispatch(unittest.TestCase):
    """Each segment type has its own declasher; sending one to the wrong handler would declash
    with the wrong notion of what is rotatable."""

    def _dispatch(self, segtype):
        t = _task()
        for name in ('declash_protein_loops', 'declash_na_loops', 'declash_glycans'):
            setattr(t, name, mock.Mock())
        t.declash_segtype({'maxcycles': 5}, segtype=segtype)
        return t

    def test_protein_goes_to_the_loop_declasher(self):
        t = self._dispatch('protein')
        t.declash_protein_loops.assert_called_once()
        t.declash_glycans.assert_not_called()

    def test_nucleic_acid_goes_to_its_own_declasher(self):
        self._dispatch('nucleicacid').declash_na_loops.assert_called_once()

    def test_glycan_goes_to_the_pendant_declasher(self):
        self._dispatch('glycan').declash_glycans.assert_called_once()

    def test_an_unknown_segment_type_is_refused(self):
        from pestifer.core.errors import PestiferBuildError
        with self.assertRaises(PestiferBuildError) as cm:
            self._dispatch('lipid')
        self.assertIn('lipid', str(cm.exception))


class TestDeclashCounts(unittest.TestCase):
    """What ``declash`` decides to run.

    Counts are per asymmetric unit but the structure being declashed is the assembly, so every
    count is multiplied by the image count.  A build of a trimer declashes three of each loop.
    """

    @staticmethod
    def _specs(maxcycles=5, glycan_cycles=5, min_loop_length=4):
        declash = {'maxcycles': maxcycles, 'include_C_termini': False, 'model_tails': False}
        return {'source': {'sequence': {
            'loops': {'min_loop_length': min_loop_length, 'declash': dict(declash)},
            'glycans': {'declash': dict(declash, maxcycles=glycan_cycles)}}}}

    def _declash(self, loop_counts=None, nglycans=0, num_images=1, **specs_kw):
        t = _task(specs=self._specs(**specs_kw))
        t.base_molecule = mock.Mock()
        t.base_molecule.loop_counts.return_value = dict(loop_counts or
                                                        {'protein': 0, 'nucleicacid': 0})
        t.base_molecule.num_images.return_value = num_images
        t.base_molecule.nglycans.return_value = nglycans
        t.declash_segtype = mock.Mock()
        t.model_terminal_tails = mock.Mock()
        t.resolve_glycan_piercings = mock.Mock()
        t.declash()
        return t

    def test_nothing_to_declash_dispatches_nothing(self):
        t = self._declash()
        t.declash_segtype.assert_not_called()

    def test_a_protein_loop_is_declashed(self):
        t = self._declash({'protein': 1, 'nucleicacid': 0})
        self.assertEqual([c.kwargs['segtype'] for c in t.declash_segtype.call_args_list],
                         ['protein'])

    def test_counts_are_multiplied_by_the_assembly_images(self):
        """A loop in the asymmetric unit appears once per image of its chain."""
        t = self._declash({'protein': 2, 'nucleicacid': 0}, num_images=3)
        self.assertEqual(t.declash_counts['protein'], 6)

    def test_glycans_are_counted_separately_and_also_per_image(self):
        t = self._declash({'protein': 0, 'nucleicacid': 0}, nglycans=4, num_images=2)
        self.assertEqual(t.declash_counts['glycan'], 8)
        self.assertEqual([c.kwargs['segtype'] for c in t.declash_segtype.call_args_list],
                         ['glycan'])

    def test_a_zero_cycle_budget_skips_the_pass_even_when_there_is_work(self):
        t = self._declash({'protein': 5, 'nucleicacid': 0}, maxcycles=0)
        t.declash_segtype.assert_not_called()

    def test_the_minimum_loop_length_reaches_the_counter(self):
        t = self._declash(min_loop_length=7)
        t.base_molecule.loop_counts.assert_called_once_with(min_loop_length=7)

    def test_tails_are_modeled_and_glycan_piercings_resolved_regardless(self):
        """Both run whether or not anything was declashed: a grafted glycan can thread a ring
        with no atomic clash at all, so the clash-count pass above cannot see it."""
        t = self._declash()
        t.model_terminal_tails.assert_called_once()
        t.resolve_glycan_piercings.assert_called_once()


class TestProteinLoopCTerminusHandoff(unittest.TestCase):
    """Who owns a C-terminal run: the phi-wiggle declasher or the free-tail modeler.

    Both would happily process it, and the result of doing both is the tail modeler's careful
    conformation being wiggled afterwards by a declasher that treats it as a loop.
    """

    def _include_c(self, include_C_termini, model_tails):
        t = _task()
        t.next_basename = mock.Mock()
        t.get_current_artifact = mock.Mock(return_value=mock.Mock())
        t.register = mock.Mock()
        t._run_loop_declash = mock.Mock()
        captured = {}

        def fake_enumerate(include_c_termini=False):
            captured['include_c'] = include_c_termini
            return []

        t._enumerate_protein_loops = fake_enumerate
        t.declash_protein_loops({'maxcycles': 5, 'include_C_termini': include_C_termini,
                                 'model_tails': model_tails})
        return captured.get('include_c')

    def test_the_tail_modeler_owns_termini_when_it_is_on(self):
        self.assertFalse(self._include_c(include_C_termini=True, model_tails=True))

    def test_the_declasher_takes_termini_when_the_modeler_is_off(self):
        self.assertTrue(self._include_c(include_C_termini=True, model_tails=False))

    def test_termini_are_left_alone_when_not_requested(self):
        self.assertFalse(self._include_c(include_C_termini=False, model_tails=False))

    def test_a_zero_cycle_budget_does_no_work_at_all(self):
        t = _task()
        t.next_basename = mock.Mock()
        t._run_loop_declash = mock.Mock()
        t.declash_protein_loops({'maxcycles': 0, 'include_C_termini': False})
        t._run_loop_declash.assert_not_called()
        t.next_basename.assert_not_called()


# --- loop enumeration ------------------------------------------------------------------------------

class TestEnumerateProteinLoops(unittest.TestCase):
    """Which model-built gaps get declashed, once per biological-assembly image.

    Two exclusions matter.  A gap shorter than ``min_loop_length`` is not worth the declash pass,
    and a C-terminal run is not a loop at all -- it is a tail, dangling into solvent rather than
    bridging two anchors, and is handled by the tail modeler instead.
    """

    def _segment(self, segname='A', segtype='protein', subsegs=(), n_residues=50):
        S = mock.Mock()
        S.segname = segname
        S.segtype = segtype
        sub = [mock.Mock(state=st, bounds=b, **{'num_items.return_value': b[1] - b[0] + 1})
               for st, b in subsegs]
        S.subsegments.data = sub
        S.subsegments.__len__ = mock.Mock(return_value=len(sub))
        S.subsegments.index = lambda x: sub.index(x)
        S.residues.data = [mock.Mock(resid=mock.Mock(resseqnum=i + 1)) for i in range(n_residues)]
        return S

    def _enumerate(self, segments, chainmaps=None, min_length=4, include_c_termini=False):
        t = _task(min_loop_length=min_length)
        ba = mock.Mock()
        ba.transforms = [mock.Mock(chainIDmap=m) for m in (chainmaps or [{}])]
        au = mock.Mock()
        au.segments.data = segments
        t.base_molecule = mock.Mock(active_biological_assembly=ba, asymmetric_unit=au)
        return t._enumerate_protein_loops(include_c_termini=include_c_termini)

    def test_a_resolved_segment_has_no_loops(self):
        seg = self._segment(subsegs=[('RESOLVED', (0, 49))])
        self.assertEqual(self._enumerate([seg]), [])

    def test_an_interior_missing_run_is_a_loop(self):
        seg = self._segment(subsegs=[('RESOLVED', (0, 9)), ('MISSING', (10, 19)),
                                     ('RESOLVED', (20, 49))])
        got = self._enumerate([seg])
        self.assertEqual(len(got), 1)
        segname, resids = got[0]
        self.assertEqual(segname, 'A')
        self.assertEqual(resids, list(range(11, 21)))

    def test_a_short_gap_is_not_worth_declashing(self):
        seg = self._segment(subsegs=[('RESOLVED', (0, 9)), ('MISSING', (10, 11)),
                                     ('RESOLVED', (12, 49))])
        self.assertEqual(self._enumerate([seg], min_length=4), [])

    def test_a_c_terminal_run_is_a_tail_not_a_loop(self):
        seg = self._segment(subsegs=[('RESOLVED', (0, 39)), ('MISSING', (40, 49))])
        self.assertEqual(self._enumerate([seg]), [])

    def test_a_c_terminal_run_can_be_included_on_request(self):
        seg = self._segment(subsegs=[('RESOLVED', (0, 39)), ('MISSING', (40, 49))])
        self.assertEqual(len(self._enumerate([seg], include_c_termini=True)), 1)

    def test_non_protein_segments_are_ignored(self):
        seg = self._segment(segtype='glycan', subsegs=[('RESOLVED', (0, 9)),
                                                       ('MISSING', (10, 19)),
                                                       ('RESOLVED', (20, 49))])
        self.assertEqual(self._enumerate([seg]), [])

    def test_each_assembly_image_gets_its_own_entry(self):
        """A loop in chain A exists in every image of A, and each image declashes separately."""
        seg = self._segment(subsegs=[('RESOLVED', (0, 9)), ('MISSING', (10, 19)),
                                     ('RESOLVED', (20, 49))])
        got = self._enumerate([seg], chainmaps=[{'A': 'A'}, {'A': 'B'}])
        self.assertEqual([g[0] for g in got], ['A', 'B'])
        self.assertEqual(got[0][1], got[1][1], 'the same residues in each image')

    def test_several_gaps_in_one_segment_are_all_reported(self):
        seg = self._segment(subsegs=[('RESOLVED', (0, 9)), ('MISSING', (10, 19)),
                                     ('RESOLVED', (20, 29)), ('MISSING', (30, 39)),
                                     ('RESOLVED', (40, 49))])
        self.assertEqual(len(self._enumerate([seg])), 2)


# --- remark stripping ------------------------------------------------------------------------------

class TestStripRemarks(unittest.TestCase):
    """psfgen writes REMARK lines into its PDB; downstream fixed-column readers do better without
    them."""

    def _strip(self, text, exists=True, name='out.pdb'):
        with tempfile.TemporaryDirectory() as d:
            cwd = os.getcwd()
            os.chdir(d)
            try:
                if exists:
                    with open(name, 'w') as f:
                        f.write(text)
                state = mock.Mock()
                state.pdb.name = name
                state.psf = mock.Mock()
                t = _task()
                t.get_current_artifact = lambda k: state
                t.register = mock.Mock()
                t.strip_remarks()
                return (open(name).read() if os.path.exists(name) else None), t
            finally:
                os.chdir(cwd)

    BODY = ('REMARK original generated coordinate pdb file\n'
            'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00\n'
            'REMARK something else\n'
            'ATOM      2  CA  ALA A   1       1.000   1.000   1.000  1.00  0.00\n'
            'END\n')

    def test_remark_lines_are_removed(self):
        out, _t = self._strip(self.BODY)
        self.assertNotIn('REMARK', out)

    def test_coordinates_are_untouched(self):
        out, _t = self._strip(self.BODY)
        self.assertEqual(len([l for l in out.splitlines() if l.startswith('ATOM')]), 2)

    def test_the_stripped_state_is_re_registered(self):
        _out, t = self._strip(self.BODY)
        t.register.assert_called_once()

    def test_a_file_with_no_remarks_is_unchanged(self):
        body = 'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00\nEND\n'
        out, _t = self._strip(body)
        self.assertEqual(out, body)

    def test_a_missing_file_warns_rather_than_raising(self):
        with self.assertLogs('pestifer.tasks.psfgen', level='WARNING') as cm:
            out, t = self._strip('', exists=False)
        self.assertIn('does not exist', ''.join(cm.output))
        t.register.assert_not_called()

    def test_a_word_beginning_with_remark_elsewhere_is_not_a_remark_line(self):
        # only line-leading REMARK counts; an atom name is not a comment
        body = 'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00  REMARKS\n'
        out, _t = self._strip(body)
        self.assertIn('ATOM', out)


if __name__ == '__main__':
    unittest.main()
