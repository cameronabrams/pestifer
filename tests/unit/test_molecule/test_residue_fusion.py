# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
The PDB deposits a lipidated amino acid as a ligand plus a LINK; CHARMM defines the two as one
residue.  Fusing them is what preserves a deposited modification's real conformation -- in 1uph
the myristate is sequestered in HIV-1 matrix's hydrophobic pocket, which no rebuild reproduces.

The failure without it is worse than a stopped build: ``MYR`` is itself a CHARMM residue (free
myristic acid), so a build can succeed and leave a DETACHED fatty acid beside the protein.
"""
import unittest
from types import SimpleNamespace

from pestifer.molecule.residue_fusion import fuse_linked_ligands


def _atom(name, resname, chainID, resid):
    return SimpleNamespace(name=name, resname=resname, chainID=chainID, resid=resid)


def _link(rn1, n1, c1, r1, rn2, n2, c2, r2):
    return SimpleNamespace(resname1=rn1, name1=n1, chainID1=c1, resid1=r1,
                           resname2=rn2, name2=n2, chainID2=c2, resid2=r2)


def _myr_gly():
    atoms = [_atom(n, 'MYR', 'A', 1) for n in ('C1', 'O1', 'C2', 'C3')]
    atoms += [_atom(n, 'MYR', 'A', 1) for n in ('H21', 'H22')]
    atoms += [_atom(n, 'GLY', 'A', 2) for n in ('N', 'CA', 'C', 'O')]
    return atoms


class TestFuseLinkedLigands(unittest.TestCase):

    def test_a_linked_myristate_becomes_part_of_the_glycine(self):
        atoms = _myr_gly()
        links = [_link('MYR', 'C1', 'A', 1, 'GLY', 'N', 'A', 2)]
        self.assertEqual(fuse_linked_ligands(atoms, links), 1)
        self.assertEqual({a.resname for a in atoms}, {'GLYM'})
        self.assertEqual({a.resid for a in atoms}, {2}, 'all atoms take the protein residue id')

    def test_the_consumed_link_is_removed(self):
        atoms = _myr_gly()
        links = [_link('MYR', 'C1', 'A', 1, 'GLY', 'N', 'A', 2)]
        fuse_linked_ligands(atoms, links)
        self.assertEqual(links, [], 'the bond is now internal; leaving the link patches a '
                                    'residue to itself')

    def test_heavy_atoms_are_carried_and_hydrogens_dropped(self):
        atoms = _myr_gly()
        fuse_linked_ligands(atoms, [_link('MYR', 'C1', 'A', 1, 'GLY', 'N', 'A', 2)])
        names = {a.name for a in atoms}
        self.assertTrue({'C1', 'O1', 'C2', 'C3'} <= names, 'heavy atoms must survive; they carry '
                                                           'the deposited conformation')
        self.assertNotIn('H21', names, 'ligand hydrogens use a different naming convention and '
                                       'are rebuilt by psfgen')

    def test_a_free_ligand_with_no_link_is_left_alone(self):
        atoms = _myr_gly()
        self.assertEqual(fuse_linked_ligands(atoms, []), 0)
        self.assertEqual({a.resname for a in atoms}, {'MYR', 'GLY'})

    def test_a_different_bond_between_the_same_residues_does_not_fuse(self):
        # the table names the bond, so an unrelated contact cannot trigger a fusion
        atoms = _myr_gly()
        links = [_link('MYR', 'O1', 'A', 1, 'GLY', 'CA', 'A', 2)]
        self.assertEqual(fuse_linked_ligands(atoms, links), 0)
        self.assertEqual({a.resname for a in atoms}, {'MYR', 'GLY'})

    def test_the_link_may_be_written_either_way_round(self):
        atoms = _myr_gly()
        links = [_link('GLY', 'N', 'A', 2, 'MYR', 'C1', 'A', 1)]
        self.assertEqual(fuse_linked_ligands(atoms, links), 1)
        self.assertEqual({a.resname for a in atoms}, {'GLYM'})

    def test_an_unlisted_pair_is_ignored(self):
        atoms = [_atom('C1', 'PLM', 'A', 1), _atom('SG', 'CYS', 'A', 2)]
        links = [_link('PLM', 'C1', 'A', 1, 'CYS', 'SG', 'A', 2)]
        self.assertEqual(fuse_linked_ligands(atoms, links), 0)


class TestFusionIsWiredIntoIngest(unittest.TestCase):
    """The tests above call the function directly, so they all pass if it is never called."""

    def test_asymmetricunit_fuses_before_grouping_and_before_segtypes(self):
        import inspect
        from pestifer.molecule import asymmetricunit as m
        src = inspect.getsource(m)
        self.assertIn('fuse_linked_ligands(atoms, links)', src)
        # order matters: fusing after grouping would require splicing Residue objects, and
        # fusing after apply_segtypes would classify the ligand on its own resname
        i_fuse = src.index('fuse_linked_ligands(atoms, links)')
        i_group = src.index('ResidueList.from_residuegrouped_atomlist(atoms)')
        i_seg = src.index('residues.apply_segtypes()')
        self.assertLess(i_fuse, i_group, 'fusion must precede residue grouping')
        self.assertLess(i_fuse, i_seg, 'fusion must precede segtype assignment')
