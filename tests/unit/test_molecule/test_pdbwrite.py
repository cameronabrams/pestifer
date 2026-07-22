import os
import tempfile
import unittest

from pestifer.molecule.atom import Atom, AtomList, Hetatm
from pestifer.molecule.pdbwrite import write_atoms_pdb
from pestifer.objs.resid import ResID
from pidibble.pdbparse import PDBParser


def _atom(cls=Atom, **kw):
    base = dict(serial=1, name='CA', altloc=' ', resname='ALA', chainID='A', resid=ResID(1),
                x=1.0, y=2.0, z=3.0, occ=1.0, beta=0.0, elem='C', charge=' ', segname='PROA')
    base.update(kw)
    return cls(**base)


class TestWriteAtomsPdb(unittest.TestCase):

    def test_coords_pinned_at_standard_columns(self):
        line = write_atoms_pdb([_atom(x=12.345, y=-6.789, z=0.5)], end=False)[0]
        # x/y/z occupy cols 31-54 (0-based [30:54]) regardless of resName width
        self.assertEqual(line[30:38], '  12.345')
        self.assertEqual(line[38:46], '  -6.789')
        self.assertEqual(line[46:54], '   0.500')

    def test_record_keyword_tracks_class(self):
        atom_line = write_atoms_pdb([_atom()], end=False)[0]
        het_line = write_atoms_pdb([_atom(cls=Hetatm)], end=False)[0]
        self.assertTrue(atom_line.startswith('ATOM  '))
        self.assertTrue(het_line.startswith('HETATM'))

    def test_wide_charmm_resname_does_not_overflow_coords(self):
        # 6-char CHARMM carb name would shift coords in the standard 4-col resName field
        line = write_atoms_pdb([_atom(cls=Hetatm, resname='BGLCNA', segname='AG01',
                                      x=218.343, y=223.129, z=159.086)], end=False)[0]
        self.assertEqual(line[17:23], 'BGLCNA')          # wide resName cols 18-23
        self.assertEqual(line[30:38], ' 218.343')         # coords still pinned
        self.assertEqual(line[38:46], ' 223.129')
        self.assertEqual(line[46:54], ' 159.086')
        self.assertEqual(line[72:76], 'AG01')             # segID cols 73-76

    def test_segid_falls_back_to_chain(self):
        line = write_atoms_pdb([_atom(segname=None, chainID='Q')], end=False)[0]
        self.assertEqual(line[72:76].strip(), 'Q')

    def test_insertion_code_written(self):
        line = write_atoms_pdb([_atom(resid=ResID(52, 'A'))], end=False)[0]
        self.assertEqual(line[28], 'A')                    # iCode col 29

    def test_file_written_with_end(self):
        fd, path = tempfile.mkstemp(suffix='.pdb')
        os.close(fd)
        try:
            write_atoms_pdb([_atom(), _atom(cls=Hetatm, serial=2)], filename=path)
            with open(path) as fh:
                lines = fh.read().splitlines()
            self.assertEqual(lines[-1], 'END')
            self.assertEqual(len(lines), 3)
        finally:
            os.remove(path)

    def test_atomlist_write_pdb_roundtrips_through_charmm_parse(self):
        al = AtomList([
            _atom(serial=1, name='N', elem='N', resname='ALA', chainID='B', resid=ResID(7),
                  x=10.0, y=20.0, z=30.0, segname='PROB'),
            _atom(cls=Hetatm, serial=2, name='C1', resname='ANE5AC', chainID='B', resid=ResID(8),
                  x=-1.5, y=2.25, z=99.999, segname='BG01'),
        ])
        fd, path = tempfile.mkstemp(suffix='.pdb')
        os.close(fd)
        try:
            al.write_pdb(path)
            parsed = PDBParser(filepath=path, dialect='charmm').parse().parsed
            back = AtomList.from_pdb(parsed)
            self.assertEqual(len(back), 2)
            a0, a1 = back[0], back[1]
            self.assertEqual((a0.resname, a0.chainID, a0.resid.resseqnum), ('ALA', 'B', 7))
            self.assertAlmostEqual(a0.x, 10.0, places=3)
            self.assertEqual(a1.resname, 'ANE5AC')          # 6-char name survives
            self.assertAlmostEqual(a1.z, 99.999, places=3)
        finally:
            os.remove(path)


if __name__ == '__main__':
    unittest.main()
