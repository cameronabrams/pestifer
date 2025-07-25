import unittest

from pestifer.molecule.atom import Atom, AtomList
from pidibble.pdbparse import PDBParser
from pestifer.util.cifutil import CIFdict, CIFload
from pathlib import Path
class TestAtom(unittest.TestCase):
    def test_atom_create(self):
        atom = Atom(
            serial=1,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resseqnum=1,
            insertion=' ',
            x=1.0,
            y=2.0,
            z=3.0,
            occ=1.0,
            beta=0.0,
            elem='C',
            charge=' '
        )
        self.assertEqual(atom.serial, 1)
        self.assertEqual(atom.name, 'C')
        self.assertEqual(atom.resname, 'ALA')
        self.assertEqual(atom.chainID, 'A')
        self.assertEqual(atom.x, 1.0)
        self.assertEqual(atom.y, 2.0)
        self.assertEqual(atom.z, 3.0)
        self.assertEqual(atom.occ, 1.0)
        self.assertEqual(atom.beta, 0.0)
        self.assertEqual(atom.elem, 'C')
        self.assertEqual(atom.charge, ' ')
        # check default values
        self.assertEqual(atom.altloc, ' ')
        self.assertEqual(atom.insertion, ' ')
        self.assertEqual(atom.segname, None)
        self.assertEqual(atom.empty, False)
        self.assertEqual(atom.link, 'None')
        self.assertEqual(atom.recordname, 'ATOM')
        self.assertEqual(atom.auth_seq_id, None)
        self.assertEqual(atom.auth_comp_id, None)
        self.assertEqual(atom.auth_asym_id, None)
        self.assertEqual(atom.auth_atom_id, None)

    def test_atom_copy(self):
        atom1 = Atom(
            serial=1,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resseqnum=1,
            insertion=' ',
            x=1.0,
            y=2.0,
            z=3.0,
            occ=1.0,
            beta=0.0,
            elem='C',
            charge=' '
        )
        atom2 = atom1.copy()
        self.assertEqual(atom2.serial, atom1.serial)
        self.assertEqual(atom2.name, atom1.name)
        self.assertEqual(atom2.resname, atom1.resname)
        self.assertEqual(atom2.chainID, atom1.chainID)
        self.assertEqual(atom2.x, atom1.x)
        self.assertEqual(atom2.y, atom1.y)
        self.assertEqual(atom2.z, atom1.z)
        self.assertEqual(atom2.occ, atom1.occ)
        self.assertEqual(atom2.beta, atom1.beta)
        self.assertEqual(atom2.elem, atom1.elem)
        self.assertEqual(atom2.charge, atom1.charge)

    def test_atom_from_pdbrecord(self):
        p=PDBParser(filepath='fixtures/data/4zmj.pdb').parse()
        atom_record = p.parsed[Atom._PDB_keyword][0]
        atom = Atom.new(atom_record)
        self.assertEqual(atom.serial, atom_record.serial)
        self.assertEqual(atom.name, atom_record.name)
        self.assertEqual(atom.resname, atom_record.residue.resName)
        self.assertEqual(atom.chainID, atom_record.residue.chainID)
        self.assertEqual(atom.resseqnum, atom_record.residue.seqNum)
        self.assertEqual(atom.insertion, atom_record.residue.iCode)
        self.assertEqual(atom.x, atom_record.x)
        self.assertEqual(atom.y, atom_record.y)
        self.assertEqual(atom.z, atom_record.z)
        self.assertEqual(atom.occ, atom_record.occupancy)
        self.assertEqual(atom.beta, atom_record.tempFactor)
        self.assertEqual(atom.elem, atom_record.element)
        self.assertEqual(atom.charge, atom_record.charge)

    def test_atom_from_cifdict(self):
        p=CIFload(Path('fixtures/data/4zmj.cif'))
        obj=p.getObj(Atom._mmCIF_name)
        d=CIFdict(obj, 0)
        atom = Atom.new(d)
        self.assertEqual(atom.serial, int(d['id']))
