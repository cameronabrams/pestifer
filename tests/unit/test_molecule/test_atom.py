import unittest

from pestifer.molecule.atom import Atom, AtomList, Hetatm
from pidibble.pdbparse import PDBParser
from pidibble.pdbrecord import PDBRecordDict
from pestifer.util.cifutil import CIFdict, CIFload
from pestifer.psfutil.psfatom import PSFAtom, PSFAtomList
from pathlib import Path
from pestifer.objs.resid import ResID

class TestAtom(unittest.TestCase):

    inputs_dir = Path(__file__).parents[2] / "inputs"

    def test_atom_create(self):
        atom = Atom(
            serial=1,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
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
        self.assertEqual(atom.resid.resid, 1)
        self.assertEqual(atom.resid.pdbresid, '1 ')
        # check default values
        self.assertEqual(atom.altloc, ' ')
        self.assertEqual(atom.segname, None)
        self.assertEqual(atom.empty, False)
        self.assertEqual(atom.link, 'None')
        self.assertEqual(atom.recordname, 'ATOM')
        self.assertEqual(atom.auth_seq_id, None)
        self.assertEqual(atom.auth_comp_id, None)
        self.assertEqual(atom.auth_asym_id, None)
        self.assertEqual(atom.auth_atom_id, None)

        self.assertEqual(atom._PDB_keyword, 'ATOM')
        self.assertEqual(atom._CIF_CategoryName, 'atom_site')
        self.assertEqual(atom.ORIGINAL_ATTRIBUTES, {})

    def test_atom_copy(self):
        atom1 = Atom(
            serial=1,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
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
        p=PDBParser(filepath=self.inputs_dir / '4zmj.pdb').parse()
        atom_record = p.parsed[Atom._PDB_keyword][0]
        atom = Atom(atom_record)
        self.assertEqual(atom.serial, atom_record.serial)
        self.assertEqual(atom.name, atom_record.name)
        self.assertEqual(atom.resname, atom_record.residue.resName)
        self.assertEqual(atom.chainID, atom_record.residue.chainID)
        self.assertEqual(atom.resid.resseqnum, atom_record.residue.seqNum)
        self.assertEqual(atom.resid.insertion, atom_record.residue.iCode)
        self.assertEqual(atom.x, atom_record.x)
        self.assertEqual(atom.y, atom_record.y)
        self.assertEqual(atom.z, atom_record.z)
        self.assertEqual(atom.occ, atom_record.occupancy)
        self.assertEqual(atom.beta, atom_record.tempFactor)
        self.assertEqual(atom.elem, atom_record.element)
        self.assertEqual(atom.charge, atom_record.charge)

    def test_atom_from_cifdict(self):
        p=CIFload(Path(self.inputs_dir / '4zmj.cif'))
        obj=p.getObj(Atom._CIF_CategoryName)
        d=CIFdict(obj, 0)
        atom = Atom(d)
        self.assertEqual(atom.serial, int(d['id']))

    def test_atom_pdb_line(self):
        atom = Atom(
            serial=1,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=1.0,
            y=2.0,
            z=3.0,
            occ=1.0,
            beta=0.0,
            elem='C',
            charge=' '
        )
        pdbline = atom.pdb_line()
        expected_line = 'ATOM      1  C   ALA A   1       1.000   2.000   3.000  1.00  0.00           C  '
        self.assertEqual(pdbline, expected_line)

    def test_atom_overwrite_position_from_Atom(self):
        atom = Atom(
            serial=1,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=1.0,
            y=2.0,
            z=3.0,
            occ=1.0,
            beta=0.0,
            elem='C',
            charge=' '
        )
        other_atom = Atom(
            serial=2,
            name='O',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=2.0,
            y=3.0,
            z=4.0,
            occ=1.0,
            beta=0.0,
            elem='O',
            charge=' '
        )
        atom.overwrite_position(other_atom)
        self.assertEqual(atom.x, other_atom.x)
        self.assertEqual(atom.y, other_atom.y)
        self.assertEqual(atom.z, other_atom.z)

    def test_atom_overwrite_position_from_dict(self):
        atom = Atom(
            serial=1,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=1.0,
            y=2.0,
            z=3.0,
            occ=1.0,
            beta=0.0,
            elem='C',
            charge=' '
        )
        position_dict = {'x': 2.0, 'y': 3.0, 'z': 4.0}
        atom.overwrite_position(position_dict)
        self.assertEqual(atom.x, position_dict['x'])
        self.assertEqual(atom.y, position_dict['y'])
        self.assertEqual(atom.z, position_dict['z'])

    def test_atom_overwrite_position_from_floats(self):
        atom = Atom(
            serial=1,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=1.0,
            y=2.0,
            z=3.0,
            occ=1.0,
            beta=0.0,
            elem='C',
            charge=' '
        )
        atom.overwrite_position(2.0, 3.0, 4.0)
        self.assertEqual(atom.x, 2.0)
        self.assertEqual(atom.y, 3.0)
        self.assertEqual(atom.z, 4.0)

class TestAtomList(unittest.TestCase):

    inputs_dir = Path(__file__).parents[2] / "inputs"

    def tearDown(self):
        for f in Path('.').glob("*.log"):
            f.unlink()

    def test_atom_list_create(self):
        atom1 = Atom(
            serial=1,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=1.0,
            y=2.0,
            z=3.0,
            occ=1.0,
            beta=0.0,
            elem='C',
            charge=' '
        )
        atom2 = Atom(
            serial=2,
            name='O',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=2.0,
            y=3.0,
            z=4.0,
            occ=1.0,
            beta=0.0,
            elem='O',
            charge=' '
        )
        atom_list = AtomList([atom1, atom2])
        self.assertEqual(len(atom_list), 2)
        self.assertEqual(atom_list[0].serial, 1)
        self.assertEqual(atom_list[1].serial, 2)

    def test_atom_list_reserialize(self):
        atom1 = Atom(
            serial=9,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=1.0,
            y=2.0,
            z=3.0,
            occ=1.0,
            beta=0.0,
            elem='C',
            charge=' '
        )
        atom2 = Atom(
            serial=10,
            name='O',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=2.0,
            y=3.0,
            z=4.0,
            occ=1.0,
            beta=0.0,
            elem='O',
            charge=' '
        )
        atom_list = AtomList([atom1, atom2])
        self.assertEqual(atom_list[0].serial, 9)
        self.assertEqual(atom_list[1].serial, 10)
        atom_list.reserialize()
        self.assertEqual(atom_list[0].serial, 1)
        self.assertEqual(atom_list[1].serial, 2)
        self.assertEqual(atom_list[0].ORIGINAL_ATTRIBUTES['serial'], 9)
        self.assertEqual(atom_list[1].ORIGINAL_ATTRIBUTES['serial'], 10)
        self.assertFalse(atom_list[0].ORIGINAL_ATTRIBUTES is atom_list[1].ORIGINAL_ATTRIBUTES)

    def test_atom_list_from_pdb(self):
        p = PDBParser(filepath=self.inputs_dir / '4zmj.pdb').parse().parsed
        self.assertIsInstance(p, PDBRecordDict)
        atom_list = AtomList.from_pdb(p)
        self.assertEqual(len(atom_list), 4856)
        self.assertIsInstance(atom_list[0], Atom)
        self.assertIsInstance(atom_list[4518], Hetatm)
        self.assertEqual(atom_list[0].serial, 1)
        self.assertEqual(atom_list[0].name, 'N')
        self.assertEqual(atom_list[0].resname, 'LEU')
        self.assertEqual(atom_list[0].chainID, 'G')
        self.assertEqual(atom_list[0].x, -0.092)
        self.assertEqual(atom_list[0].y, 99.33)
        self.assertEqual(atom_list[0].z, 57.967)

    def test_atom_list_from_cif(self):
        p = CIFload(Path(self.inputs_dir / '4zmj.cif'))
        atom_list = AtomList.from_cif(p)
        self.assertEqual(len(atom_list), 4856)
        self.assertIsInstance(atom_list[0], Atom)
        self.assertEqual(atom_list[0].serial, 1)
        self.assertEqual(atom_list[0].name, 'N')
        self.assertEqual(atom_list[0].resname, 'LEU')
        self.assertEqual(atom_list[0].chainID, 'A')
        self.assertEqual(atom_list[0].auth_asym_id, 'G')
        self.assertEqual(atom_list[0].x, -0.092)
        self.assertEqual(atom_list[0].y, 99.33)
        self.assertEqual(atom_list[0].z, 57.967)

    def test_atom_list_overwrite_positions(self):
        atom1 = Atom(
            serial=1,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=1.0,
            y=2.0,
            z=3.0,
            occ=1.0,
            beta=0.0,
            elem='C',
            charge=' '
        )
        atom2 = Atom(
            serial=2,
            name='O',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=2.0,
            y=3.0,
            z=4.0,
            occ=1.0,
            beta=0.0,
            elem='O',
            charge=' '
        )
        atom_list1 = AtomList([atom1])
        atom_list2 = AtomList([atom2])
        atom_list1.overwrite_positions(atom_list2)
        self.assertEqual(atom_list1[0].x, 2.0)
        self.assertEqual(atom_list1[0].y, 3.0)
        self.assertEqual(atom_list1[0].z, 4.0)

    def test_atom_list_apply_psf_resnames(self):
        atom1 = Atom(
            serial=1,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=1.0,
            y=2.0,
            z=3.0,
            occ=1.0,
            beta=0.0,
            elem='C',
            charge=' '
        )
        atom2 = Atom(
            serial=2,
            name='O',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=2.0,
            y=3.0,
            z=4.0,
            occ=1.0,
            beta=0.0,
            elem='O',
            charge=' '
        )
        atom_list = AtomList([atom1, atom2])
        psf_atom1 = PSFAtom('1 C 3A ARG C C 0.00 12.01 dum')
        psf_atom2 = PSFAtom('2 C 4  ARG O O 0.00 16.00 dum')
        psf_atom_list = PSFAtomList([psf_atom1, psf_atom2])
        atom_list.apply_psf_attributes(psf_atom_list)
        self.assertEqual(atom_list[0].resname, 'ARG')
        self.assertEqual(atom_list[1].resname, 'ARG')