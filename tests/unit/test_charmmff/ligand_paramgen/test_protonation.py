# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Unit tests for :mod:`pestifer.charmmff.ligand_paramgen.protonation`.

These tests skip cleanly when the optional dependencies (``rdkit``,
``dimorphite_dl``) are not installed.
"""
import unittest

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import dimorphite_dl  # noqa: F401
    _HAVE_DEPS = True
except ImportError:
    _HAVE_DEPS = False


@unittest.skipUnless(_HAVE_DEPS, "rdkit + dimorphite_dl required")
class TestProtonateLigand(unittest.TestCase):
    """Round-trip a synthetic residue through ``protonate_ligand``.

    For each case we embed the SMILES to 3D, manufacture a pestifer Residue
    whose heavy atoms carry PDB-style names, and feed it back through the
    function. We verify the invariants the downstream CGenFF stages rely
    on:

    1. heavy-atom order matches the input PDB residue
    2. PDB atom names are stamped on every heavy atom via ``_pdb_name``
    3. Hs were placed at non-origin 3D positions
    4. total formal charge falls in a chemically plausible range for pH 7.4
    """

    def _build_residue(self, smiles: str, resname: str):
        """Embed ``smiles`` and wrap the heavy atoms as a pestifer Residue."""
        from pestifer.molecule.atom import AtomList, Hetatm
        from pestifer.molecule.residue import Residue
        from pestifer.objs.resid import ResID

        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        mol = Chem.RemoveHs(mol)

        atoms = AtomList()
        per_elem: dict[str, int] = {}
        conf = mol.GetConformer()
        for i, a in enumerate(mol.GetAtoms()):
            e = a.GetSymbol()
            per_elem[e] = per_elem.get(e, 0) + 1
            pos = conf.GetAtomPosition(i)
            atoms.append(Hetatm(
                serial=i + 1, name=f"{e}{per_elem[e]}", altloc=" ",
                resname=resname, chainID="A",
                resid=ResID(resseqnum=1, insertion=" "),
                x=pos.x, y=pos.y, z=pos.z, occ=1.0, beta=0.0,
                elem=e, charge="0",
            ))
        return Residue(
            resname=resname, resid=ResID(resseqnum=1, insertion=" "),
            chainID="A", asym_chainID="A", segtype="ligand", segname="A",
            auth_asym_id="A", auth_comp_id=resname, auth_seq_id=1,
            atoms=atoms, resolved=True, recordname="HETATM",
        )

    def _assert_invariants(
        self, residue, mol, charge_min: int, charge_max: int
    ):
        from rdkit import Chem  # local for type stability under skip

        n_heavy_in = len(residue.atoms.data)
        self.assertEqual(mol.GetNumHeavyAtoms(), n_heavy_in,
                         "heavy-atom count must match input")
        expected_names = [a.name for a in residue.atoms.data]
        actual_names = [mol.GetAtomWithIdx(i).GetProp("_pdb_name")
                        for i in range(n_heavy_in)]
        self.assertEqual(actual_names, expected_names,
                         "PDB names must round-trip in input order")

        # Hs must have real 3D coordinates (AddHs(addCoords=True))
        n_h = mol.GetNumAtoms() - n_heavy_in
        if n_h > 0:
            conf = mol.GetConformer()
            for atom in mol.GetAtoms():
                if atom.GetSymbol() != "H":
                    continue
                pos = conf.GetAtomPosition(atom.GetIdx())
                self.assertGreater(
                    pos.x * pos.x + pos.y * pos.y + pos.z * pos.z, 0.0,
                    f"H {atom.GetIdx()} sits at origin"
                )

        total_q = sum(a.GetFormalCharge() for a in mol.GetAtoms())
        self.assertTrue(
            charge_min <= total_q <= charge_max,
            f"charge {total_q} outside [{charge_min}, {charge_max}]"
        )

    def test_methylamine_protonates(self):
        from pestifer.charmmff.ligand_paramgen import protonate_ligand
        residue = self._build_residue("CN", "MAM")
        mol = protonate_ligand(residue, "CN", ph=7.4)
        self._assert_invariants(residue, mol, +1, +1)

    def test_acetic_acid_deprotonates(self):
        from pestifer.charmmff.ligand_paramgen import protonate_ligand
        residue = self._build_residue("CC(=O)O", "ACE")
        mol = protonate_ligand(residue, "CC(=O)O", ph=7.4)
        self._assert_invariants(residue, mol, -1, -1)

    def test_glycine_zwitterion(self):
        from pestifer.charmmff.ligand_paramgen import protonate_ligand
        residue = self._build_residue("NCC(=O)O", "GLY")
        mol = protonate_ligand(residue, "NCC(=O)O", ph=7.4)
        # NH3+ ... COO- — net zero at pH 7.4
        self._assert_invariants(residue, mol, 0, 0)

    def test_phosphate_deprotonates(self):
        from pestifer.charmmff.ligand_paramgen import protonate_ligand
        residue = self._build_residue("OP(=O)(O)O", "PHO")
        mol = protonate_ligand(residue, "OP(=O)(O)O", ph=7.4)
        # pKa2 = 7.21; precision=0.0 picks dominant -- accept HPO4^2- or PO4^3-
        self._assert_invariants(residue, mol, -3, -2)

    def test_mismatched_heavy_atom_count_raises(self):
        from pestifer.charmmff.ligand_paramgen import (
            LigandProtonationError, protonate_ligand,
        )
        residue = self._build_residue("CN", "MAM")
        # Use a SMILES with the wrong number of heavy atoms
        with self.assertRaises(LigandProtonationError):
            protonate_ligand(residue, "CCO", ph=7.4)

    def test_unparseable_smiles_raises(self):
        from pestifer.charmmff.ligand_paramgen import (
            LigandProtonationError, protonate_ligand,
        )
        residue = self._build_residue("CN", "MAM")
        with self.assertRaises(LigandProtonationError):
            protonate_ligand(residue, "not-a-smiles!!!", ph=7.4)


if __name__ == "__main__":
    unittest.main()
