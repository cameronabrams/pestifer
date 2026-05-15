# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Unit tests for :mod:`pestifer.charmmff.ligand_paramgen.mol2_writer`.

The mol2 writer shells out to ``obabel``; tests are skipped when neither
RDKit nor the obabel binary are available.
"""
import shutil
import tempfile
import unittest
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    _HAVE_RDKIT = True
except ImportError:
    _HAVE_RDKIT = False

_HAVE_OBABEL = shutil.which("obabel") is not None


def _read_mol2_atoms(path: Path) -> list[tuple[str, str]]:
    """Return [(atom_name, atom_type)] from a Tripos mol2 file."""
    out = []
    in_atoms = False
    for line in path.read_text().splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if line.startswith("@<TRIPOS>"):
            in_atoms = False
            continue
        if in_atoms and line.strip():
            toks = line.split()
            # Cols: serial name x y z type subst_id subst_name charge
            out.append((toks[1], toks[5]))
    return out


def _read_mol2_title(path: Path) -> str:
    """Return the molecule-name line from a Tripos mol2 file."""
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("@<TRIPOS>MOLECULE"):
            return lines[i + 1].strip()
    raise AssertionError("MOLECULE record not found")


@unittest.skipUnless(_HAVE_RDKIT and _HAVE_OBABEL,
                     "rdkit + obabel binary required")
class TestWriteMol2(unittest.TestCase):

    def _build_named_mol(self, smiles: str, heavy_names: list[str]):
        """Build a 3D RDKit Mol with each heavy atom carrying a ``_pdb_name``."""
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        # Tag heavy atoms with the supplied names
        heavy = [a for a in mol.GetAtoms() if a.GetSymbol() != "H"]
        self.assertEqual(len(heavy), len(heavy_names),
                         "test set-up: heavy_names length mismatch")
        for atom, name in zip(heavy, heavy_names):
            atom.SetProp("_pdb_name", name)
        return mol

    def test_title_is_resname_not_tempfile_path(self):
        """Regression: the @<TRIPOS>MOLECULE line was leaking the tempfile
        path, which CGenFF then picked up as the residue name (e.g.
        ``RESI /tmp/AP5_xyz.pdb``)."""
        from pestifer.charmmff.ligand_paramgen.mol2_writer import write_mol2

        mol = self._build_named_mol("CC(=O)O", ["CA", "C", "O", "OXT"])
        with tempfile.TemporaryDirectory() as td:
            outpath = Path(td) / "ACE.mol2"
            write_mol2(mol, "ACE", outpath)
            title = _read_mol2_title(outpath)
            self.assertEqual(title, "ACE")
            self.assertNotIn("/tmp/", title)
            self.assertNotIn(".pdb", title)

    def test_heavy_atom_names_propagate(self):
        from pestifer.charmmff.ligand_paramgen.mol2_writer import write_mol2

        heavy_names = ["CA", "C", "O", "OXT"]
        mol = self._build_named_mol("CC(=O)O", heavy_names)
        with tempfile.TemporaryDirectory() as td:
            outpath = Path(td) / "ACE.mol2"
            write_mol2(mol, "ACE", outpath)
            atoms = _read_mol2_atoms(outpath)
            # Heavy atoms come first in PDB-block order
            heavy_atoms_in_mol2 = [n for (n, t) in atoms if not t.startswith("H")]
            self.assertEqual(heavy_atoms_in_mol2, heavy_names)

    def test_hydrogens_get_sequential_names(self):
        from pestifer.charmmff.ligand_paramgen.mol2_writer import write_mol2

        mol = self._build_named_mol("CN", ["N1", "C1"])
        with tempfile.TemporaryDirectory() as td:
            outpath = Path(td) / "TST.mol2"
            write_mol2(mol, "TST", outpath)
            atoms = _read_mol2_atoms(outpath)
            h_names = [n for (n, t) in atoms if t.startswith("H")]
            # Should be H1, H2, H3, ... in order
            self.assertEqual(h_names, [f"H{i}" for i in range(1, len(h_names) + 1)])

    def test_no_conformer_raises(self):
        from pestifer.charmmff.ligand_paramgen.mol2_writer import (
            Mol2WriteError, write_mol2,
        )
        mol = Chem.MolFromSmiles("CC")
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaises(Mol2WriteError):
                write_mol2(mol, "ETH", Path(td) / "out.mol2")


if __name__ == "__main__":
    unittest.main()
