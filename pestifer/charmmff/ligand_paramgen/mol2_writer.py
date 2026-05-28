# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Write a protonated RDKit Mol to a Tripos mol2 file via Open Babel.

We stamp PDB-style atom names onto every atom (heavy atoms keep the names
attached by :func:`protonate_ligand`; hydrogens get sequential ``H<n>``
names), emit a PDB block, then shell out to ``obabel`` for format
conversion. Open Babel only sees the already-protonated, already-named
structure, so it never mangles atom names or alters the H count.
"""
from __future__ import annotations

import logging
import os
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

from ...core.command import Command

if TYPE_CHECKING:
    from rdkit import Chem

logger = logging.getLogger(__name__)


class Mol2WriteError(RuntimeError):
    """Raised when mol2 generation fails."""


def _import_rdkit_chem():
    """Lazily import ``rdkit.Chem`` (an optional ``ligand-paramgen`` dependency)."""
    try:
        from rdkit import Chem
    except ModuleNotFoundError as exc:
        raise Mol2WriteError(
            "RDKit is required for ligand parameterization. Install the "
            "optional dependencies with `pip install pestifer[ligand-paramgen]`."
        ) from exc
    return Chem


def write_mol2(
    mol: Chem.Mol,
    resname: str,
    outpath: str | os.PathLike,
) -> Path:
    """
    Write ``mol`` to ``outpath`` as a mol2 file via Open Babel.

    Atom names on the result are:

    - heavy atoms: the value of the RDKit property ``_pdb_name`` if set,
      otherwise the atom's element symbol.
    - hydrogens: sequential ``H1``, ``H2``, ... in RDKit atom-index order.

    Parameters
    ----------
    mol
        RDKit Mol with 3D coordinates and explicit Hs (as returned by
        :func:`pestifer.charmmff.ligand_paramgen.protonate_ligand`).
    resname
        Three-letter residue name to stamp on every atom (used by CGenFF
        and by Open Babel's PDB reader).
    outpath
        Destination mol2 path.

    Returns
    -------
    pathlib.Path
        The path that was written.
    """
    Chem = _import_rdkit_chem()
    if mol.GetNumConformers() == 0:
        raise Mol2WriteError("Mol has no conformer; cannot write 3D mol2.")

    _stamp_pdb_names(mol, resname)
    pdb_block = Chem.MolToPDBBlock(mol, flavor=0)

    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        suffix=".pdb", prefix=f"{resname}_", mode="w", delete=False,
    ) as tmp:
        tmp.write(pdb_block)
        tmp_pdb = tmp.name
    try:
        # `obabel <in> -O <out>` infers formats from extensions. `--title`
        # sets the @<TRIPOS>MOLECULE name; without it obabel writes the
        # input filename, which CGenFF then picks up as the RESI name —
        # so the resulting .str registers a residue named like
        # "/tmp/AP5_xyz.pdb" that pestifer can never match.
        c = Command(f"obabel {tmp_pdb} -O {outpath} --title {resname}")
        rc = c.run(quiet=True)
        if rc != 0 or not outpath.exists():
            raise Mol2WriteError(
                f"obabel failed (rc={rc}) for {resname}: {c.stderr or c.stdout}"
            )
    finally:
        try:
            os.unlink(tmp_pdb)
        except OSError:
            pass
    return outpath


def _stamp_pdb_names(mol: Chem.Mol, resname: str) -> None:
    """Attach AtomPDBResidueInfo to every atom so RDKit's PDB writer emits
    a stable name for each one."""
    Chem = _import_rdkit_chem()
    resname3 = resname.strip()[:3].ljust(3)
    h_idx = 0
    used_names: set[str] = set()
    for atom in mol.GetAtoms():
        if atom.HasProp("_pdb_name"):
            name = atom.GetProp("_pdb_name")
        elif atom.GetSymbol() == "H":
            h_idx += 1
            name = f"H{h_idx}"
        else:
            name = atom.GetSymbol()
        name = name.strip()[:4]
        # Avoid duplicate atom names within a residue (PDB requires uniqueness).
        if name in used_names:
            base, suffix = name, 1
            while f"{base}{suffix}" in used_names:
                suffix += 1
            name = f"{base}{suffix}"[:4]
        used_names.add(name)
        info = Chem.AtomPDBResidueInfo()
        info.SetName(name.ljust(4))
        info.SetResidueName(resname3)
        info.SetChainId("A")
        info.SetResidueNumber(1)
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
