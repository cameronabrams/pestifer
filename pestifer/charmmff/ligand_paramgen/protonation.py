# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Protonate a ligand residue at a target pH while preserving its PDB
heavy-atom order, names, and coordinates.

Used as the first stage of the on-the-fly CGenFF parameterization
pipeline: the protonated 3D Mol returned here feeds the mol2 writer
that drives the ``cgenff`` binary.
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import dimorphite_dl
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

if TYPE_CHECKING:
    from ...molecule.residue import Residue

logger = logging.getLogger(__name__)


class LigandProtonationError(RuntimeError):
    """Raised when a ligand cannot be protonated for CGenFF parameterization."""


def protonate_ligand(
    residue: "Residue",
    smiles: str,
    ph: float = 7.4,
) -> Chem.Mol:
    """
    Protonate a ligand residue at the given pH using Dimorphite-DL.

    The input ``residue`` supplies heavy-atom coordinates and PDB atom names;
    the ``smiles`` supplies bond orders and formal charges. Heavy-atom order
    in the returned Mol matches the input residue's heavy-atom order, and
    each heavy atom carries its PDB name as the RDKit property ``_pdb_name``.
    Hydrogens are placed in 3D at the protonation state appropriate for ``ph``.

    Parameters
    ----------
    residue
        A pestifer :class:`~pestifer.molecule.residue.Residue` for the ligand.
        Hydrogens in ``residue.atoms`` are ignored; only heavy atoms are used.
    smiles
        SMILES describing the ligand's heavy-atom skeleton and bond orders.
        Heavy-atom count must match the number of heavy atoms in ``residue``.
    ph
        pH at which to protonate (default 7.4). ``ph_min`` and ``ph_max``
        passed to Dimorphite-DL are both set to this value.
    """
    heavy_pdb_atoms = [
        a for a in residue.atoms.data
        if (a.elem or "").strip().upper() != "H"
    ]
    if not heavy_pdb_atoms:
        raise LigandProtonationError(f"Residue {residue} has no heavy atoms.")

    coord_mol = _build_coord_mol(heavy_pdb_atoms)
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise LigandProtonationError(f"Could not parse SMILES: {smiles!r}")
    if template.GetNumHeavyAtoms() != len(heavy_pdb_atoms):
        raise LigandProtonationError(
            f"SMILES has {template.GetNumHeavyAtoms()} heavy atoms but residue "
            f"{residue} has {len(heavy_pdb_atoms)}."
        )

    # Substructure match maps each template atom -> its position in coord_mol
    # (which is PDB heavy-atom order). Stamp the template's atoms with 1-indexed
    # map numbers that encode the PDB order; this survives the SMILES round-trip
    # through Dimorphite-DL.
    # coord_mol has all-single bonds from DetermineConnectivity, so we relax
    # bond-order strictness on the template's match.
    match = coord_mol.GetSubstructMatch(template)
    if not match:
        params = Chem.AdjustQueryParameters.NoAdjustments()
        params.makeBondsGeneric = True
        relaxed = Chem.AdjustQueryProperties(template, params)
        match = coord_mol.GetSubstructMatch(relaxed)
    if not match:
        raise LigandProtonationError(
            f"SMILES skeleton does not match {residue}'s heavy-atom graph."
        )
    for t_idx, m_idx in enumerate(match):
        template.GetAtomWithIdx(t_idx).SetAtomMapNum(m_idx + 1)

    mapped_smiles = Chem.MolToSmiles(template, canonical=False)
    logger.debug("Dimorphite input SMILES: %s", mapped_smiles)
    # precision=0.0 collapses Dimorphite-DL's pKa uncertainty window so it
    # returns a single dominant microspecies at the target pH; with the default
    # precision Dimorphite-DL widens to pKa +/- 1 and can return both forms
    # near borderline groups (e.g. amines whose pKa is ~3 above pH).
    protonated_smiles = dimorphite_dl.protonate_smiles(
        mapped_smiles,
        ph_min=ph,
        ph_max=ph,
        precision=0.0,
        max_variants=1,
    )
    if not protonated_smiles:
        raise LigandProtonationError(
            f"Dimorphite-DL returned no protomers for {residue} at pH {ph}."
        )
    logger.debug("Dimorphite output SMILES: %s", protonated_smiles[0])

    protonated = Chem.MolFromSmiles(protonated_smiles[0])
    if protonated is None:
        raise LigandProtonationError(
            f"Could not parse Dimorphite output SMILES: {protonated_smiles[0]!r}"
        )

    protonated = Chem.RenumberAtoms(protonated, _heavy_atom_map_order(protonated))

    n_heavy = len(heavy_pdb_atoms)
    if protonated.GetNumAtoms() < n_heavy:
        raise LigandProtonationError(
            f"Protonated SMILES has {protonated.GetNumAtoms()} atoms but expected "
            f"at least {n_heavy} heavy atoms."
        )

    # Copy heavy-atom 3D coords from coord_mol (which is in PDB order); stamp
    # PDB names on the heavy atoms.
    conf_src = coord_mol.GetConformer()
    conf_dst = Chem.Conformer(protonated.GetNumAtoms())
    for i in range(n_heavy):
        conf_dst.SetAtomPosition(i, conf_src.GetAtomPosition(i))
        protonated.GetAtomWithIdx(i).SetProp("_pdb_name", heavy_pdb_atoms[i].name)
        protonated.GetAtomWithIdx(i).SetAtomMapNum(0)
    protonated.RemoveAllConformers()
    protonated.AddConformer(conf_dst, assignId=True)

    return Chem.AddHs(protonated, addCoords=True)


def _build_coord_mol(heavy_atoms) -> Chem.Mol:
    """Build a single-conformer RDKit Mol from heavy atoms in input order.

    Connectivity is inferred from interatomic distances by
    :func:`rdDetermineBonds.DetermineConnectivity`; all bonds are single
    and bond orders are not meaningful. This Mol is used solely as the
    target of a (bond-order-relaxed) substructure match against the
    SMILES template, so we can stamp PDB-order atom map numbers onto
    the template's atoms.
    """
    rwmol = Chem.RWMol()
    conf = Chem.Conformer(len(heavy_atoms))
    for i, a in enumerate(heavy_atoms):
        sym = (a.elem or "").strip()
        if not sym:
            raise LigandProtonationError(
                f"Atom {a.name} in residue has no element symbol; "
                f"cannot build RDKit Mol."
            )
        idx = rwmol.AddAtom(Chem.Atom(sym))
        conf.SetAtomPosition(idx, (float(a.x), float(a.y), float(a.z)))
    rwmol.AddConformer(conf, assignId=True)
    mol = rwmol.GetMol()
    try:
        rdDetermineBonds.DetermineConnectivity(mol)
    except ValueError as exc:
        raise LigandProtonationError(
            f"RDKit could not infer connectivity from heavy-atom coordinates: {exc}"
        ) from exc
    return mol


def _heavy_atom_map_order(mol: Chem.Mol) -> list[int]:
    """Permutation that reorders ``mol`` so mapped atoms come first in
    map-number order, with any unmapped atoms (e.g., explicit Hs) trailing."""
    mapped = sorted(
        (a.GetAtomMapNum(), a.GetIdx())
        for a in mol.GetAtoms()
        if a.GetAtomMapNum() > 0
    )
    unmapped = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum() == 0]
    return [idx for _, idx in mapped] + unmapped
