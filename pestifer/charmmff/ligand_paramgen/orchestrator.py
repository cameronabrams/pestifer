# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Orchestrates the per-PDB pipeline for the ``make-ligand-mol2`` subcommand:
fetch or load the PDB, enumerate HETATM residues, drop any whose resname
is already covered by the CHARMM force field, and for the remaining
ligands run SMILES lookup -> protonation -> mol2 emission.

Same internals will be reused by the build workflow when it needs to
either drive a local ``cgenff`` binary or pause for the user to upload
the mol2 to the CGenFF web tool.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from pidibble.pdbparse import PDBParser

from ...molecule.atom import AtomList, Hetatm
from ...molecule.residue import Residue
from .mol2_writer import Mol2WriteError, write_mol2
from .protonation import LigandProtonationError, protonate_ligand
from .rcsb import RCSBLookupError, fetch_ligand_smiles

logger = logging.getLogger(__name__)


@dataclass
class LigandResult:
    resname: str
    status: str  # 'ok' | 'no_smiles' | 'protonation_failed' | 'mol2_failed'
    path: Path | None = None
    smiles: str | None = None
    message: str = ""


@dataclass
class LigandGenSummary:
    successes: list[LigandResult] = field(default_factory=list)
    failures: list[LigandResult] = field(default_factory=list)
    skipped_known: list[str] = field(default_factory=list)


def fetch_or_load_pdb(source: str):
    """
    Resolve ``source`` to a parsed pidibble :class:`PDBRecordDict`.

    A 4-character alphanumeric token is treated as an RCSB PDB code and
    fetched; anything else is treated as a path to a local ``.pdb`` file.
    """
    src = source.strip()
    path = Path(src)
    if path.is_file():
        parser = PDBParser(filepath=path)
    elif len(src) == 4 and src.isalnum():
        parser = PDBParser(source_db="rcsb", source_id=src.upper())
        parser.fetch()
    else:
        raise ValueError(
            f"Could not interpret {source!r} as a 4-letter PDB code or a "
            f"local .pdb file."
        )
    parser.parse()
    return parser.parsed


def find_unknown_ligand_residues(
    parsed,
    charmmff_content,
) -> tuple[dict[str, Residue], list[str]]:
    """
    Group HETATM records by resname; for each resname not covered by the
    CHARMM force field, return a representative pestifer :class:`Residue`
    built from the atoms of its first appearance in the structure.

    Returns
    -------
    (unknown_by_resname, skipped_known)
        ``unknown_by_resname`` maps resname -> a single Residue suitable
        for :func:`protonate_ligand`. ``skipped_known`` lists resnames
        that were filtered out because CHARMM already defines them.
    """
    records = parsed.get(Hetatm._PDB_keyword, [])
    if not records:
        return {}, []

    atoms_by_resname: dict[str, list[Hetatm]] = {}
    for rec in records:
        atom = Hetatm(rec)
        atoms_by_resname.setdefault(atom.resname, []).append(atom)

    unknowns: dict[str, Residue] = {}
    skipped: list[str] = []
    for resname, atoms in atoms_by_resname.items():
        if charmmff_content.get_topfile_of_resname(resname) is not None:
            skipped.append(resname)
            continue
        # Pick the first (chainID, resid) group as the representative
        # residue, in case the same resname appears multiple times.
        first = atoms[0]
        key = (first.chainID, first.resid.resid)
        residue_atoms = [
            a for a in atoms if (a.chainID, a.resid.resid) == key
        ]
        residue_atoms = _drop_altlocs(residue_atoms)
        unknowns[resname] = _residue_from_hetatms(residue_atoms)
    return unknowns, skipped


def generate_ligand_mol2s(
    parsed,
    charmmff_content,
    outdir: Path | str,
    ph: float = 7.4,
    smiles_overrides: dict[str, str] | None = None,
) -> LigandGenSummary:
    """
    For each unknown ligand resname in the parsed PDB, fetch a SMILES,
    protonate at ``ph``, and write a mol2 file under ``outdir``.

    On per-ligand failure the function logs a warning and continues; the
    returned summary lists what succeeded and what failed and why.
    """
    overrides = smiles_overrides or {}
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary = LigandGenSummary()

    unknowns, skipped = find_unknown_ligand_residues(parsed, charmmff_content)
    summary.skipped_known = sorted(skipped)
    if not unknowns:
        logger.info("No unknown ligands found.")
        return summary

    for resname in sorted(unknowns):
        residue = unknowns[resname]
        smi = overrides.get(resname)
        if smi is None:
            try:
                smi = fetch_ligand_smiles(resname)
            except RCSBLookupError as exc:
                logger.warning("%s: no SMILES (%s)", resname, exc)
                summary.failures.append(LigandResult(
                    resname, "no_smiles", message=str(exc)
                ))
                continue
        try:
            mol = protonate_ligand(residue, smi, ph=ph)
        except LigandProtonationError as exc:
            logger.warning("%s: protonation failed (%s)", resname, exc)
            summary.failures.append(LigandResult(
                resname, "protonation_failed", smiles=smi, message=str(exc)
            ))
            continue
        outpath = outdir / f"{resname}.mol2"
        try:
            write_mol2(mol, resname, outpath)
        except Mol2WriteError as exc:
            logger.warning("%s: mol2 write failed (%s)", resname, exc)
            summary.failures.append(LigandResult(
                resname, "mol2_failed", smiles=smi, message=str(exc)
            ))
            continue
        logger.info("%s: wrote %s", resname, outpath)
        summary.successes.append(LigandResult(
            resname, "ok", path=outpath, smiles=smi
        ))
    return summary


def _drop_altlocs(atoms: list[Hetatm]) -> list[Hetatm]:
    """Deduplicate alternate locations.

    For each unique atom name, keep one record: the one with the highest
    occupancy. Ties are broken by preferring no altloc, then ``A``, then
    ``B``, etc. (the standard PDB convention).

    Some PDB entries (e.g., 1AKE's AP5) carry disordered ligand atoms as
    multiple HETATM records with the same atom name but different altlocs;
    feeding all of them to the protonation step would make the heavy-atom
    count exceed the SMILES skeleton.
    """
    def priority(a: Hetatm) -> tuple[float, int]:
        # Higher tuple wins. occ first; then altloc (blank highest, then A>B>...).
        altloc = (a.altloc or "").strip()
        altloc_rank = 999 if not altloc else -ord(altloc[:1])
        return (a.occ, altloc_rank)

    best: dict[str, Hetatm] = {}
    for a in atoms:
        name = a.name.strip()
        prev = best.get(name)
        if prev is None or priority(a) > priority(prev):
            best[name] = a
    chosen = set(map(id, best.values()))
    return [a for a in atoms if id(a) in chosen]


def _residue_from_hetatms(atoms: Iterable[Hetatm]) -> Residue:
    """Build a minimal pestifer :class:`Residue` from a list of HETATM
    atoms, bypassing :data:`Labels.segtype_of_resname` (which would
    KeyError for novel resnames)."""
    atom_list = list(atoms)
    first = atom_list[0]
    return Residue(
        resname=first.resname,
        resid=first.resid,
        chainID=first.chainID,
        asym_chainID=first.chainID,
        segtype="ligand",
        segname=first.chainID,
        auth_asym_id=getattr(first, "auth_asym_id", None),
        auth_comp_id=getattr(first, "auth_comp_id", None),
        auth_seq_id=getattr(first, "auth_seq_id", None),
        atoms=AtomList(atom_list),
        resolved=True,
        recordname="HETATM",
    )
