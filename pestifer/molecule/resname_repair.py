# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Restore CHARMM residue names that a PDB writer cut to four characters -- by what the residue
contains, not by what the stub is called.

Several writers keep only four columns of a residue name.  CHARMM's longer names collide when cut,
and the four-character stub is usually a real residue too: BGLCNA (beta-GlcNAc), BGLCA (beta-glucuronic
acid) and BGLCN all become BGLC, which is beta-glucose.  pestifer used to paper over one case with a
psfgen alias, ``BGLC -> BGLCNA``.  psfgen applies an alias to every residue of that name, so a genuine
beta-glucose -- rebuilding from a CHARMM-named PDB of example 31, for instance -- was silently built
as GlcNAc, with its acetyl atoms invented.

Here each residue is matched against the candidates in ``truncated_resnames.json`` (generated from the
force field): a residue keeps its name if its heavy atoms fit that residue, and otherwise takes the one
longer candidate they fit.
"""
from __future__ import annotations

import logging
from collections import OrderedDict

from ..core.labels import Labels, load_truncated_resname_table, _atom_aliases, _residue_aliases

logger = logging.getLogger(__name__)


def _is_heavy(atom) -> bool:
    elem = (getattr(atom, 'elem', '') or '').strip().upper()
    if elem:
        return elem != 'H'
    return not atom.name.upper().startswith('H')


def _pdb_to_charmm_atom_names(resname: str) -> dict:
    """PDB-style atom name -> CHARMM atom name, from ``resname``'s atom aliases (C7 -> C for BGLCNA)."""
    out = {}
    for alias in _atom_aliases:
        parts = alias.split()
        if len(parts) >= 3 and parts[0] == resname:
            out[parts[1]] = parts[2]
    return out


def choose_resname(stem: str, heavy: set, table: dict | None = None) -> tuple[str, str]:
    """
    Return ``(resname, reason)`` for a residue named ``stem`` with heavy atoms ``heavy``.
    ``resname`` is ``stem`` unless a longer CHARMM residue is the one unambiguous fit.

    Each candidate is judged on the residue's atom names translated into that candidate's own names
    (its atom aliases), against the candidate's real heavy atoms: the residue must be a subset, and
    among those that fit, fewer absent atoms is better.
    """
    table = load_truncated_resname_table() if table is None else table
    cands = table.get(stem)
    if not cands:
        return stem, 'not a truncation stem'
    missing = {}
    for resname, atoms in cands.items():
        names = set(atoms)
        rename = _pdb_to_charmm_atom_names(resname)
        translated = {rename.get(n, n) for n in heavy}
        if translated <= names:
            missing[resname] = len(names - translated)
    if stem in missing:
        return stem, 'fits its own residue'
    if not missing:
        return stem, 'fits no candidate'
    fewest = min(missing.values())
    best = sorted(r for r, m in missing.items() if m == fewest)
    if len(best) > 1:
        built = {a.split()[1] for a in _residue_aliases}      # what pestifer builds from PDB codes
        preferred = [r for r in best if r in built]
        if len(preferred) == 1:
            return preferred[0], f'tied with {[r for r in best if r not in preferred]}; preferred as the residue pestifer builds'
        return stem, f'ambiguous between {best}'
    return best[0], 'unique fit'


def _writable_name(charmm: str) -> str | None:
    """
    The name to give a residue so it reaches psfgen as ``charmm``.

    pestifer writes the segment PDBs psfgen reads in standard columns, so a residue name keeps at
    most four characters: a residue renamed BGLCNA would be written as BGLC and built as glucose.
    A name of four or fewer characters is written as is.  A longer one is given the PDB code that a
    residue alias maps to it (BGLCNA -> NAG), which is exactly how a deposit's own sugars reach psfgen.
    """
    if len(charmm) <= 4:
        return charmm
    codes = sorted((a.split()[0] for a in _residue_aliases if a.split()[1] == charmm), key=len)
    codes = [c for c in codes if len(c) <= 4]
    return codes[0] if codes else None


def repair_truncated_resnames(atoms, table: dict | None = None) -> int:
    """Rename, in place, the atoms of every residue whose four-character name hides a longer CHARMM
    residue.  Returns the number of residues renamed."""
    table = load_truncated_resname_table() if table is None else table
    if not table:
        return 0
    groups: OrderedDict = OrderedDict()
    for a in atoms:
        if a.resname in table:
            key = (a.chainID, str(a.resid), a.resname)
            groups.setdefault(key, []).append(a)
    renamed = 0
    for (chain, resid, stem), members in groups.items():
        heavy = {a.name for a in members if _is_heavy(a)}
        resname, reason = choose_resname(stem, heavy, table)
        if resname != stem:
            written = _writable_name(resname)
            if written is None:
                logger.warning(f'residue {stem} {chain}:{resid} looks like {resname} ({reason}), but '
                               f'{resname} has no PDB code pestifer aliases to it, so the name cannot '
                               f'survive the 4-column segment PDB psfgen reads; left as {stem} -- '
                               f'check this residue, or rename it in the input to a code aliased to {resname}')
                continue
            for a in members:
                a.resname = written
            renamed += 1
            shown = resname if written == resname else f'{resname} (written as {written})'
            logger.info(f'residue {stem} {chain}:{resid} identified as {shown} from its atoms ({reason})')
        elif reason.startswith('ambiguous'):
            logger.warning(f'residue {stem} {chain}:{resid} could be a truncated name, but is {reason}; '
                           f'left as {stem}')
    return renamed


def psfgen_segment_resname(resname: str) -> str:
    """
    The residue name to put in a segment PDB that psfgen reads, so that it builds ``resname``.

    Segment PDBs are written in standard columns, which keep four characters of a name: BGLCNA is
    written BGLC and psfgen builds glucose.  A longer name is written as the PDB code a residue alias
    maps to it (BGLCNA -> NAG), which psfgen then aliases back.  A long name with no such code is
    returned unchanged -- it is truncated as before, which cannot be made right from here.
    """
    if len(resname) <= 4:
        return resname
    return _writable_name(resname) or resname
