# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Parsing and writing CHARMM parameter files (.prm and the parameter sections of .str files).

The main class is :class:`CharmmParamFile`, which can be used to parse, merge, filter,
and write CHARMM parameter files.  The typical workflow for generating a minimal
consolidated parameter file for a simulation system is:

1. Parse all the parameter files that were used in the build.
2. Merge them into a single :class:`CharmmParamFile`.
3. Extract the subset of parameters needed for the atom types present in the PSF.
4. Write the result as a single ``.prm`` file.

Usage example::

    from pestifer.charmmff.charmmffprm import CharmmParamFile

    combined = CharmmParamFile()
    for fname in ['par_all36m_prot.prm', 'toppar_water_ions.str']:
        combined.merge(CharmmParamFile.from_file(fname))

    atom_types = {'C', 'NH1', 'CT1', ...}   # from PSF
    minimal = combined.extract_for_atomtypes(atom_types)
    minimal.write('my_system_minimal.prm')
"""

import logging
import re

from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Parameter record dataclasses
# ---------------------------------------------------------------------------

@dataclass
class CharmmBondParam:
    """A CHARMM bond parameter record."""
    type1: str
    type2: str
    Kb: float
    b0: float
    comment: str = ''


@dataclass
class CharmmAngleParam:
    """A CHARMM angle parameter record (with optional Urey-Bradley terms)."""
    type1: str
    type2: str
    type3: str
    Ktheta: float
    Theta0: float
    Kub: float | None = None
    S0: float | None = None
    comment: str = ''


@dataclass
class CharmmDihedralParam:
    """A CHARMM dihedral parameter record.

    Note that multiple records with the same four atom types are allowed
    (different multiplicities).  The wildcard atom type ``X`` may appear
    in position 1 or 4 (or both).
    """
    type1: str
    type2: str
    type3: str
    type4: str
    Kchi: float
    n: int
    delta: float
    comment: str = ''


@dataclass
class CharmmImproperParam:
    """A CHARMM improper dihedral parameter record."""
    type1: str
    type2: str
    type3: str
    type4: str
    Kpsi: float
    n: int
    psi0: float
    comment: str = ''


@dataclass
class CharmmNonbondedParam:
    """A CHARMM Lennard-Jones (nonbonded) parameter record for one atom type."""
    atomtype: str
    ignored: float
    epsilon: float
    Rmin_half: float
    ignored14: float | None = None
    epsilon14: float | None = None
    Rmin_half14: float | None = None
    comment: str = ''


@dataclass
class CharmmNBFixParam:
    """A CHARMM NBFIX record (pairwise LJ override)."""
    type1: str
    type2: str
    Emin: float
    Rmin: float
    comment: str = ''


@dataclass
class CharmmCMAPParam:
    """A CHARMM CMAP correction record."""
    types: list
    grid_size: int
    data: list
    comment: str = ''


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

_PARAM_SECTION_KEYWORDS = frozenset([
    'BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPER',
    'NONBONDED', 'NBFIX', 'CMAP',
])
_END_SECTION_KEYWORDS = frozenset(['ATOMS', 'HBOND', 'END'])


class CharmmParamFile:
    """
    Parses, merges, filters, and writes CHARMM parameter files.

    Supports both standalone ``.prm`` files and the parameter sections
    embedded in ``.str`` (stream) files.

    Parameters
    ----------
    None – use the class-methods :meth:`from_file` or :meth:`from_text`
    to construct an instance from existing data.
    """

    def __init__(self):
        self.bonds: list[CharmmBondParam] = []
        self.angles: list[CharmmAngleParam] = []
        self.dihedrals: list[CharmmDihedralParam] = []
        self.impropers: list[CharmmImproperParam] = []
        self.nonbonded: dict[str, CharmmNonbondedParam] = {}
        self.nbfix: list[CharmmNBFixParam] = []
        self.cmaps: list[CharmmCMAPParam] = []
        # Preserved header line for NONBONDED section (controls cutoffs, etc.)
        self.nonbonded_header: str = (
            'NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -\n'
            'cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5'
        )

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_file(cls, filename: str) -> 'CharmmParamFile':
        """Parse a CHARMM parameter file (or stream file) from disk."""
        with open(filename, 'r') as f:
            text = f.read()
        return cls.from_text(text)

    @classmethod
    def from_text(cls, text: str) -> 'CharmmParamFile':
        """Parse CHARMM parameter content from a string.

        Works for both ``.prm`` files (bare section headers at top level) and
        ``.str`` files (sections inside ``read para … END`` blocks).
        """
        obj = cls()
        sections = obj._extract_param_sections(text)
        for section_text in sections:
            obj._parse_param_section(section_text)
        return obj

    # ------------------------------------------------------------------
    # Internal parsing helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_param_sections(text: str) -> list[str]:
        """Return a list of text blocks that contain parameter data.

        For ``.prm``-style files the whole text is returned as a single block.
        For ``.str``-style files each ``read para … END`` block is extracted.
        """
        lines = text.splitlines()
        # Detect .str-style by looking for 'read para'
        has_read_para = any(re.match(r'^\s*read\s+para', l, re.IGNORECASE) for l in lines)
        if not has_read_para:
            return [text]

        blocks = []
        in_param = False
        current: list[str] = []
        for line in lines:
            stripped = line.strip()
            if re.match(r'^read\s+para', stripped, re.IGNORECASE):
                in_param = True
                current = []
            elif in_param and stripped.upper() == 'END':
                blocks.append('\n'.join(current))
                in_param = False
                current = []
            elif in_param:
                current.append(line)
        # If a block was opened but never closed, keep it
        if in_param and current:
            blocks.append('\n'.join(current))
        return blocks if blocks else [text]

    def _parse_param_section(self, text: str):
        """Parse one parameter block and append records to this instance."""
        section = None
        nb_header_lines: list[str] = []
        nb_header_complete = False

        # CMAP reading state
        cmap_entry: CharmmCMAPParam | None = None
        cmap_remaining: int = 0

        for raw_line in text.splitlines():
            # Split comment from data
            if '!' in raw_line:
                excl_pos = raw_line.index('!')
                comment = raw_line[excl_pos + 1:].strip()
                data_part = raw_line[:excl_pos].strip()
            else:
                comment = ''
                data_part = raw_line.strip()

            if not data_part:
                continue

            upper = data_part.upper()

            # ----------------------------------------------------------
            # Section header detection
            # ----------------------------------------------------------
            if upper.startswith('BOND'):
                section = 'BONDS'
                continue
            elif upper.startswith('ANGLE'):
                section = 'ANGLES'
                continue
            elif upper.startswith('DIHEDRAL'):
                section = 'DIHEDRALS'
                continue
            elif upper.startswith('IMPROPER'):
                section = 'IMPROPER'
                continue
            elif upper.startswith('NONBONDED'):
                section = 'NONBONDED'
                nb_header_lines = [raw_line]
                nb_header_complete = False
                # If the line does NOT end with '-' the header is one line
                if not data_part.rstrip().endswith('-'):
                    nb_header_complete = True
                    self.nonbonded_header = raw_line
                continue
            elif upper.startswith('NBFIX'):
                section = 'NBFIX'
                continue
            elif upper.startswith('CMAP'):
                section = 'CMAP'
                continue
            elif upper.startswith('HBOND') or any(upper.startswith(k) for k in _END_SECTION_KEYWORDS):
                section = None
                continue

            if section is None:
                continue

            tokens = data_part.split()

            # ----------------------------------------------------------
            # NONBONDED header continuation (lines after NONBONDED until
            # we get to actual atom entries)
            # ----------------------------------------------------------
            if section == 'NONBONDED' and not nb_header_complete:
                # A continuation of the keyword header: either cutnb... or
                # any line that is not an atom type + 3 floats
                try:
                    float(tokens[1])
                    float(tokens[2])
                    float(tokens[3])
                    # This looks like a proper atom entry
                    nb_header_complete = True
                    self.nonbonded_header = '\n'.join(nb_header_lines)
                except (IndexError, ValueError):
                    # Still header continuation
                    nb_header_lines.append(raw_line)
                    if not data_part.rstrip().endswith('-'):
                        nb_header_complete = True
                        self.nonbonded_header = '\n'.join(nb_header_lines)
                    continue

            # ----------------------------------------------------------
            # CMAP data reading
            # ----------------------------------------------------------
            if section == 'CMAP':
                if cmap_remaining > 0:
                    # Accumulate grid data
                    try:
                        vals = [float(t) for t in tokens]
                    except ValueError:
                        continue  # skip any spurious non-numeric line
                    cmap_entry.data.extend(vals)
                    cmap_remaining -= len(vals)
                    if cmap_remaining <= 0:
                        self.cmaps.append(cmap_entry)
                        cmap_entry = None
                        cmap_remaining = 0
                elif len(tokens) >= 9:
                    # New CMAP entry header: 8 types + grid_size
                    try:
                        grid_size = int(tokens[8])
                    except (ValueError, IndexError):
                        continue
                    cmap_entry = CharmmCMAPParam(
                        types=tokens[:8],
                        grid_size=grid_size,
                        data=[],
                        comment=comment,
                    )
                    cmap_remaining = grid_size * grid_size
                continue

            # ----------------------------------------------------------
            # Standard section parsing
            # ----------------------------------------------------------
            try:
                if section == 'BONDS' and len(tokens) >= 4:
                    self.bonds.append(CharmmBondParam(
                        type1=tokens[0], type2=tokens[1],
                        Kb=float(tokens[2]), b0=float(tokens[3]),
                        comment=comment,
                    ))

                elif section == 'ANGLES' and len(tokens) >= 5:
                    Kub = float(tokens[5]) if len(tokens) >= 7 else None
                    S0  = float(tokens[6]) if len(tokens) >= 7 else None
                    self.angles.append(CharmmAngleParam(
                        type1=tokens[0], type2=tokens[1], type3=tokens[2],
                        Ktheta=float(tokens[3]), Theta0=float(tokens[4]),
                        Kub=Kub, S0=S0,
                        comment=comment,
                    ))

                elif section == 'DIHEDRALS' and len(tokens) >= 7:
                    self.dihedrals.append(CharmmDihedralParam(
                        type1=tokens[0], type2=tokens[1],
                        type3=tokens[2], type4=tokens[3],
                        Kchi=float(tokens[4]), n=int(tokens[5]),
                        delta=float(tokens[6]),
                        comment=comment,
                    ))

                elif section == 'IMPROPER' and len(tokens) >= 7:
                    self.impropers.append(CharmmImproperParam(
                        type1=tokens[0], type2=tokens[1],
                        type3=tokens[2], type4=tokens[3],
                        Kpsi=float(tokens[4]), n=int(tokens[5]),
                        psi0=float(tokens[6]),
                        comment=comment,
                    ))

                elif section == 'NONBONDED' and len(tokens) >= 4:
                    ignored14 = epsilon14 = Rmin_half14 = None
                    if len(tokens) >= 7:
                        ignored14    = float(tokens[4])
                        epsilon14    = float(tokens[5])
                        Rmin_half14  = float(tokens[6])
                    self.nonbonded[tokens[0]] = CharmmNonbondedParam(
                        atomtype=tokens[0],
                        ignored=float(tokens[1]),
                        epsilon=float(tokens[2]),
                        Rmin_half=float(tokens[3]),
                        ignored14=ignored14,
                        epsilon14=epsilon14,
                        Rmin_half14=Rmin_half14,
                        comment=comment,
                    )

                elif section == 'NBFIX' and len(tokens) >= 4:
                    self.nbfix.append(CharmmNBFixParam(
                        type1=tokens[0], type2=tokens[1],
                        Emin=float(tokens[2]), Rmin=float(tokens[3]),
                        comment=comment,
                    ))

            except (ValueError, IndexError) as exc:
                logger.debug(f'Skipping malformed parameter line: {repr(raw_line)} ({exc})')

    # ------------------------------------------------------------------
    # Merging
    # ------------------------------------------------------------------

    def merge(self, other: 'CharmmParamFile') -> None:
        """Merge *other* into this instance (in-place).

        For NONBONDED entries, later values override earlier ones (same as
        CHARMM's ``READ PARAM APPEND`` semantics).  For all list-based
        sections (BONDS, ANGLES, …) entries are simply appended; duplicates
        are not removed, which mirrors the behaviour of CHARMM when reading
        multiple parameter files.
        """
        self.bonds.extend(other.bonds)
        self.angles.extend(other.angles)
        self.dihedrals.extend(other.dihedrals)
        self.impropers.extend(other.impropers)
        self.nonbonded.update(other.nonbonded)
        self.nbfix.extend(other.nbfix)
        self.cmaps.extend(other.cmaps)

    # ------------------------------------------------------------------
    # Filtering
    # ------------------------------------------------------------------

    def extract_for_atomtypes(self, atomtypes: set[str]) -> 'CharmmParamFile':
        """Return a new :class:`CharmmParamFile` containing only the records
        relevant to *atomtypes*.

        A record is included when every non-wildcard (``X``) atom type in
        that record appears in *atomtypes*.  This is a conservative filter:
        it may include a handful of wildcard entries that are not strictly
        needed, but it never omits a record that is needed.

        Parameters
        ----------
        atomtypes : set of str
            The set of atom types present in the target PSF.

        Returns
        -------
        CharmmParamFile
            A new instance containing the filtered records.
        """
        def match(types: list[str]) -> bool:
            return all(t == 'X' or t in atomtypes for t in types)

        result = CharmmParamFile()
        result.nonbonded_header = self.nonbonded_header
        result.bonds     = [b for b in self.bonds     if match([b.type1, b.type2])]
        result.angles    = [a for a in self.angles    if match([a.type1, a.type2, a.type3])]
        result.dihedrals = [d for d in self.dihedrals if match([d.type1, d.type2, d.type3, d.type4])]
        result.impropers = [i for i in self.impropers if match([i.type1, i.type2, i.type3, i.type4])]
        result.nonbonded = {t: v for t, v in self.nonbonded.items() if t in atomtypes}
        result.nbfix     = [n for n in self.nbfix     if match([n.type1, n.type2])]
        result.cmaps     = [c for c in self.cmaps     if match(c.types)]
        return result

    # ------------------------------------------------------------------
    # Writing
    # ------------------------------------------------------------------

    def write(self, filename: str, title: str = '') -> None:
        """Write this parameter set to *filename* in CHARMM format.

        Parameters
        ----------
        filename : str
            Output file path.
        title : str, optional
            A short description written into the file header.
        """
        lines: list[str] = []

        # Header
        if title:
            lines.append(f'* {title}')
        else:
            lines.append('* Minimal CHARMM parameter file generated by pestifer')
        lines.append('*')
        lines.append('')

        def _comment(c: str) -> str:
            return f' ! {c}' if c else ''

        # BONDS
        if self.bonds:
            lines.append('BONDS')
            for b in self.bonds:
                lines.append(
                    f'{b.type1:<6s} {b.type2:<6s}  {b.Kb:10.3f}  {b.b0:8.4f}{_comment(b.comment)}'
                )
            lines.append('')

        # ANGLES
        if self.angles:
            lines.append('ANGLES')
            for a in self.angles:
                ub = (f'  {a.Kub:8.3f}  {a.S0:8.4f}' if a.Kub is not None else '')
                lines.append(
                    f'{a.type1:<6s} {a.type2:<6s} {a.type3:<6s}  {a.Ktheta:10.3f}  {a.Theta0:8.4f}'
                    f'{ub}{_comment(a.comment)}'
                )
            lines.append('')

        # DIHEDRALS
        if self.dihedrals:
            lines.append('DIHEDRALS')
            for d in self.dihedrals:
                lines.append(
                    f'{d.type1:<6s} {d.type2:<6s} {d.type3:<6s} {d.type4:<6s}  '
                    f'{d.Kchi:10.4f}  {d.n:1d}  {d.delta:8.2f}{_comment(d.comment)}'
                )
            lines.append('')

        # IMPROPER
        if self.impropers:
            lines.append('IMPROPER')
            for i in self.impropers:
                lines.append(
                    f'{i.type1:<6s} {i.type2:<6s} {i.type3:<6s} {i.type4:<6s}  '
                    f'{i.Kpsi:10.4f}  {i.n:1d}  {i.psi0:8.2f}{_comment(i.comment)}'
                )
            lines.append('')

        # CMAP
        if self.cmaps:
            lines.append('CMAP')
            for c in self.cmaps:
                lines.append(' '.join(c.types) + f'  {c.grid_size}')
                # Write grid data in lines of 5 values (CHARMM convention)
                for i in range(0, len(c.data), 5):
                    chunk = c.data[i:i + 5]
                    lines.append(''.join(f'{v:9.2f}' for v in chunk))
                lines.append('')

        # NONBONDED
        if self.nonbonded:
            lines.append(self.nonbonded_header)
            for t, nb in self.nonbonded.items():
                nb14 = ''
                if nb.ignored14 is not None:
                    nb14 = (
                        f'  {nb.ignored14:8.4f}  {nb.epsilon14:12.6f}  {nb.Rmin_half14:10.6f}'
                    )
                lines.append(
                    f'{nb.atomtype:<8s}  {nb.ignored:8.4f}  {nb.epsilon:12.6f}  '
                    f'{nb.Rmin_half:10.6f}{nb14}{_comment(nb.comment)}'
                )
            lines.append('')

        # NBFIX
        if self.nbfix:
            lines.append('NBFIX')
            for n in self.nbfix:
                lines.append(
                    f'{n.type1:<8s} {n.type2:<8s}  {n.Emin:12.6f}  {n.Rmin:10.6f}{_comment(n.comment)}'
                )
            lines.append('')

        lines.append('END')
        lines.append('')

        with open(filename, 'w') as f:
            f.write('\n'.join(lines))

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------

    def summary(self) -> str:
        """Return a one-line summary of the parameter counts."""
        return (
            f'bonds={len(self.bonds)}, angles={len(self.angles)}, '
            f'dihedrals={len(self.dihedrals)}, impropers={len(self.impropers)}, '
            f'nonbonded={len(self.nonbonded)}, nbfix={len(self.nbfix)}, '
            f'cmap={len(self.cmaps)}'
        )
