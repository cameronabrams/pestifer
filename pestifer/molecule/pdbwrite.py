# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Author PDB coordinate files directly from pestifer's :class:`~pestifer.molecule.atom.Atom`
model, offloading the fixed-column line formatting to pidibble's writer.

pestifer's ``Atom`` is built *from* pidibble records, so a pestifer-own PDB writer would be a
second, independent encoding of the PDB column layout that must agree with pidibble's forever --
exactly the drift that produced the wide-resname column corruption (see
:func:`pestifer.util.coord.standardize_pdb_columns`).  Instead we reconstruct lightweight record
objects from our ``Atom``s and hand them to :class:`pidibble.pdbwrite.PDBWriter`, so there is a
single source of truth for the column model.

We default to pidibble's **CHARMM dialect**, which widens ``resName`` to six columns (so
``BGLCNA``/``ANE5AC`` do not overflow), writes the authoritative segID in columns 73-76, and pins
x/y/z at columns 31-54 regardless of resName width.  The result is a psfgen-ready coordinate PDB
(``coordpdb``/``readpdb``) that needs no column re-anchoring afterward.
"""
import logging

from functools import lru_cache
from types import SimpleNamespace

logger = logging.getLogger(__name__)


@lru_cache(maxsize=None)
def _dialect_formats(dialect: str):
    """
    Return ``(record_formats, custom_formats)`` for a pidibble write dialect.

    A bare :class:`pidibble.pdbparse.PDBParser` loads and dialect-resolves the format tables
    without fetching or parsing anything, so this is a cheap, cached way to obtain the exact
    tables the parser would read/write with.
    """
    from pidibble.pdbparse import PDBParser
    p = PDBParser(dialect=dialect)
    return p.record_formats, p.pdb_format_dict['custom_formats']


def _atom_record(atom) -> SimpleNamespace:
    """
    Build a pidibble-writer record shim from a pestifer :class:`~pestifer.molecule.atom.Atom`.

    :meth:`pidibble.pdbwrite.PDBWriter.emit` reads a record's attributes by the field names in
    the (CHARMM) ATOM/HETATM format, and the composite ``residue`` field by the sub-field names
    in ``ResidueCharmm`` (``resName``/``chainID``/``seqNum``/``iCode``).  A plain attribute bag
    satisfies both.
    """
    resid = atom.resid
    residue = SimpleNamespace(
        resName=atom.resname,
        chainID=atom.chainID,
        seqNum=resid.resseqnum,
        iCode=resid.insertion or '',
    )
    return SimpleNamespace(
        serial=atom.serial,
        name=atom.name,
        altLoc=atom.altloc,
        residue=residue,
        x=atom.x,
        y=atom.y,
        z=atom.z,
        occupancy=atom.occ,
        tempFactor=atom.beta,
        element=atom.elem,
        charge=atom.charge,
        # segID is the authoritative CHARMM segment; fall back to the chain when unset
        segID=atom.segname if atom.segname else atom.chainID,
    )


def write_atoms_pdb(atoms, filename: str = None, dialect: str = 'charmm', end: bool = True) -> list[str]:
    """
    Format an iterable of :class:`~pestifer.molecule.atom.Atom` objects as PDB coordinate lines.

    Each atom is emitted through :class:`pidibble.pdbwrite.PDBWriter` under its own record
    keyword (``ATOM`` or ``HETATM``, from ``atom._PDB_keyword``), so the class distinction is
    preserved.  The atoms are written in the order given -- reserialize/renumber the model
    beforehand if a particular ordering or serial run is required.

    Parameters
    ----------
    atoms : iterable of Atom
        The atoms to write (e.g. an :class:`~pestifer.molecule.atom.AtomList`).
    filename : str, optional
        Destination path.  If omitted, nothing is written and the lines are only returned.
    dialect : str, optional
        pidibble write dialect, ``'charmm'`` (default) or ``'standard'``.
    end : bool, optional
        Append a terminal ``END`` record (default True).

    Returns
    -------
    list of str
        The formatted record lines (without trailing newlines).
    """
    from pidibble.pdbwrite import PDBWriter
    record_formats, custom_formats = _dialect_formats(dialect)
    w = PDBWriter(record_formats, custom_formats)
    lines = [w.emit(_atom_record(a), a._PDB_keyword) for a in atoms]
    if end:
        lines.append('END')
    if filename:
        with open(filename, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')
    return lines
