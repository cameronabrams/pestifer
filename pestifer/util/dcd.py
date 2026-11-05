# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Minimal reader for NAMD DCD trajectories: header, per-frame steps, unit cell and coordinates.

NAMD writes a CHARMM-format DCD with 32-bit Fortran record markers.  The header's ``ISTART`` is
the step of the FIRST frame (``firsttimestep + dcdfreq``, not ``firsttimestep``) and ``NSAVC``
is ``dcdfreq``, so frame ``k`` was written at step ``ISTART + k*NSAVC``.  That is what lets a DCD
frame be matched to a ``PRESSUREPROFILE`` record or an ``.xst`` row by step rather than by
position.  ``NSET`` is not trusted: a run that was killed leaves it stale, so the frame count is
taken from the file size.

With ``DCDUnitCell yes`` every frame starts with a 6-double record ``[A, gamma, B, beta, alpha,
C]``; only the three lengths are used here.  The cell ORIGIN is not in a DCD -- read it from the
``.xst``.
"""
from __future__ import annotations

import os
import struct
from dataclasses import dataclass

import numpy as np


@dataclass
class DCDHeader:
    path: str
    natoms: int
    nframes: int
    istart: int
    nsavc: int
    has_cell: bool
    endian: str
    first_frame_offset: int
    frame_bytes: int

    def steps(self) -> np.ndarray:
        """The MD step each frame was written at."""
        return self.istart + self.nsavc * np.arange(self.nframes, dtype=np.int64)


def _record(f, endian):
    raw = f.read(4)
    if len(raw) < 4:
        raise EOFError
    n = struct.unpack(endian + 'i', raw)[0]
    body = f.read(n)
    if struct.unpack(endian + 'i', f.read(4))[0] != n:
        raise ValueError('mismatched Fortran record markers')
    return body


def read_header(path: str) -> DCDHeader:
    with open(path, 'rb') as f:
        raw = f.read(4)
        endian = '<' if struct.unpack('<i', raw)[0] == 84 else '>'
        if struct.unpack(endian + 'i', raw)[0] != 84:
            raise ValueError(f'{path}: not a DCD file (first record is not 84 bytes)')
        f.seek(0)
        hdr = _record(f, endian)
        if hdr[:4] != b'CORD':
            raise ValueError(f'{path}: not a coordinate DCD (no CORD tag)')
        icntrl = struct.unpack(endian + '20i', hdr[4:84])
        istart, nsavc = icntrl[1], icntrl[2]
        has_cell = bool(icntrl[10])
        if icntrl[8]:
            raise ValueError(f'{path}: fixed atoms in a DCD are not supported')
        _record(f, endian)                              # title block
        natoms = struct.unpack(endian + 'i', _record(f, endian))[0]
        first = f.tell()
    frame_bytes = 3 * (8 + 4 * natoms) + (8 + 48 if has_cell else 0)
    nframes = (os.path.getsize(path) - first) // frame_bytes
    return DCDHeader(path, natoms, int(nframes), istart, nsavc, has_cell, endian, first,
                     frame_bytes)


def iter_frames(path: str, atoms: np.ndarray | None = None, frames=None):
    """Yield ``(k, cell_lengths, xyz)`` for each frame ``k``.

    ``cell_lengths`` is ``(a, b, c)`` or ``None`` without unit-cell records; ``xyz`` is
    ``(n, 3)`` float64 for ``atoms`` (all atoms if ``None``).  ``frames`` restricts to those
    frame indices.
    """
    h = read_header(path)
    sel = slice(None) if atoms is None else np.asarray(atoms)
    want = range(h.nframes) if frames is None else frames
    dt = np.dtype(h.endian + 'f4')
    with open(path, 'rb') as f:
        for k in want:
            f.seek(h.first_frame_offset + k * h.frame_bytes)
            cell = None
            if h.has_cell:
                c6 = np.frombuffer(_record(f, h.endian), dtype=h.endian + 'f8')
                cell = (float(c6[0]), float(c6[2]), float(c6[5]))
            xyz = np.empty((h.natoms if atoms is None else len(atoms), 3))
            for d in range(3):
                xyz[:, d] = np.frombuffer(_record(f, h.endian), dtype=dt)[sel]
            yield k, cell, xyz
