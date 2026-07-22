# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Post-build cleanup of psfgen-written PSF topology.

Several CHARMM ``modify_res.str`` residues -- every GFP-type chromophore (CRO,
GYS, CR2, CRQ, DYG, NRQ, ...), the beta-amino acids, SUI, MAL, and others --
declare an explicit *backward* peptide bond (``BOND N -C``) in addition to the
standard forward link (``BOND C +N``).  Standard residues declare only the
forward link; the backward link is implied by the preceding residue's forward
link.  When psfgen builds one of these residues in the middle of a chain, the
preceding residue's ``+N`` and the modified residue's ``-C`` describe the *same*
peptide bond, so the bond is created twice.  ``regenerate angles dihedrals``
then turns that duplicate bond into a degenerate angle (e.g. ``C-NG2S1-C`` with
the same carbon on both ends), which NAMD rejects with
``UNABLE TO FIND ANGLE PARAMETERS``.

:func:`dedupe_psf_topology` rewrites the finished PSF, removing duplicate bonds
and dropping any angle/dihedral/improper that references an atom twice (which
can only arise from such a malformed bond list).  A clean PSF is left byte-for-
byte untouched.
"""

import logging
import os

logger = logging.getLogger(__name__)

# number of atom indices per record in each counted section we clean
_SECTION_K = {'NBOND': 2, 'NTHETA': 3, 'NPHI': 4, 'NIMPHI': 4}
# records per line psfgen writes for each section
_GROUPS_PER_LINE = {'NBOND': 4, 'NTHETA': 3, 'NPHI': 2, 'NIMPHI': 2}


def _canon(name, t):
    """Order-independent key so a record and its reverse count as one."""
    if name == 'NBOND':
        return tuple(sorted(t))
    if name == 'NTHETA':
        i, j, k = t
        return (min(i, k), j, max(i, k))
    # dihedrals / impropers: a-b-c-d is the same term as d-c-b-a
    return min(t, t[::-1])


def _filter_section(name, records):
    """Drop degenerate (repeated-index) records and de-duplicate the rest."""
    k = _SECTION_K[name]
    seen = set()
    kept = []
    n_degenerate = n_duplicate = 0
    for t in records:
        if len(set(t)) < k:
            n_degenerate += 1
            continue
        key = _canon(name, t)
        if key in seen:
            n_duplicate += 1
            continue
        seen.add(key)
        kept.append(t)
    return kept, n_degenerate, n_duplicate


def dedupe_psf_topology(psf_path):
    """
    Remove duplicate bonds and degenerate angle/dihedral/improper terms from a
    psfgen-written PSF in place.  The file is rewritten only if something was
    actually removed, so clean structures are left unchanged.

    Parameters
    ----------
    psf_path : str
        Path to the PSF file to clean.

    Returns
    -------
    dict
        Per-section ``{name: (n_degenerate, n_duplicate)}`` for the terms
        removed; empty if the PSF was already clean.
    """
    if not psf_path or not os.path.exists(psf_path):
        return {}
    with open(psf_path) as fh:
        lines = fh.read().splitlines()
    if not lines or not lines[0].startswith('PSF'):
        return {}
    width = 10 if 'EXT' in lines[0] else 8

    out = []
    removed = {}
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        name = next((s for s in _SECTION_K if f'!{s}' in line), None)
        if name is None:
            out.append(line)
            i += 1
            continue
        count = int(line.split()[0])
        k = _SECTION_K[name]
        need = count * k
        toks = []
        j = i + 1
        while len(toks) < need and j < n:
            toks.extend(int(x) for x in lines[j].split())
            j += 1
        records = [tuple(toks[p:p + k]) for p in range(0, need, k)]
        kept, n_deg, n_dup = _filter_section(name, records)
        if n_deg or n_dup:
            removed[name] = (n_deg, n_dup)
        label = line.split('!', 1)[1]
        out.append(f'{len(kept):>{width}} !{label}')
        gpl = _GROUPS_PER_LINE[name]
        buf = []
        for idx, t in enumerate(kept):
            buf.extend(t)
            if (idx + 1) % gpl == 0:
                out.append(''.join(f'{v:>{width}}' for v in buf))
                buf = []
        if buf:
            out.append(''.join(f'{v:>{width}}' for v in buf))
        i = j

    if removed:
        with open(psf_path, 'w') as fh:
            fh.write('\n'.join(out) + '\n')
        detail = ', '.join(f'{name} (-{d} degenerate, -{u} duplicate)'
                           for name, (d, u) in removed.items())
        logger.info(f'dedupe_psf_topology: cleaned {os.path.basename(psf_path)}: {detail}')
    else:
        logger.debug(f'dedupe_psf_topology: {os.path.basename(psf_path)} already clean')
    return removed
