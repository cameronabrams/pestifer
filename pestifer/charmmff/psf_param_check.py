# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Force-field-consistency preflight for an incoming (foreign) PSF.

Verifies that a build's CHARMM release resolves **every atom type and every
bond/angle/dihedral/improper type-tuple** present in a PSF, so a topology built
elsewhere (e.g. CHARMM-GUI, against a different CHARMM version) fails **up front**
-- naming the exact residue and term -- instead of dying cryptically deep in a
later NAMD run (``DIDN'T FIND vdW PARAMETER FOR ATOM TYPE ...`` or a missing bonded
parameter fatal).

This module is deliberately NAMD-free and depends only on a parsed
:class:`~pestifer.psfutil.psfcontents.PSFContents` and a merged
:class:`~pestifer.charmmff.charmmffprm.CharmmParamFile`, so it is fully unit-testable.

Matching mirrors CHARMM's own parameter lookup:

- **atom types** (vdW/nonbonded): exact membership, no wildcards.
- **bonds**: unordered pair (``b0``/``Kb`` are symmetric).
- **angles**: the triple, canonicalised by reversibility.
- **dihedrals**: the quartet, canonicalised by reversibility; the wildcard atom
  type ``X`` may appear in positions 1 and/or 4 (the common ``X A B X`` central-pair
  form is indexed for speed; rarer partial-wildcard entries are matched directly).
- **impropers**: matched against each parameter's quartet with ``X`` wildcards in
  any position, in either direction (CHARMM improper matching is permissive).

**Atom types are compared case-insensitively**, because CHARMM's are.  A single force-field
file routinely spells one type three ways -- ``toppar_all36_prot_na_combined.str`` declares
``MASS -1 ON2B``, writes ``ATOM OH ON2b`` in the ``TP1`` patch that uses it, gives every bonded
parameter as ``ON2b``, and gives the vdW record as ``ON2B``.  psfgen writes the **MASS**
spelling, so a real phosphotyrosine PSF carries ``ON2B`` (verified against a built PTR system);
its vdW lookup then succeeds while every one of its bonded terms fails a case-sensitive match.
Every comparison here therefore runs over upper-cased tuples, while the PSF's own spelling is
what gets reported.

- **CMAP cross-terms**: the 8-type tuple, matched exactly.  CHARMM's CMAP records carry no
  wildcards (the shipped release defines six, all fully specified) and the tuple is a phi
  quartet followed by a psi quartet, so it is directional -- reversing it would swap the two
  angles and mean something else.  Checked against real builds: every cross-term in each of
  them matched a parameter record exactly, with no reversal needed.

Terms are deduplicated by their atom-type tuple, so the scan is cheap even on a
multi-million-atom system: each distinct type-tuple is looked up once.
"""
import logging

from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass
class MissingParameters:
    """The unresolved atom types and bonded terms found in a PSF.

    Each list holds ``(term, residue_label)`` tuples, where ``term`` is a string
    (an atom type, or a ``-``-joined type-tuple) and ``residue_label`` names a
    representative residue (``'<resname> <segid><resid>'``) carrying that term.
    """
    atomtypes: list = field(default_factory=list)
    bonds: list = field(default_factory=list)
    angles: list = field(default_factory=list)
    dihedrals: list = field(default_factory=list)
    impropers: list = field(default_factory=list)
    cmaps: list = field(default_factory=list)

    def any(self) -> bool:
        return bool(self.atomtypes or self.bonds or self.angles
                    or self.dihedrals or self.impropers or self.cmaps)


def _serial_chunks(lines, n):
    """Yield successive ``n``-tuples of integer atom serials from raw PSF section lines."""
    for line in lines:
        vals = [int(x) for x in line.split()]
        for i in range(0, len(vals) - n + 1, n):
            yield tuple(vals[i:i + n])


def _u(types) -> tuple:
    """Upper-case a tuple of atom types, for case-insensitive comparison."""
    return tuple(t.upper() for t in types)


def _quartet_match(pattern, quartet) -> bool:
    """True if a parameter ``pattern`` (may contain ``'X'``) matches ``quartet`` in either direction."""
    def m(p, q):
        return all(pi == 'X' or pi == qi for pi, qi in zip(p, q))
    return m(pattern, quartet) or m(pattern, quartet[::-1])


def _index_params(param):
    """Build fast lookup structures from a merged :class:`CharmmParamFile`."""
    bondset = {tuple(sorted(_u((b.type1, b.type2)))) for b in param.bonds}
    angleset = set()
    for a in param.angles:
        t = _u((a.type1, a.type2, a.type3))
        angleset.add(min(t, t[::-1]))
    dih_exact = set()
    dih_wild_mid = set()   # central pair (b,c) for the common 'X b c X' form
    dih_wild_other = []    # rarer partial-wildcard quartets, matched directly
    for d in param.dihedrals:
        q = _u((d.type1, d.type2, d.type3, d.type4))
        if 'X' not in q:
            dih_exact.add(min(q, q[::-1]))
        elif q[0] == 'X' and q[3] == 'X' and q[1] != 'X' and q[2] != 'X':
            dih_wild_mid.add(min((q[1], q[2]), (q[2], q[1])))
        else:
            dih_wild_other.append(q)
    improper_quartets = [_u((i.type1, i.type2, i.type3, i.type4)) for i in param.impropers]
    cmapset = {_u(tuple(c.types)) for c in param.cmaps}
    return (bondset, angleset, dih_exact, dih_wild_mid, dih_wild_other,
            improper_quartets, cmapset)


def check_psf_parameters(psf, param) -> MissingParameters:
    """Return the atom types and bonded terms in *psf* that *param* does not resolve.

    Parameters
    ----------
    psf : PSFContents
        A parsed PSF (atoms + raw token sections; no ``parse_topology`` required).
    param : CharmmParamFile
        The merged parameter set for the build's CHARMM release.
    """
    (bondset, angleset, dih_exact, dih_wild_mid,
     dih_wild_other, improper_quartets, cmapset) = _index_params(param)

    serial_to_atom = {int(a.serial): a for a in psf.atoms}

    def typ(s):
        return serial_to_atom[s].atomtype

    def label(s):
        a = serial_to_atom[s]
        return f'{a.resname} {a.segname}{a.resid.resid}'

    missing = MissingParameters()

    # --- atom types (vdW / nonbonded): exact membership, no wildcards ---
    seen_at = set()
    nonbonded = {k.upper() for k in param.nonbonded}
    for a in psf.atoms:
        at = a.atomtype
        if at.upper() in seen_at:
            continue
        seen_at.add(at.upper())
        if at.upper() not in nonbonded:
            missing.atomtypes.append((at, f'{a.resname} {a.segname}{a.resid.resid}'))

    # --- bonded terms: dedupe by type-tuple so each distinct term is checked once ---
    def scan(section, n, key_fn, ok_fn, out):
        seen = set()
        for serials in _serial_chunks(psf.token_lines.get(section, []), n):
            try:
                types = tuple(typ(s) for s in serials)
            except KeyError:
                continue  # a term referencing an atom serial not in ATOM (malformed PSF) -- skip
            utypes = _u(types)
            k = key_fn(utypes)
            if k in seen:
                continue
            seen.add(k)
            if not ok_fn(utypes):
                out.append(('-'.join(types), label(serials[0])))

    scan('BOND', 2, lambda t: tuple(sorted(t)),
         lambda t: tuple(sorted(t)) in bondset, missing.bonds)
    scan('THETA', 3, lambda t: min((t[0], t[1], t[2]), (t[2], t[1], t[0])),
         lambda t: min((t[0], t[1], t[2]), (t[2], t[1], t[0])) in angleset, missing.angles)

    def dihedral_ok(t):
        a, b, c, d = t
        if min(t, t[::-1]) in dih_exact:
            return True
        if min((b, c), (c, b)) in dih_wild_mid:
            return True
        return any(_quartet_match(w, t) for w in dih_wild_other)

    scan('PHI', 4, lambda t: min(t, t[::-1]), dihedral_ok, missing.dihedrals)

    def improper_ok(t):
        return any(_quartet_match(w, t) for w in improper_quartets)

    scan('IMPHI', 4, lambda t: min(t, t[::-1]), improper_ok, missing.impropers)

    # CMAP: exact 8-tuple, no wildcards, direction significant (phi quartet then psi quartet)
    scan('CRTERM', 8, lambda t: t, lambda t: t in cmapset, missing.cmaps)

    return missing


def subtract(a: MissingParameters, b: MissingParameters) -> MissingParameters:
    """Return the entries of *a* whose term does not also appear in *b*.

    Used to split "this term is in none of the loaded parameter files" from "this term was
    loaded but the consolidated file does not carry it" -- two failures that look identical
    from NAMD but have opposite remedies.
    """
    out = MissingParameters()
    for kind in ('atomtypes', 'bonds', 'angles', 'dihedrals', 'impropers', 'cmaps'):
        seen = {term for term, _ in getattr(b, kind)}
        setattr(out, kind, [(t, w) for t, w in getattr(a, kind) if t not in seen])
    return out


def format_missing(missing: MissingParameters, release: str,
                   header: str = '', advice: str = '') -> str:
    """Render a :class:`MissingParameters` into an actionable multi-line error message.

    *header* and *advice* replace the first and last lines, for callers other than the
    incoming-PSF preflight (the build path has a different cause and a different remedy).
    """
    lines = [header or
             f'The incoming PSF references force-field terms that charmmff {release} does not '
             f'resolve. This usually means the PSF was built against a different CHARMM version.']
    if missing.atomtypes:
        lines.append(f'  Unresolved atom type(s) [{len(missing.atomtypes)}] (no vdW parameter):')
        for at, where in missing.atomtypes:
            lines.append(f"    '{at}'  (e.g. residue {where})")
    for label, items in (('bond term', missing.bonds), ('angle term', missing.angles),
                         ('dihedral term', missing.dihedrals),
                         ('improper term', missing.impropers),
                         ('CMAP cross-term', missing.cmaps)):
        if items:
            lines.append(f'  Unresolved {label}(s) [{len(items)}]:')
            for term, where in items:
                lines.append(f"    '{term}'  (e.g. residue {where})")
    lines.append(advice or
                 'Rebuild the topology against this release, add the missing parameters to your '
                 'build, or bring the correct parameter/stream files into the run directory.')
    return '\n'.join(lines)
