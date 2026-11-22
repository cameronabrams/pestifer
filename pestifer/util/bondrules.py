# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Topological bond rules shared by every code path that rotates a molecule about a bond.

Three places in pestifer turn a branch rigidly about a bridge of the bond graph: the
athermal-MC conformer sampler (``charmmff.athermal_mc``), the ring-piercing resolver's hinge
search (``psfutil.psfring.RingChecker.pendant_axes``), and psfgen's pendant declasher
(``tasks.psfgen.PsfgenTask._rotatable_bonds``).  A C=C bond is a bridge, so each of them
offered it as a hinge.  A turn about it leaves a cis double bond at an arbitrary angle,
minimization drops it into the nearest CHARMM well, and ``X CEL1 CEL1 X`` has multiplicity 2,
so trans is a minimum too and no MD restores cis.

The sampler was fixed first (2026-11-20) with a private rule; the resolver kept flipping
double bonds during membrane builds (found 2026-11-22 in example 17: 14 oleoyl bonds trans and
2 sphingosine bonds cis in one quilt, from a collection with none).  The rule lives here so the
next caller that rotates about a bond gets it by importing, not by remembering.
"""

_MASSES = (('H', 1.008), ('C', 12.011), ('N', 14.007), ('O', 15.999),
           ('P', 30.974), ('S', 32.06))


def element_of_mass(mass: float) -> str:
    """Coarse element symbol from atomic weight (enough to tell C/H/O/N/P/S apart)."""
    for sym, m in _MASSES:
        if abs(mass - m) < 1.0:
            return sym
    return '?'


def is_multiple_bond(elem_a: str, degree_a: int, elem_b: str, degree_b: int) -> bool:
    """True when a bond joins two under-coordinated carbons, i.e. is a C=C (or C#C) bond.

    ``degree_*`` is the atom's FULL coordination, hydrogens included -- count it on the whole
    molecule, not on a fragment cut out of it, or an atom at the cut reads as unsaturated.  A
    carbon with fewer than four bonded neighbors is sp2 or sp, so a bond between two of them is
    a double or triple bond, about which there is no torsion to turn.

    The test is topological on purpose: it needs no bond orders or atom types, which a PSF does
    not carry.  It also excludes the single bond in the middle of a conjugated diene (both ends
    are sp2), which loses one legitimate degree of freedom; losing a torsion is the safe
    direction of error.  It is meaningful only when hydrogens are explicit; a caller that may
    see a hydrogen-free graph must decide that case itself.
    """
    return elem_a == 'C' and elem_b == 'C' and degree_a < 4 and degree_b < 4
