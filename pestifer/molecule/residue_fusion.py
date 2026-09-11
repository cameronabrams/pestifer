# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Fuse a ligand the PDB deposits separately into the polymer residue CHARMM defines it as part of.

The PDB writes a lipidated amino acid as two residues joined by a ``LINK`` -- a ``MYR`` ligand
bonded to a ``GLY``'s backbone nitrogen -- while CHARMM defines a single ``GLYM`` whose heavy
atoms are exactly the union of the two, name for name.  Nothing reconciles the two
representations, and the failure is not merely a stopped build: ``MYR`` is itself a CHARMM
residue (free myristic acid, in the detergent stream), so a build can succeed and produce a
**detached fatty acid** sitting beside the protein it was covalently attached to.

The fusion runs on ATOMS, before they are grouped into residues.  Relabelling the ligand's atoms
with the partner's ``resid`` and the combined ``resname`` lets the existing grouping perform the
merge, so no Residue objects have to be spliced together.

The ``LINK`` is the trigger, never the ligand alone.  A free fatty acid that happens to be in a
structure is left exactly where it is.
"""
import logging

from ..core.labels import Labels

logger = logging.getLogger(__name__)


def fuse_linked_ligands(atoms, links) -> int:
    """Fuse every ligand-plus-LINK pair that CHARMM defines as one residue.

    Parameters
    ----------
    atoms : AtomList
        All atoms read from the structure.  Modified in place: the ligand's atoms take the
        partner residue's ``resid`` and the combined ``resname``, and the partner's atoms take
        the combined ``resname``.
    links : LinkList
        Links read from the structure.  A link that has been consumed by a fusion is removed --
        it describes a bond that is now internal to one residue, and leaving it would ask psfgen
        to patch a residue to itself.

    Returns
    -------
    int
        How many ligands were fused.
    """
    table = Labels.fusible_ligands
    if not table or not links:
        return 0
    consumed, fused = [], 0
    for link in list(links):
        for (lig, partner), (combined, lig_atom, partner_atom) in table.items():
            # a link is stored in either order, so try both
            for a_res, a_name, a_cid, a_rid, b_res, b_name, b_cid, b_rid in (
                (link.resname1, link.name1, link.chainID1, link.resid1,
                 link.resname2, link.name2, link.chainID2, link.resid2),
                (link.resname2, link.name2, link.chainID2, link.resid2,
                 link.resname1, link.name1, link.chainID1, link.resid1)):
                if (a_res, b_res) != (lig, partner):
                    continue
                if (a_name, b_name) != (lig_atom, partner_atom):
                    # same residue pair, different bond -- not this modification
                    continue
                lig_atoms = [x for x in atoms if x.resname == lig and x.chainID == a_cid
                             and x.resid == a_rid]
                par_atoms = [x for x in atoms if x.resname == partner and x.chainID == b_cid
                             and x.resid == b_rid]
                if not lig_atoms or not par_atoms:
                    logger.debug(f'{lig}/{partner} link found but atoms are missing; not fusing')
                    continue
                # the ligand's hydrogens are named on a different convention (H21/H22 against
                # CHARMM's H2A/H2B); psfgen rebuilds them, so drop them rather than alias 27 names
                dropped = [x for x in lig_atoms if x.name.startswith('H')]
                for x in dropped:
                    atoms.remove(x)
                for x in (x for x in lig_atoms if x not in dropped):
                    x.resname, x.resid, x.chainID = combined, b_rid, b_cid
                for x in par_atoms:
                    x.resname = combined
                consumed.append(link)
                fused += 1
                logger.info(f'fused {lig} {a_cid}:{a_rid} into {partner} {b_cid}:{b_rid} as '
                            f'{combined} ({len(lig_atoms) - len(dropped)} heavy atoms carried, '
                            f'{len(dropped)} hydrogens dropped for psfgen to rebuild)')
                break
            else:
                continue
            break
    for link in consumed:
        links.remove(link)
    return fused
