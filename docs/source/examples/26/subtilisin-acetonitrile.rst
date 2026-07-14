.. _example subtilisin-acetonitrile:

Example 26: Subtilisin Carlsberg in acetonitrile
------------------------------------------------

`PDB ID 1scd <https://www.rcsb.org/structure/1scd>`_ is the structure of the serine protease subtilisin Carlsberg.  Like :ref:`Example 23 <example subtilisin-dmso>` (DMSO) and :ref:`Example 24 <example subtilisin-acetone>` (acetone), this example builds the enzyme in a **non-aqueous solvent** -- here acetonitrile, another classic medium for the study of enzyme catalysis in organic solvents.

As with acetone, pestifer ships no pre-equilibrated acetonitrile box, so the ``solvate`` task **builds one on demand**: because ``ACN`` is defined in the CGenFF force field, pestifer packs a periodic box, minimizes and NPT-equilibrates it, and **caches** the result under ``~/.pestifer/pdbrepository/<release>/solvent/ACN/`` for reuse.  The first build pays a one-time cost to generate the box; from then on it is as fast as a shipped solvent.

What makes acetonitrile special is a subtlety in the force field.  Acetonitrile's heavy-atom skeleton (H\ :sub:`3`\ C--C≡N) is **linear**, so the H--C--C≡N dihedral is *geometrically degenerate* -- there is no defined torsion angle and no torsional energy.  CGenFF therefore ships **no parameter** for it, yet psfgen's ``autogenerate`` still emits the term from the bond graph, and NAMD refuses to run a term with no parameter (``UNABLE TO FIND DIHEDRAL PARAMETERS FOR HGA3 CG331 CG1N1 NG1T1``).  Pestifer supplies the physically-correct value -- a reviewed ``k = 0`` entry in ``custom/toppar_pestifer_dihedral_fills.prm`` -- so the box equilibrates cleanly.  This is *not* a blanket zero-fill: any other unrecognized missing parameter is still surfaced as an actionable error naming the exact residue and atom-type term, never silently masked.

The rest of the build mirrors the DMSO and acetone examples: the equilibrated acetonitrile box is tiled to fill the cell, the system's net charge is neutralized by replacing a few acetonitrile molecules with counter-ions (VMD's ``autoionize`` only works for water, so pestifer does this by solvent replacement), and a staged, progressively longer NPT schedule relaxes the box to its equilibrium density.

.. literalinclude:: ../../../../pestifer/resources/examples/26/inputs/subtilisin-acetonitrile.yaml
    :language: yaml

.. task-table:: ../../../../pestifer/resources/examples/26/inputs/subtilisin-acetonitrile.yaml
