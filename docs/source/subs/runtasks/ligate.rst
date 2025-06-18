.. _subs_runtasks_ligate:

ligate 
------

A ``ligate`` task is used to build missing protein loops that are specified in the input structure file.  It is very often the case that a structure file declares "missing" or "zero-occupancy" residues that are part of the physical molecule that was analyzed but for whatever reason are not resolved.  Pestifer by default will build these in as unstructured loops.  The algorithm that does this is as follows:

1. For each missing residue, include the appropriate ``residue`` command in the ``segment`` block in the automatically-generated psfen input script.
2. When ``guesscoords`` is run in psfgen, these missing residues are inserted with default internal coordinates.  More than one or two of them, and you have a situation in which the carbonyl carbon of the last inserted residue can be very far away to the amide nitrogen atom to which it is bonded.
3. The appropriate NTER patch is applied to this residue, and a CTER patch is applied to the last inserted residue.
4. A steered MD simulation is performed to pull the CTER to the NTER for each gap.
5. A ``psfgen`` run is peformed in which the CTER and NTER patches are reversed and a LINK patch is applied to make the bond.

For insertions of two or fewer residues, no steering is performed.

This is illustrated in :ref:`example 6`.