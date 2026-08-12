.. _importing_systems:

Importing a pre-built system
============================

Not every system starts from a PDB entry.  You may already have a PSF and coordinates built
somewhere else — CHARMM-GUI, another toolchain, or an earlier pestifer run — and want to run it,
solvate it further, embed it in a membrane, or edit it.  Pestifer takes such a system in as
**pipeline state**, and everything downstream treats it exactly like a system it built itself.

This page is the end-to-end walkthrough.  The task pages carry the full option lists:
:ref:`continuation <subs_buildtasks_continuation>` for the import itself and
:ref:`psfgen <subs_buildtasks_psfgen>` for the edits.

State enters through ``continuation``
-------------------------------------

There is deliberately **no** ``psfgen`` source key for a foreign PSF.  A pre-built system is a
*state*, not a coordinate source, so it enters through :ref:`continuation <subs_buildtasks_continuation>`
(or through :ref:`merge <subs_buildtasks_merge>`, which combines several pre-built systems into one):

.. code-block:: yaml

    tasks:
      - continuation:
          psf: from_charmm_gui.psf
          pdb: from_charmm_gui.pdb
          xsc: from_charmm_gui.xsc
      - md:
          ensemble: NPT
          nsteps: 10000
      - terminate:
          basename: my_system

A ``psf`` plus either a ``pdb`` or a binary ``coor`` is the minimum; ``xsc`` (the periodic cell) and
``vel`` are optional and carried forward when given.  From here the system flows into MD, membrane
embedding, solvation, and packaging like any other.

The force-field preflight
-------------------------

A foreign PSF may have been built against a different CHARMM release than the one pestifer is using,
and the failure mode for that is ugly: a missing parameter surfaces as a cryptic
``DIDN'T FIND vdW PARAMETER`` abort from NAMD, potentially hours into a build.

So ``continuation`` runs a consistency check first, controlled by ``verify_parameters`` (default
``true``).  It inventories the structure — a per-segid composition summary, with a warning for any
residue name it cannot map to a known segment type — and then verifies that **every** atom type and
every bond, angle, dihedral, and improper type-tuple in the PSF resolves against the build's CHARMM
release.  Anything unresolved stops the build immediately, naming the residue and the offending term.

Editing an imported system
--------------------------

A ``psfgen`` task placed after the import edits it.  You do not tell ``psfgen`` which mode to use —
it infers that from what precedes it: after a ``fetch`` it builds a topology from fetched
coordinates, and after a state with no fetched source it edits the topology it was handed.

Which edits are possible depends on whether they disturb the segment topology, and the distinction
matters because it decides how much of the incoming system survives.

**Additive edits keep the incoming topology verbatim.**  Pestifer ``readpsf``\ s the structure and
never re-derives its segments, so patches, custom residues, and an already-solvated, equilibrated box
all come through untouched.  Patches, disulfides, links, backbone and side-chain rotations, and
glycan grafts that *extend* a terminal residue are all in this family.

**Re-segmenting edits rebuild.**  A mutation, insertion, or deletion changes a segment's sequence,
and a segment loaded by ``readpsf`` is immutable — so these reconstruct the molecule from the
incoming PSF plus coordinates and re-segment it.  Sequence and coordinates are preserved, standard
links and disulfides are re-derived, and the periodic cell is carried forward.  The cost is that
topology is re-derived from residue templates, so any non-standard connectivity that is not
expressible as a recognized patch or link is **not** carried through.  Pestifer guards this up front:
if any residue in the incoming PSF has no ``RESI`` in the build's release, it stops and names it
rather than rebuilding a system it cannot reproduce.

.. code-block:: yaml

    tasks:
      - continuation:
          psf: from_charmm_gui.psf
          pdb: from_charmm_gui.pdb
          xsc: from_charmm_gui.xsc
      - psfgen:
          mods:
            ssbonds:
              - A_42-A_58
      - md:
          ensemble: minimize
      - terminate:
          basename: my_edited_system

What will not work
------------------

Mods that read *source metadata* have nothing to read on a foreign PSF, because it carries none:
``biological_assembly``, ``SEQADV``-driven mutations, ``REMARK 465`` loop healing, and
``terminal_tails`` all hard-error rather than silently doing nothing.  Those belong to a build that
starts from a fetched structure.

.. raw:: html

    <div class="autogen-footer">
        <p>Author: Cameron F. Abrams &nbsp;&nbsp;&nbsp; Contact: <a href="mailto:cfa22@drexel.edu">cfa22@drexel.edu</a></p>
    </div>
