.. _subs_buildtasks_psfgen_mods_patches:

patches
-------

A ``patches`` directive is used to specify patches, which are changes to single amino acid residues enabled by the CHARMM36 topologies.  The syntax for specifying a patch is as follows:

``<patchname>:<chain>:<seqnum>``

``<patchname>`` is the CHARMM36 patch name, ``<chain>`` is the chain ID, and ``<seqnum>`` is the residue number.  For example, if you wanted to protonate an aspartate at position 12 in chain A, your ``patches`` mod directive would look like this:

.. code-block:: yaml

    patches:
      - ASPP:A:12

``ASPP`` is the CHARMM36 patch name for protonating an aspartate, ``A`` is the chain ID, and ``12`` is the residue number.  This will apply the specified patch to the residue at that position in the specified chain.

A patch is not a modified residue
+++++++++++++++++++++++++++++++++

``<patchname>`` must be a CHARMM ``PRES``.  Most post-translational modifications are *not*
patches: CHARMM defines them as whole ``RESI`` residues, which are installed with
:ref:`mutations <subs_buildtasks_psfgen_mods_mutations>` instead.  Handing one to ``patches`` is
reported rather than silently ignored, because psfgen would otherwise skip the unrecognized name
and produce an unmodified system that looks like a successful build:

.. code-block:: text

    ERROR> unknown patch name(s) in this task's `patches`: 'SEP' is a CHARMM residue (RESI),
           not a patch (PRES); a whole modified residue is applied as a mutation, not a patch

Patches are still the right tool for protonation states (``ASPP``, ``GLUP``, ``LSN``, ``HS2``),
for linkages, and for the phospho patches when you want explicit control of charge state --
``SP1``/``SP2`` on a serine rather than the monoanionic ``RESI SEP``.  See
:ref:`Post-translational modifications <post_translational_modifications>`.

