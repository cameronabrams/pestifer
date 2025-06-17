.. _subs_runtasks_psfgen_mods_patches:

patches
-------

A ``patches`` directive is used to specify patches, which are changes to single amino acid residues enabled by the CHARMM36 topologies.  The syntax for specifying a patch is as follows:

``<patchname>:<chain>:<seqnum>``

``<patchname>`` is the CHARMM36 patch name, ``<chain>`` is the chain ID, and ``<seqnum>`` is the residue number.  For example, if you wanted to protonate an aspartate at position 12 in chain A, your ``patches`` mod directive would look like this:

.. code-block:: yaml

    patches:
      - ASPP:A:12

``ASPP`` is the CHARMM36 patch name for protonating an aspartate, ``A`` is the chain ID, and ``12`` is the residue number.  This will apply the specified patch to the residue at that position in the specified chain.