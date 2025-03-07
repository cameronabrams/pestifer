.. _example 3:

Example 3: BPTI with One Reduced Disulfide and Some Point Mutations
-------------------------------------------------------------------

Building on Example 2, here we show how to introduce point mutations and how to undo disulfides.  Both of these actions are specified in the ``psfgen`` task under the ``mods`` subdirective:

.. literalinclude:: ../../../pestifer/resources/examples/03-bpti-mutant-deldisulf.yaml
    :language: yaml

First, note the ``mutations`` list.  Each element specifies one particular point mutation using a *shortcode*.  There are two allowable shortcodes for a point mutation:

1. ``CHAIN``:``OLRCRESIDOLRC``
2. ``CHAIN``:``TLRC,RESID,TLRC``

``CHAIN`` is the chain ID, ``OLRC`` is a one-letter residue code, and ``RESID`` is the residue sequence number; ``TLRC`` is a three-letter residue code.  Note that both formats are showcased here.

Second, note the ``ssbondsdelete`` list.  Again, a shortcode is used to identify a disulfide to reduce; I think you can see that we are reducing the disulfide between residues 5 and 55.

.. list-table::

    * - .. figure:: pics/yes_disu.png

           BPTI with the 5-55 disulfide intact, showing 
           sidechains for residues T11, P13, K15, and M52.

      - .. figure:: pics/no_disu.png

           BPTI with the 5-55 disulfide reduced, and 
           point mutations T11A, P13A, K15R, and M52L.
