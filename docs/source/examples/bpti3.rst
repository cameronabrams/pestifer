Example 3: BPTI with One Reduced Disulfide and Some Point Mutations
-------------------------------------------------------------------

Building on Example 2, here we show how to introduce point mutations and how to undo disulfides.  Both of these actions are specified in the ``psfgen`` task under the ``mods`` subdirective:

.. code-block:: yaml

  title: BPTI, no phosphate, some random mutations plus deletion of one disulfide
  tasks:
    - psfgen:
        source:
          id: 6pti
          exclude:
            resnames:
              - PO4
        mods:
          mutations: # showcasing the two shortcode formats
            - A:T11A # threonine to alanine at position 11
            - A:PRO,13,ALA # proline to alanine at position 13
            - A:K15R
            - A:MET,52,LEU
          ssbondsdelete:
            - A_5-A_55
    - md:
        ensemble: minimize
    - solvate:
    - md:
        ensemble: minimize
    - md:
        ensemble: NVT
    - md:
        ensemble: NPT
        nsteps: 200
    - md:
        ensemble: NPT
        nsteps: 400
    - md:
        ensemble: NPT
        nsteps: 800
    - md:
        ensemble: NPT
        nsteps: 1600
    - terminate:
        basename: my_6pti
        package:
          ensemble: NPT
          basename: prod_6pti

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
