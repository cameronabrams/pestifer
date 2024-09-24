Example 4: BPTI with an Extra Disulfide Added
---------------------------------------------

Using the ``mods`` subdirective, one can introduce new disulfides into an existing structure.  This example introduces a disulfide linking residues 11 and 34.

.. code-block::  yaml

  title: BPTI, no phosphate, introducing a disulfide via mutations
  tasks:
    - psfgen:
        source:
          id: 6pti
          exclude:
            resnames:
              - PO4 # we don't need no stinkin phosphates
        mods:
          mutations: # get me two cysteines, stat!
            - A:T11C
            - A:V34C
          ssbonds:
            - A_11-A_34  # now ligation!
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

Note that this required first mutating the residues at positions 11 and 34 to cysteines, and *then* introducing the disulfide mod.
