:orphan:

.. note::
   The ``domainswap`` task is **retired** and no longer available in a build config.
   This page is kept for reference; the code remains in the source tree but is not
   registered as a usable task.

.. _subs_buildtasks_domainswap:

domainswap
----------

The ``domainswap`` task performs a targeted (steered) molecular-dynamics simulation that swaps a defined structural domain between chains — for example, exchanging the position of a domain on chain A with the equivalent domain on chain B.  A domain, selected by the VMD atomselect string ``swap_domain_def``, is biased toward its partner chain's location while anchor regions (``anchor_domain_def``) are held in place, using collective-variable harmonic restraints in NAMD.  The pairing of chains to swap is given by ``chain_directional_swaps`` (e.g. ``[[A, B], [B, A]]`` swaps A and B).

This is an advanced task; the complete parameter list is in :ref:`config_ref tasks domainswap`, and the underlying Tcl is the :ref:`tcl-domainswap` script.

.. code-block:: yaml

   tasks:
     - ...
     - domainswap:
         swap_domain_def: "protein and resid 1 to 100"
         anchor_domain_def: "protein and resid 200 to 300"
         chain_directional_swaps: [[A, B], [B, A]]
         force_constant: 200.0
         nsteps: 10000
