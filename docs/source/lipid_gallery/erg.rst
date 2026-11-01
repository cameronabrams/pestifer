.. _lipid_gallery erg:

ERG
===

**Ergosterol**

==============  ======================================
CHARMM resname  ``ERG``
Net charge      0.00
Atoms           73
Topology file   ``toppar_all36_lipid_cholesterol.str``
Conformers      1 fluid
==============  ======================================

Structure
---------

.. image:: img/erg-structure.svg
   :width: 340px
   :alt: Two-dimensional structure of ERG

Conformers
----------

A sterol is dominated by its rigid fused rings, so the collection carries one minimized structure for it rather than a sampled ensemble.

.. image:: img/erg-conformers.png
   :width: 188px
   :alt: 1 ERG conformer(s) at a shared scale

=========  =============  ================  ================
Conformer  Head-tail (Å)  Max internal (Å)  Min chain extent
=========  =============  ================  ================
``ERG-0``  15.65          17.96             -0.02
=========  =============  ================  ================

*Min chain extent* is how far the least-extended acyl chain descends below its own attachment carbon, as a fraction of that chain's all-trans contour; an all-trans chain scores 1.0, and a conformer under 0.25 is rejected at generation as folded.

Provenance
----------

How the ensemble was generated, read from the entry's own ``info.yaml``:

=======  ==========
Setting  ``ERG``
=======  ==========
Sampler  ``single``
=======  ==========
