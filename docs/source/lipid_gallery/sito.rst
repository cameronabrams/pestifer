.. _lipid_gallery sito:

SITO
====

**beta-Sitosterol**

==============  ======================================
CHARMM resname  ``SITO``
Net charge      0.00
Atoms           80
Topology file   ``toppar_all36_lipid_cholesterol.str``
Conformers      1 fluid
==============  ======================================

Structure
---------

.. image:: img/sito-structure.svg
   :width: 340px
   :alt: Two-dimensional structure of SITO

Conformers
----------

A sterol is dominated by its rigid fused rings, so the collection carries one minimized structure for it rather than a sampled ensemble.

.. image:: img/sito-conformers.png
   :width: 188px
   :alt: 1 SITO conformer(s) at a shared scale

==========  =============  ================  ================
Conformer   Head-tail (Å)  Max internal (Å)  Min chain extent
==========  =============  ================  ================
``SITO-0``  15.23          18.30             -0.01
==========  =============  ================  ================

*Min chain extent* is how far the least-extended acyl chain descends below its own attachment carbon, as a fraction of that chain's all-trans contour; an all-trans chain scores 1.0, and a conformer under 0.25 is rejected at generation as folded.

Provenance
----------

How the ensemble was generated, read from the entry's own ``info.yaml``:

=======  ==========
Setting  ``SITO``
=======  ==========
Sampler  ``single``
=======  ==========
