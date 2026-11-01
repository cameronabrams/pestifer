.. _lipid_gallery popc:

POPC
====

.. note::

   **Prototype page.** One lipid, built by hand, to decide whether a full gallery of all 277
   shipped entries is worth generating. See :ref:`lipid_gallery index` for what a full build
   would cost and what is still unsettled.

**3-palmitoyl-2-oleoyl-D-glycero-1-Phosphatidylcholine** — a 16:0/18:1 phosphatidylcholine, and
one of the most-used entries in the shipped collection.

=====================  ==========================================
CHARMM resname         ``POPC``
Net charge             0.00
Atoms                  134
Topology file          ``top_all36_lipid.rtf``
Stream                 ``lipid``
Conformers             10 fluid, 10 liquid-ordered (``POPC__Lo``)
=====================  ==========================================

Structure
---------

.. image:: img/popc-structure.svg
   :width: 340px
   :alt: Two-dimensional structure of POPC, headgroup up and both acyl tails hanging down

Drawn the way lipids conventionally are: headgroup at the top, both tails hanging down, no
explicit hydrogens. That orientation is derived from the molecule rather than set by hand, so it
holds for every non-sterol lipid — see :ref:`lipid_gallery index`.

Conformers
----------

The grid packer draws a conformer per lipid from an ensemble, so what the ensemble contains is
what a built bilayer is made of. POPC ships two: a fluid one, used by default, and a
liquid-ordered one, used for a leaflet declared ``phase: Lo``. Every conformer on this page, in
both ensembles, is rendered **at one shared scale**, so differences in extension are real rather
than an artifact of per-image framing.

Fluid (``POPC``)
~~~~~~~~~~~~~~~~

.. image:: img/popc-conformers.png
   :width: 100%
   :alt: Ten POPC conformers rendered at a shared scale

Read left to right, top row then bottom. The spread is what the ensemble is for: the packer draws
one shape per lipid, so a bilayer is built out of this variety rather than out of one repeated
conformer.

================  ===============  ==================  ==================
Conformer         Head-tail (Å)    Max internal (Å)    Min chain extent
================  ===============  ==================  ==================
``POPC-00``       16.11            29.42               0.83
``POPC-01``       14.10            30.63               0.72
``POPC-02``       16.99            28.86               0.83
``POPC-03``       14.99            28.23               0.75
``POPC-04``       6.86             25.95               0.34
``POPC-05``       15.60            30.33               0.80
``POPC-06``       14.86            28.53               0.76
``POPC-07``       15.42            28.05               0.74
``POPC-08``       12.04            26.77               0.62
``POPC-09``       16.66            28.94               0.79
================  ===============  ==================  ==================

*Min chain extent* is how far the least-extended acyl chain descends below its own attachment
carbon, as a fraction of that chain's all-trans contour. An all-trans chain scores 1.0; anything
under 0.25 is rejected at generation as folded. Every conformer here passes, ``POPC-04`` least
comfortably — its sn-2 chain descends 7.4 Å of a 21.6 Å contour.

Liquid-ordered (``POPC__Lo``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generated the same way with a bias toward trans chain torsions, which is what makes the tails
straighter and more nearly parallel.

.. image:: img/popc-lo-conformers.png
   :width: 100%
   :alt: Ten liquid-ordered POPC conformers rendered at the same scale as the fluid ones

================  ===============  ==================  ==================
Conformer         Head-tail (Å)    Max internal (Å)    Min chain extent
================  ===============  ==================  ==================
``POPC-00``       16.87            30.18               0.87
``POPC-01``       18.33            30.50               0.88
``POPC-02``       15.34            27.57               0.73
``POPC-03``       9.42             30.67               0.48
``POPC-04``       17.54            30.27               0.90
``POPC-05``       15.32            29.47               0.79
``POPC-06``       17.21            30.63               0.89
``POPC-07``       16.97            30.41               0.87
``POPC-08``       15.97            29.08               0.76
``POPC-09``       17.75            29.53               0.84
================  ===============  ==================  ==================

Against the fluid ensemble, the mean head-tail length is 16.1 Å rather than 14.4 Å and the mean
min chain extent 0.80 rather than 0.72. ``POPC-03`` is the least extended, with one tail bent
sideways near its end.

Provenance
----------

How each ensemble was generated, read from the entry's own ``info.yaml``:

====================  ==============  ==============
Setting               ``POPC``        ``POPC__Lo``
====================  ==============  ==============
Sampler               ``mc``          ``mc``
Phase                 none (fluid)    Lo
Trans bias            0.0             2.5
Chain order           0.212           0.286
MC seed               0               1
Chains                2               2
Cylinder APL          60.0 Å²         60.0 Å²
Cylinder inflation    1.9             1.9
Equilibration steps   20000           20000
Decorrelation steps   3000            3000
Max MC angle          30°             30°
Radius scale          0.8             0.8
Axial bias            1.0             1.0
====================  ==============  ==============

Every value above is read from the shipped entry, not chosen for the page: the collection records
its own generation settings so an entry can be reproduced or audited later.
