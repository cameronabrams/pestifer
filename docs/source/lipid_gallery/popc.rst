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
Conformers shipped     10
=====================  ==========================================

Structure
---------

.. image:: img/popc-structure.svg
   :width: 320px
   :alt: Two-dimensional structure of POPC

Conformers
----------

The grid packer draws a conformer per lipid from this ensemble, so what the ensemble contains is
what a built bilayer is made of. All ten are rendered **at one shared scale**, so differences in
extension are real rather than an artifact of per-image framing.

.. image:: img/popc-conformers.png
   :width: 100%
   :alt: Ten POPC conformers rendered at a shared scale

Read left to right, top row then bottom. The first three are visibly curled; the rest are
extended and — as the table below shows — mostly identical.

================  ===============  ==================  ==================
Conformer         Head-tail (Å)    Max internal (Å)    Min chain extent
================  ===============  ==================  ==================
``POPC-00``       14.99            28.23               0.75
``POPC-01``       6.86             25.95               0.34
``POPC-02``       14.86            28.53               0.76
``POPC-03``       16.66            28.94               0.79
``POPC-04``       16.66            28.94               0.79
``POPC-05``       16.66            28.94               0.79
``POPC-06``       16.66            28.94               0.79
``POPC-07``       16.66            28.94               0.79
``POPC-08``       16.66            28.94               0.79
``POPC-09``       16.66            28.94               0.79
================  ===============  ==================  ==================

*Min chain extent* is how far the least-extended acyl chain descends below its own attachment
carbon, as a fraction of that chain's all-trans contour. An all-trans chain scores 1.0; anything
under 0.25 is rejected at generation as folded. Every conformer here passes.

.. warning::

   **POPC ships four distinct conformers, not ten.** ``POPC-03`` through ``POPC-09`` are
   byte-identical — one draw repeated seven times, confirmed by hashing the coordinates, not
   inferred from the table. So 70% of POPC placements in any built bilayer receive the same
   shape.

   This is pre-existing and not caused by any recent change, but its size was not known until
   this page was drawn: the duplication was on record as "4–5 distinct of 10" for POPE, and
   nobody had looked at POPC. It is the strongest argument for building the rest of the gallery
   — a contact sheet makes it obvious at a glance, where a column of numbers did not.

Provenance
----------

How this ensemble was generated, read from the entry's own ``info.yaml``:

====================  =============
Sampler               ``mc``
Phase                 none (fluid)
Trans bias            0.0
Chain order           0.249
MC seed               0
Chains                2
Cylinder APL          60.0 Å²
Cylinder inflation    1.9
Equilibration steps   20000
Decorrelation steps   3000
Max MC angle          30°
Radius scale          0.8
Axial bias            1.0
====================  =============

Every value above is read from the shipped entry, not chosen for the page: the collection records
its own generation settings so an entry can be reproduced or audited later.
