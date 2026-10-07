.. _lipid_gallery index:

Lipid gallery (prototype)
=========================

.. note::

   **This is a one-lipid experiment, not a feature.** It exists to show what a gallery page can
   look like and what building 277 of them would cost, so that the full project can be judged
   before it is started. Only :ref:`POPC <lipid_gallery popc>` is here.

.. toctree::
   :maxdepth: 1

   popc

What a page carries
-------------------

Four things, in the order a reader wants them:

1. **Identity** — CHARMM resname, the synonym from the topology, charge, atom count, and which
   topology file defines it. All read from the force field, none transcribed.
2. **Structure** — a 2D diagram, generated from the conformer coordinates with Open Babel and
   embedded as SVG so it scales.
3. **Conformers** — every shipped conformer rendered at **one shared scale**, with the
   per-conformer geometry beside them. The shared scale is the point: framing each image to its
   own molecule would hide exactly the differences in extension the ensemble exists to provide.
4. **Provenance** — the generation settings the entry records about itself, so a reader can
   reproduce or audit it.

What the prototype already found
--------------------------------

POPC ships **four distinct conformers of ten** — one draw repeated seven times, byte-identical.
The duplication was known in the abstract for POPE and recorded as "4–5 distinct of 10"; nobody
had measured POPC, which is among the most-used entries in the collection. A contact sheet made
it obvious on sight.

That is the argument for the gallery, and it is worth more than the documentation value: a
rendering is a different instrument from a table of numbers, and it caught something the numbers
had not.

Cost and open questions, if this is built out
---------------------------------------------

**Rendering.** About 4 s per conformer with the headless VMD path, so roughly 2700 conformers is
under three hours, trivially parallel. The Tcl renderer must be fed to VMD on **stdin**; with
``-e`` every image comes back containing nothing but the corner axes. The existing
``scripts/pestifer-snapshot`` wrapper exists for that reason and the gallery renderer follows it.

**Repository size is the real constraint.** A contact sheet plus an SVG is ~240 KB per entry,
which is **~66 MB for 277 entries** — too much to commit. Three options, none yet chosen:

- generate the images at docs-build time, which makes the build need VMD and Open Babel;
- commit smaller sheets (fewer conformers shown, or lower resolution), losing the detail that
  made the duplication visible;
- host the images outside the repository and reference them, which breaks offline builds.

**Per-entry scale.** The prototype uses one hand-calibrated scale for POPC. A real build needs a
scale computed per entry from its most-extended conformer — shared within an entry, not across
entries, since a 24-carbon lipid and a 12-carbon one cannot usefully share a frame.

**2D structures are imperfect.** Open Babel lays the diagram out from 3D coordinates, so bond
crossings are common and it warns about blank atom names in pestifer's PDBs. Acceptable for a
prototype; a real gallery may want templates or a different layout engine.
