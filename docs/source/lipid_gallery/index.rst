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
2. **Structure** — a 2D diagram in the conventional lipid orientation: headgroup at the top,
   both tails hanging down, no explicit hydrogens. Embedded as SVG so it scales.
3. **Conformers** — every shipped conformer rendered at **one shared scale**, with the
   per-conformer geometry beside them. The shared scale is the point: framing each image to its
   own molecule would hide exactly the differences in extension the ensemble exists to provide.
4. **Provenance** — the generation settings the entry records about itself, so a reader can
   reproduce or audit it.

Cost and open questions, if this is built out
---------------------------------------------

**Rendering.** About 4 s per conformer with the headless VMD path, so roughly 2700 conformers is
under three hours, trivially parallel. The Tcl renderer is fed to VMD on stdin, as
``scripts/pestifer-snapshot`` does.

**Images are committed, as 64-color palette PNGs.** A contact sheet saved that way is about
37 KB with no visible loss (RMSE 0.006 against the full-color render), so a full build is about
16 MB: a sheet and a ~27 KB structure drawing for each of the 219 base lipids, and a sheet for
each of the 58 ``__Lo`` ensembles. Docs are excluded from the wheel and the sdist, so this costs
clone size only. Generating the images while Read the Docs builds is not possible (VMD cannot be
installed there, and rendering takes hours); hosting them outside the repository becomes worth
its plumbing only if the conformer collection is regenerated often, since each full re-render
adds another ~16 MB to history.

**How a sheet is made.** Each conformer is rendered by ``lipidsnap.tcl`` straight at tile size,
``-size 188x332 -scale 0.0954``, tiled five across with ``magick montage -tile 5x2 -geometry
+0+0 -background white``, and reduced with ``magick ... -depth 8 -colors 64``. POPC's two
ensembles share that one scale.

**Per-entry scale.** That scale was calibrated by hand for POPC. A full build needs a scale
computed per entry from its most-extended conformer — shared within an entry, not across
entries, since a 24-carbon lipid and a 12-carbon one cannot usefully share a frame.

Drawing the 2D structures
-------------------------

``lipid2d.py`` draws them: Open Babel perceives bond orders from a 3D conformer, and RDKit does
the depiction. The orientation is derived from the molecule, never set per species, so it holds
across the collection:

- **Stand it up** so the headgroup sits above the chain attachments. The headgroup reference is
  phosphorus if the lipid has any (every glycerophospholipid and sphingomyelin; a cardiolipin has
  two), otherwise the whole polar N/O region, which is what a ceramide's head is.
- **Swing each tail** rigidly about its last substituted carbon until it hangs downward, fanned
  slightly so two chains leaving nearby attachments do not draw on top of each other. A rigid
  rotation leaves every bond length and angle inside the chain untouched.

No stereochemistry is drawn: these are identity and shape references, and stereochemistry
perceived from one conformer would mark every chain CH₂ as a stereocenter.

``batch2d.py`` runs it over the whole collection and draws **276 of 277 entries** in the same
orientation — every PC, PE, PS, PG, PA, PI, sphingomyelin, ceramide, lyso lipid and four-chain
cardiolipin, and the sterols in their own conventional pose, 3-OH up and the iso-octyl tail
down. The one it declines is ``CHM1``, a twelve-atom ring model compound with no acyl chain.

It needs RDKit, which pestifer carries as the optional ``ligand-paramgen`` extra, plus Open Babel.
