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

Drawing the 2D structures
-------------------------

``lipid2d.py`` draws them, using each tool for what it is good at. Open Babel perceives bond
orders from the 3D conformer — in POPC it finds the oleoyl C=C, both ester C=O and the P=O — and
RDKit does the depiction, which can then be rotated deliberately. Open Babel alone lays out from
the 3D coordinates and offers no control over orientation, so a POPC came out lying on its side.

The orientation is derived from the molecule, never set per species, which is what makes it hold
across the collection:

- **Stand it up** so the headgroup sits above the chain attachments. The headgroup reference is
  chosen by element — phosphorus if the lipid has any (every glycerophospholipid and
  sphingomyelin; a cardiolipin has two), otherwise the whole polar N/O region, which is what a
  ceramide's head actually is.
- **Swing each tail** rigidly about its own carbonyl until it hangs downward, fanned slightly so
  two chains leaving nearby attachments do not draw on top of each other. Rotating about the
  carbonyl leaves every bond length and angle inside the chain untouched.

Verified across six lipid classes — PC, PE, PS, PI, sphingomyelin, ceramide and a four-chain
cardiolipin — all of which come out in the same orientation.

Two things it has to get right that are easy to miss. CHARMM names acyl carbons **two** ways:
``C21``/``C318``/``CA16`` puts the chain tag first, while sphingolipids write ``C10S``/``C16F``
with the position first — read only the first and every sphingomyelin and ceramide reports "no
acyl chains". And pestifer's conformer PDBs leave the element column blank, so a reader guessing
from the atom name reads a cardiolipin's ``CA16`` as **calcium**; the element column is filled in
before bond perception runs rather than corrected afterwards.

It needs RDKit, which pestifer carries as the optional ``ligand-paramgen`` extra, plus Open Babel.

**Remaining rough edge:** RDKit's depiction decides where the backbone and headgroup sit, and for
the larger lipids that can be sprawling — a cardiolipin's four tails hang correctly but one curls
around the headgroup. Fixing that needs a template layout for the backbone, not more rotation.
