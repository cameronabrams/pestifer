.. _build_provenance:

What a build records: environment and citations
===============================================

A pestifer build is a chain of other people's programs — VMD, psfgen, NAMD — driven by a Python
process against a particular CHARMM force-field release, starting from coordinates somebody
deposited.  Two questions follow from that, and a build answers both without being asked:

* **What actually ran?**  Reproducing a result means knowing every version in the chain.
* **Who gets the credit?**  A single YAML input can pull in half a dozen papers, and no user
  should be expected to work out which.

Both answers go into the log at ``INFO``, so they reach the console and the diagnostics log
alike, and every line carries a uniform prefix so either block can be recovered from a log with
one ``grep``.

The environment record
----------------------

Emitted **before the first task**, prefixed ``environment:``:

.. code-block:: console

    $ pestifer run bpti1 2>&1 | grep '^INFO> environment:'
    INFO> environment: ---- build environment ----
    INFO> environment: python      CPython 3.13.13  (/home/you/pestifer/.venv/bin/python3)
    INFO> environment: platform    Linux-7.1.8-1-default-x86_64-with-glibc2.43  on yourhost
    INFO> environment: charmmff    feb26
    INFO> environment: namd mode   cpu
    INFO> environment: catdcd      [/usr/local/bin/catdcd]
    INFO> environment: charmrun    [/usr/local/bin/charmrun]
    INFO> environment: namd3       NAMD 3.0.2 for Linux-x86_64-multicore  [/usr/local/bin/namd3]
    INFO> environment: vmd         2.0.0 (March 25, 2026)  [/usr/local/bin/vmd]
    INFO> environment: pestifer    3.17.1
    INFO> environment: packages    matplotlib 3.11.0, networkx 3.6.1, numpy 2.5.1, ...
    INFO> environment: ---------------------------

What it captures, and why each part matters:

``python``, ``platform``
    Interpreter implementation, version and path, plus the machine.

``charmmff``
    The force-field release directory the build resolved to.

``namd mode``, executables
    The **resolved path and reported version of each program the build may invoke**.  This is the
    part a log could not previously reconstruct at all: a machine may carry several VMD and NAMD
    builds, and nothing recorded which one ran.  NAMD is probed by handing it an empty config, so
    it announces itself before objecting to it; a CUDA build additionally reports the CUDA version
    it was compiled against, which a GPU-resident result depends on.  ``namd3`` is always probed,
    because even a GPU build invokes it for tasks carrying ``cpu-override``, while ``namd3gpu`` is
    probed only in GPU mode, since probing the CUDA build binds a GPU device a CPU run never
    touches.

``pestifer``, ``packages``
    Pestifer's version and the resolved version of every runtime dependency it declares.  The
    dependency names come from installed metadata rather than a hand-kept list, so the record
    cannot drift out of step with ``pyproject.toml``.  Pestifer's own version comes from the
    source tree, not from ``importlib.metadata``, which records the version at install time and so
    goes stale on an editable install the moment a release is cut.

No probe can fail a build.  Anything that does not answer is recorded as ``unknown``.

The citation report
-------------------

Emitted **after the run**, so it sits beside the results it applies to, prefixed ``citations:``.
It has two halves: the software the build used, and the coordinates it was given.

Software citations
~~~~~~~~~~~~~~~~~~

NAMD, VMD, psfgen and the CHARMM36 protein set are owed by every build.  The rest depend on what
the input asked for, and **each conditional entry states why it applies**, so the report is a
derivation you can check rather than a list to be taken on faith:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Added when
     - Citations
   * - always
     - NAMD; VMD; psfgen; CHARMM36m (protein); CHARMM force-field overview
   * - ``namd.namd-version`` is 3 (the default)
     - NAMD 3
   * - a ``make_membrane_system`` or ``membrane_equilibrate`` task
     - CHARMM36 lipids
   * - the ``psfgen`` task has non-empty ``mods.grafts``
     - CHARMM36 carbohydrates
   * - ``solvate`` names a non-water solvent, or ``charmmff.user_custom`` supplies files
     - CGenFF (plus both CGenFF-program papers in the ``user_custom`` case)
   * - a ``pdb2pqr`` task
     - PDB2PQR; APBS; PROPKA

Coordinate citations
~~~~~~~~~~~~~~~~~~~~

Coordinates enter a build in three places, and all three are reported with the role they played:
the ``fetch`` task's structure, graft donors (``mods.grafts``), and C-terminal fusion donors
(``mods.Cfusions``).  Grafts and fusions are read with pestifer's own shortcode parsers, so this
cannot drift from the syntax a build accepts.

Identifiers are read from the structure files the build **already downloaded** — no network call
is made on the citation report's account, and the citation therefore belongs to the coordinates
actually used.  Because the report is emitted after the run, and the
:ref:`terminate <subs_buildtasks_terminate>` task tidies the run directory by sweeping every
intermediate into ``<basename>-artifacts.tar.gz``, the originals are usually gone by then; the
report reads them out of that archive instead.  Either way it is the file the build used.

They are reported as **DOI and PMID rather than as reference strings**, deliberately.  The PDB
format stores author names and titles in upper case and the original capitalization is not
recoverable from them: ``A.B.MCDERMOTT`` cannot be turned back into ``McDermott, A.B.``, and
guessing sentence case from an ALL-CAPS title mangles things like ``HIV-1``.  An identifier
cannot be subtly wrong in that way.

Example: the SARS-CoV-2 BA.2 spike build (:ref:`Example 15 <example sars cov2 spike ba2>`), which
starts from one entry and grafts glycans from three others:

.. code-block:: text

    citations: --- and the coordinates this build was given ---
    citations: [PDB 7XIX] doi:10.1038/s41586-022-04980-y pmid:35714668
    citations:     ^ the structure this build starts from
    citations: [PDB 4B7I] doi:10.1021/ja306068g pmid:23025485
    citations:     ^ coordinates were grafted from it
    citations: [PDB 2WAH] doi:10.1016/j.jmb.2009.02.033 pmid:19236877
    citations:     ^ coordinates were grafted from it
    citations: [PDB 4BYH] doi:10.1073/pnas.1310657110 pmid:23929778
    citations:     ^ coordinates were grafted from it

Resolving identifiers requires ``pidibble`` 1.10.0 or newer, which is the declared floor.  Where
they cannot be resolved, the structures are **named anyway**, with a note saying why — a report
that claims to be complete has to say when it is not.  Coordinates with no deposited citation of
their own, such as an AlphaFold model or a donor file you built yourself, are likewise named
rather than dropped:

.. code-block:: text

    citations: note: P22033 (alphafold) has no deposited citation of its own -- the structure
                     this build starts from

What the report does not know
-----------------------------

The software half is derived from your **input config**, so it covers what you asked for.  It
cannot see what your starting structure already contained: a PDB entry that arrives carrying
glycans, nucleic acids or a bound ligand may owe further CHARMM parameter-set citations that no
reading of the config could reveal.  The report says so on its last lines.  Check the force-field
files your build pulled in — the ``environment:`` block records the release — if you need to be
certain.

A worked example
----------------

Every conditional path in one place, from the AlphaFold methylmalonyl-CoA mutase build
(:ref:`Example 20 <example methylmalonyl-coa-mutase>`), which fetches a model rather than a PDB
entry and protonates it with ``pdb2pqr``:

.. code-block:: text

    citations: ---- please cite ----
    citations: Pestifer orchestrates other people's software.  Work that uses a system
    citations: prepared by this build should cite the following:
    citations: [NAMD] Phillips JC, et al. Scalable molecular dynamics with NAMD. ... doi:10.1002/jcc.20289
    citations: [VMD] Humphrey W, Dalke A, Schulten K. VMD: Visual molecular dynamics. ... doi:10.1016/0263-7855(96)00018-5
    citations: [psfgen] Ribeiro J, Radak B, Stone J, Gullingsrud J, Saam J, Phillips J. psfgen User's Guide, v2.0 (2020).
    citations: [CHARMM36m (protein)] Huang J, et al. ... doi:10.1038/nmeth.4067
    citations: [CHARMM force field] Vanommeslaeghe K, MacKerell AD. ... doi:10.1016/j.bbagen.2014.08.004
    citations: [NAMD 3] Phillips JC, et al. Scalable molecular dynamics on CPU and GPU architectures with NAMD. ... doi:10.1063/5.0014475
    citations:     ^ included because this build ran NAMD 3
    citations: [PDB2PQR] Dolinsky TJ, et al. ... Nucleic Acids Res 35, W522-W525 (2007).
    citations:     ^ included because this build has a pdb2pqr task
    citations: [APBS/PDB2PQR] Jurrus E, et al. ... Protein Sci 27, 112-128 (2018).
    citations:     ^ included because this build has a pdb2pqr task
    citations: [PROPKA] Olsson MHM, Sondergaard CR, Rostkowski M, Jensen JH. ... J Chem Theory Comput 7, 525-537 (2011).
    citations:     ^ included because this build has a pdb2pqr task, which assigns pKa with PROPKA
    citations: --- and the coordinates this build was given ---
    citations: note: P22033 (alphafold) has no deposited citation of its own -- the structure this build starts from
    citations: Software citations are derived from your input config, so they cover what
    citations: you asked for.  A structure that already carries glycans, nucleic acids or
    citations: ligands may owe further CHARMM parameter-set citations.
    citations: ------------------------

The PDB2PQR and PROPKA reference strings are the citations those packages declare for themselves;
the rest match the bibliography of the pestifer paper.

Recovering either block from a log
----------------------------------

Both blocks are in the build log and in the diagnostics log:

.. code-block:: bash

    grep 'environment:' build.log     # what ran
    grep 'citations:'   build.log     # who to cite
