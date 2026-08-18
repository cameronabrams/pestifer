.. _subs_report_methods:

report-methods
--------------

Drafts a Methods section, and its bibliography, from what a build actually did.  The analogue of
``gmx report-methods``, with more to say: a pestifer build knows not only the final system but how
it was produced — which entry it started from, what was modelled in, which force field, what
equilibration actually ran, under which software versions, from which seed.

.. code-block:: console

   $ cd my-build-directory
   $ pestifer report-methods
   Drafted a Methods section from 1 run record:
     ./methods.tex
     ./methods.bib

``methods.tex`` is a **fragment** — no preamble, ``\input``-able straight into a manuscript.
``methods.bib`` carries every citation the build owes.  Add ``--format md`` if you are not writing
LaTeX.

Where the facts come from
=========================

Every completed build writes :ref:`run-record.json <build_provenance>`, a machine-readable record
of what it did.  ``report-methods`` reads that; it does not re-derive anything from your config.

That distinction is the point.  ``density_equilibrate`` and ``membrane_equilibrate`` run until
their convergence criterion is met, so **the number of steps they take is decided at run time and
appears nowhere in your input file** — the same config legitimately produces different step counts
on different days.  A Methods paragraph written from the config would state numbers the run never
produced.

An aborted build writes no record, so there is nothing to report on; that is deliberate.

Replicas
========

Pass several run directories and they are described as one set:

.. code-block:: console

   $ pestifer report-methods example-01/rep-01 example-01/rep-02 example-01/rep-03

.. code-block:: text

   3 independent replicas were prepared from the same input, differing only in their
   random number seeds (27021972, 27021973, 27021974).

Facts every replica agrees on are stated once.  A fact they disagree on — replicas equilibrate to
slightly different cells — is reported as a range rather than dropped:

.. code-block:: text

   Final box dimensions across replicas ranged from 55.945 $\times$ 49.378 $\times$ 49.998
   to 56.231 $\times$ 49.63 $\times$ 50.253 \AA.

See :ref:`build_provenance` for how to produce replicas in the first place.

What it will not do
===================

The output is a **draft**, and says so in a banner you have to delete deliberately.  Three rules
constrain it:

* **Every sentence traces to a recorded fact from the run that produced it.**  Nothing is inferred
  from your config where the run could have diverged from it.
* **Pestifer never describes what it did not do.**  It prepared a system; it did not run your
  production simulation, choose the ensemble for your science, or analyse anything.  Those arrive
  as visible ``\TODO`` markers, never as prose.
* **Nothing is stated that the record does not support.**  A build that stopped at its step ceiling
  is reported as having stopped at its step ceiling, with a ``\TODO`` telling you to say so or
  re-run — it is never described as having converged.

The register is **deliberately machine-generated**: short declarative sentences, no connective
flourish, no rationale.  That is a design decision.  Fluent prose reads as authoritative and
invites being pasted unread; text that visibly came from a tool invites the editing it needs.  You
supply the scientific voice — this supplies the facts you would otherwise transcribe by hand and
occasionally get wrong.

For the same reason, it is **not** run automatically at the end of a build.  ``terminate`` records
the facts; turning them into prose is a separate, explicit act.

The bibliography
================

Software citations are curated BibTeX entries maintained inside pestifer, so they are correct by
construction and need no network.

Citations for the coordinates you started from are emitted as ``@misc`` entries carrying the DOI
and nothing else.  This is deliberate: a PDB-format file stores author names and titles in upper
case, and the original capitalization is not recoverable — ``A.B.MCDERMOTT`` cannot be turned back
into ``McDermott, A.B.``, and guessing sentence case from an ALL-CAPS title mangles things like
``HIV-1``.  Your reference manager resolves a DOI correctly; invented metadata would look complete
and be wrong.

Using the output
================

The fragment expects a ``\TODO`` command in your preamble:

.. code-block:: latex

   \newcommand{\TODO}[1]{\textbf{[TODO: #1]}}
   ...
   \input{methods}
   \bibliographystyle{plain}
   \bibliography{methods}

Options
=======

.. code-block:: bash

  runs                  run directories or run-record.json files (default: the current
                        directory); give several to describe replicas as one set
  --format {tex,bib,md} output format, repeatable (default: tex and bib)
  --output-dir DIR      where to write the files (default: the current directory)
