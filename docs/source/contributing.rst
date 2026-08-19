.. _contributing:

Contributing: running and writing tests
=======================================

Most of what makes pestifer hard to test is that it drives other programs.  VMD, psfgen and NAMD
are subprocesses, their outputs are files, and the failures that matter are usually *silent* — a
system that builds, runs, and is quietly wrong.  The test suite is shaped around that, and this
page is about the shape.

.. important::

   **Running ``pytest`` on its own gives you a partial suite, and says so only in passing.**
   Slow tests are skipped unless you ask for them, and tests needing VMD/NAMD skip themselves when
   those programs are absent.  A green run on a machine without the toolchain has not exercised a
   single build.  See :ref:`contributing_markers`.

Setting up
----------

The project is managed with `uv <https://docs.astral.sh/uv/>`_.  Nothing needs installing
first — ``uv run`` resolves the environment on demand:

.. code-block:: bash

   uv run --extra test pytest tests/unit -q

To work against an editable install (what most development looks like):

.. code-block:: bash

   uv sync --extra test

Building a system also needs ``vmd``, ``namd3`` and ``catdcd`` on ``PATH``.  Without them the
pure-Python tests still run; the build tests skip.

The two tiers
-------------

``tests/unit``
    Everything that can be asserted without running a build: schema handling, parsers,
    coordinate and geometry code, the convergence criterion, config resolution, the CLI's
    plumbing.  Fast (~7 minutes for ~1200 tests) and hermetic.  A few of these still invoke
    VMD or NAMD for a single short step and are marked accordingly.

``tests/integration``
    Builds real systems and asserts on their geometry and topology.  Ten tests, about two
    minutes for the gated set.  These are the only tests that can catch a wrong *structure*,
    as opposed to a wrong function return.

.. _contributing_markers:

Markers, and what they hide
---------------------------

Three markers gate what actually runs.  They are registered in ``tests/conftest.py``.

``slow``
    Requires ``--runslow``.  Every integration test carries it.  Without the flag they are
    silently skipped.

``needs_tools``
    Auto-skips when ``vmd`` or ``namd3`` is not on ``PATH``, rather than failing.  This is what
    lets the pure-Python core run anywhere, including CI — and it is also why a green suite on a
    toolchain-free machine proves less than it appears to.

``expensive``
    Minutes-long even by integration standards (currently only the membrane-patch build, ~3.5
    minutes).  Excluded from the release gate; run deliberately.

So:

.. code-block:: bash

   pytest tests/unit -q                                  # the fast core
   pytest tests/integration -q --runslow                 # builds, including the expensive one
   pytest tests/integration -q --runslow -m 'not expensive'   # what the release gate runs
   pytest tests/integration -q --runslow -m expensive         # only the expensive ones

Why CI runs only ``tests/unit``
-------------------------------

GitHub-hosted runners have no VMD and no NAMD.  Adding ``tests/integration`` to
``.github/workflows/tests.yaml`` would not fail — it would **skip**, silently, and report green.
That is worse than not running them, because it looks like coverage.

The integration tests therefore run in two places instead:

* on a developer machine that has the toolchain, and
* in ``scripts/release.sh``, which runs the gated set before tagging and refuses to release if it
  fails.  Set ``SKIP_INTEGRATION_TESTS=1`` to bypass that deliberately.

If you want them in CI, the answer is a self-hosted runner with VMD and NAMD installed, not a
workflow edit.

Writing an integration test
---------------------------

Two things make these tests worth their runtime.

**Assert on invariants, not snapshots.**  A stored reference file breaks whenever anything
unrelated changes and tells you nothing about *why*.  A physical invariant survives refactors and
states what is actually required.  The existing tests are the model: a BIOMT image must superpose
onto its asymmetric-unit chain (it is a rigid copy) *and* sit far from it in place (the transform
moved it); a modelled loop's junction must be a real peptide bond; every intra-glycan bond must be
a plausible covalent distance.

**Assert that the test exercised its path.**  A test that silently checks nothing is worse than no
test.  ``test_glycan_bonds`` fails outright if no branched glycan is found, because an unbranched
one would not exercise the reordering bug it guards.  ``test_loop_closure`` fails if it walked no
backbone bonds.  Add that guard whenever a test's reach depends on the input.

The shared battery
~~~~~~~~~~~~~~~~~~

``tests/integration/helpers.py`` provides ``assert_psf_sane``, which every build should pass
regardless of what it was exercising:

.. code-block:: python

   from tests.integration.helpers import assert_psf_sane, parse_psf_pdb

   assert_psf_sane(psf, pdb, context='my new test')

It checks that the PSF and PDB describe the same system, that no bond is duplicated, self-bonded
or implausibly long, that nothing sits at the origin, that the total charge is integral, and that
no two non-bonded atoms occupy the same space.  Each check corresponds to a bug that has actually
shipped here.

Two parameters matter:

``unminimized=True``
    For a build that stops at psfgen with no relaxation.  psfgen *guesses* coordinates for
    hydrogens it adds, and an unrelaxed guess can legitimately sit 3 Å from its parent or on top
    of a neighbour.  This mode skips the contact check and applies the bond-length check to heavy
    atoms only — which is what mis-wiring corrupts, so the check that matters is retained.

``max_bond``
    Defaults to 2.5 Å, not the ~1.5 Å of a covalent single bond: a disulfide S–S is ~2.05 Å and is
    correct.  A test checking a narrower population — intra-glycan bonds, at ~1.4 Å — should pass a
    tighter value.

The helpers parse raw text rather than using pestifer's own readers, deliberately.  A test that
checks the writer with the writer's reader cannot catch a bug in the representation they share,
and several of the bugs these tests guard are exactly that shape.

Making sure a test can fail
---------------------------

Before trusting a new test, break the thing it checks and confirm it goes red.  This is not
ceremony — it has caught real problems here twice:

* the generated-LaTeX compile test passed while a version string with an underscore made the
  fragment fail to compile, because the brace-balance and escaping assertions only covered the
  case that had been thought of;
* an equilibration-loop stub keyed its simulated crash to a chunk number that never advanced, so
  the retry path it was meant to exercise ran forever instead.

Coverage
--------

.. code-block:: bash

   uv run --extra test --with pytest-cov pytest tests/unit -q --cov=pestifer --cov-report=term

Note that this measures the unit tier only; the integration tests exercise substantially more of
the task layer than the numbers suggest.  Coverage is a floor, not a goal: the modules with the
highest line coverage here are the ones that were easiest to test, which is not the same as the
ones most likely to be wrong.

Building the docs
-----------------

.. code-block:: bash

   UV_OVERRIDE=docs/overrides.txt uv run --with-requirements docs/requirements.txt \
       sphinx-build -b html docs/source /tmp/pestifer-docs

The override is needed because ``pdb2pqr`` declares a stale ``docutils<0.18`` pin it never uses.

Two things to know when editing the schema:

* every new key in ``pestifer/schema/base.yaml`` needs a ``text:`` entry, or the generated
  configuration reference fails to build;
* run ``yclept check-spec pestifer/schema/base.yaml`` afterwards — ycleptic silently ignores some
  malformed declarations rather than complaining.
