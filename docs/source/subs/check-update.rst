.. _subs_check_update:

check-update
------------

The ``check-update`` subcommand reports whether a newer pestifer has been released on PyPI.

.. code-block:: console

   $ pestifer check-update
   Running:  3.24.2
   Latest:   3.25.0 (PyPI)

   A newer pestifer is available.
     upgrade:   pip install -U pestifer
     changelog: https://github.com/cameronabrams/pestifer/blob/main/CHANGELOG.md

pestifer also checks on its own, before running whatever command you asked for.  That check
writes to stderr beside the banner, so it appears in a redirected run as well -- a SLURM job, a
batch sweep, ``pestifer build system.yaml > run.log 2>&1`` -- which is where a stale install
would otherwise go unnoticed indefinitely.  It is bounded in three ways:

* **At most once a day.**  The answer is cached in ``~/.pestifer/update-check.json``.  A check
  that could not reach PyPI is cached the same way, so a cluster node with no route out pays one
  two-second timeout a day rather than one per build.
* **Never from a source checkout.**  A working tree is routinely ahead of the latest release.
* **Never fatal.**  Every failure is silent, including a bug in the check itself.  Nothing about
  a build depends on it.

``check-update`` ignores all three and asks every time.

.. note::

   Because the notice reaches redirected runs, a build log can differ between runs by this line:
   whether it appears depends on what PyPI said, which is outside the build.  Within one sweep
   the builds still agree with each other -- the first fetches and the rest read the same cached
   answer -- but two sweeps run on different days can differ.  **Where log comparability is the
   point, turn the check off rather than reasoning about it**: put
   ``export PESTIFER_NO_UPDATE_CHECK=1`` in the job script.

Turning it off
~~~~~~~~~~~~~~

Any one of these suppresses the automatic check:

.. code-block:: console

   $ pestifer check-update --disable      # persistent, recorded in ~/.pestifer/update-check.json
   $ export PESTIFER_NO_UPDATE_CHECK=1    # per shell or per job script
   $ pestifer --no-update-check build system.yaml   # one invocation

``pestifer check-update --no-disable`` turns the automatic check back on.  ``check-update`` itself
still works while the automatic check is disabled -- the setting governs only the check pestifer
makes on its own.

.. note::

   Nothing about a build depends on this check.  It contacts ``pypi.org`` and reads one version
   string; it sends nothing about you or your system, and a failure of any kind is silent.
