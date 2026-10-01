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

pestifer also checks on its own, before running whatever command you asked for, but it does so
under three restrictions that make this command worth having:

* **Only on a terminal.**  A redirected run -- ``pestifer build system.yaml > run.log 2>&1``, the
  usual form in a cluster job or a batch sweep -- makes no network request and prints nothing.
  This keeps build logs identical between runs, so a log diff still means something.
* **At most once a day.**  The answer is cached in ``~/.pestifer/update-check.json``.  A check
  that could not reach PyPI is cached in the same way, so a machine with no route out does not
  pay a timeout on every invocation.
* **Never from a source checkout.**  A working tree is routinely ahead of the latest release.

``check-update`` ignores all three and asks every time.

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
