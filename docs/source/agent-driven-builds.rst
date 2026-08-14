.. _agent_driven_builds:

Driving pestifer with an AI agent
=================================

It is increasingly common to install a tool and then hand its operation to a coding agent --
Claude Code, or something like it.  Pestifer is a reasonable fit for that: its input is a
declarative YAML file, its schema is introspectable, and its failures are mostly detectable
before anything runs.  But a few of its properties will defeat an agent that discovers them the
hard way, and this page is about those.

Install the skill first
-----------------------

.. code-block:: console

   $ pestifer setup-claude

This installs a skill describing everything below into ``~/.claude/skills/pestifer/``.  See
:ref:`subs_setup_claude`.  The rest of this page explains *why* the skill says what it says --
useful if you are adapting it, writing your own agent instructions, or driving pestifer from a
script rather than an agent.

Check, then build
-----------------

The single most important habit is to run :ref:`--check <subs_build>` before every build:

.. code-block:: console

   $ pestifer build config.yaml --check --json

An agent works by iterating, and the cost of one iteration governs how well it does.  Building
blind, a malformed config costs minutes of wall-clock before it reports a problem that was
knowable at parse time -- and an agent that has waited eleven minutes for an error is an agent
that has burned its context on a log.  ``--check`` makes that iteration cost about a second, and
its JSON form gives back the resolved toolchain, the task plan, and explicit ``errors`` and
``warnings`` arrays.

Schema errors are written to be self-correcting: a mistyped directive is reported together with
the list of directives that *are* valid at that point, which is usually enough for the next
attempt to succeed.

Prefer the examples over invention
----------------------------------

Pestifer ships 27 worked configurations, and they are the best available specification of what a
good config looks like for a given kind of system:

.. code-block:: console

   $ pestifer show-resources examples    # the catalog: ID, PDB ID, name, title
   $ pestifer fetch-example 7            # copy one here without building it

An agent adapting example 16 for a new membrane system will outperform the same agent composing
a membrane config from the reference documentation, because the example encodes decisions the
reference describes only one directive at a time.  For a structure unlike anything in the
catalog, :ref:`new-system --inspect <subs new-system>` writes a starting config annotated with
what is actually in the structure.

Know which commands expect a human
----------------------------------

Two subcommands prompt on standard input, and both will stall a non-interactive caller:

- ``pestifer new-system <id> --interactive`` walks through each detected structural feature.  The
  agent-appropriate form is ``--inspect``, which writes the same findings into the config as
  commented stubs to be edited.
- ``pestifer config-help`` defaults to an interactive query loop.  Pass ``--no-interactive``,
  with the directive path as positional arguments, to get a single answer and exit:

  .. code-block:: console

     $ pestifer config-help tasks psfgen mods --no-interactive

Builds are long-running jobs
----------------------------

A modest solvated protein -- 60,000 atoms, say -- takes on the order of ten minutes, and the
great majority of that is the final density equilibration.  Membrane systems and large
glycosylated assemblies take considerably longer.  A build is therefore a job to launch and poll,
not a command to wait on:

.. code-block:: console

   $ setsid nohup pestifer build config.yaml > build.log 2>&1 < /dev/null &

Each task logs when it begins and ends, so ``tail build.log`` is enough to follow progress.

Give every build its own directory
----------------------------------

Builds write hundreds of intermediate files under predictable names, so two builds sharing a
directory will interfere.  Create a fresh directory per build; ``--check`` warns when the working
directory already has other files in it.

Failures are resumable
----------------------

Pestifer records each completed task in ``.pestifer-manifest.json``, so an interrupted or failed
build does not have to start over:

.. code-block:: console

   $ pestifer build config.yaml --restart

This matters more for an agent than for a person, because an agent that must re-run a build from
scratch to test a one-line change will usually give up on the change instead.  Task specs are
hashed before execution, so editing the config resumes from the first task that actually changed.
See :ref:`subs_build` for ``--from`` and ``--fresh``.

The outputs are inside tarballs
-------------------------------

This surprises people and agents equally.  When a :ref:`terminate <subs_buildtasks_terminate>`
task carries a ``package:`` block, it sweeps the working directory once the package is written.
The run-ready files are then inside ``prod_<name>.tar.gz``, not loose on disk, and every
intermediate file is inside ``<basename>-artifacts.tar.gz``.  An agent looking for ``my_system.psf``
in the build directory will not find it and may conclude, wrongly, that the build failed.
