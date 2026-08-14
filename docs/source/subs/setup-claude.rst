.. _subs_setup_claude:

setup-claude
------------

The ``setup-claude`` subcommand installs pestifer's bundled `Claude Code
<https://claude.com/claude-code>`_ skill, so that an AI agent asked to prepare a system with
pestifer knows how to do it properly.

.. code-block:: console

   $ pestifer setup-claude

This copies the skill that ships inside the pestifer package to
``~/.claude/skills/pestifer/SKILL.md``, where Claude Code discovers it automatically on its
next session in any directory.  Nothing else is modified.

.. code-block:: console

   $ pestifer setup-claude --skills-dir ./.claude/skills   # scope it to one project
   $ pestifer setup-claude --force                         # overwrite an installed copy

Without ``--force`` an existing installed skill is left alone.  Re-run the subcommand after
upgrading pestifer to pick up that version's skill.

What the skill contains
=======================

The skill is not a copy of this documentation.  It carries the operating knowledge that makes
the difference between an agent that drives pestifer competently and one that flounders:

- **Check before building.**  ``pestifer build config.yaml --check`` costs a second and catches
  what would otherwise surface minutes into a run.  See :ref:`subs_build`.
- **Start from a worked example** rather than composing YAML from scratch, and find the right
  one with ``pestifer show-resources examples``.
- **Which subcommands need a terminal.**  ``new-system --interactive`` and ``config-help``
  without ``--no-interactive`` both prompt on standard input and will stall an agent; the skill
  directs it to ``new-system --inspect`` and ``config-help --no-interactive`` instead.
- **Builds take minutes to hours**, so they must be launched in the background and polled,
  never waited on.
- **One clean directory per build**, and how to recover a failed one from the run manifest.
- **Where the outputs actually are** -- inside the package tarball, not loose in the directory.

The reasoning behind each of these, and the workflow they add up to, is in
:ref:`agent_driven_builds`.
