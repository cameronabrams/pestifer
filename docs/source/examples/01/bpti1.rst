.. _example bpti1:

Example 1: Bovine Pancreatic Trypsin Inhibitor (BPTI)
-----------------------------------------------------

The ``psfgen`` `user manual <https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/ug.pdf>`_ is a necessary resource for learning how to use ``psfgen`` to generate PDB and PSF input files for NAMD.  A simple example in that manual is a solvation of bovine pancreatic trypsin inhibitor (BPTI) starting from its PDB coordinates (`PDB ID 6pti <https://www.rcsb.org/structure/6PTI>`_).  ``pestifer`` can reproduce this solvation via the input YAML-format configuration shown below:

.. literalinclude:: ../../../../pestifer/resources/examples/01/inputs/bpti1.yaml
    :language: yaml

.. task-table:: ../../../../pestifer/resources/examples/01/inputs/bpti1.yaml


You can check the :ref:`config_ref` for a complete reference to Pestifer config files.

There are three ways to launch this build.  **They are alternatives — pick one, and
run it in a clean directory.**  All three run the same configuration; they differ only
in how the config file gets into that directory, and whether you get a chance to edit
it first.

**Option 1 — build the packaged example in one step**

.. code-block:: bash

   $ pestifer build-example 1

``build-example`` copies the example's config file into the current directory as
``bpti1.yaml`` — exactly what you see above — and immediately builds it.  Nothing to
edit, nothing to name.

**Option 2 — supply your own config file**

Paste the YAML above into a local file and build that instead.  This is the route you
will use for your own systems:

.. code-block:: bash

   $ pestifer build myconfig.yaml

**Option 3 — fetch the config first, then build it**

``fetch-example`` copies the example's config file *without* building, so you can read
or edit it before committing to a run:

.. code-block:: bash

   $ pestifer fetch-example 1
   $ ls
   bpti1.yaml
   $ pestifer build bpti1

(If there is no extension on the argument of ``build``, pestifer assumes one of
``.yaml``, ``.yml``, or ``.ym``.)

``bpti1.yaml`` is a YAML-format text file, and the keywords (of course) have particular meanings.  This is also an example of a "minimal" configuration file; ``pestifer`` has many more controls that can be set in a configuration file than are shown here.  Here, this configuration file contains two topmost directives: ``title`` and :ref:`config_ref tasks`.  The value of ``title`` is the string ``Bovine Pancreatic Trypsin Inhibitor (BPTI)`` and the value of ``tasks`` is a *list*.  Each element in the list of tasks is itself a directive describing a task, and ``pestifer`` in general executes tasks in the order they appear in the ``tasks`` list.

Pestifer's configuration schema documents itself: every task and subdirective can be
explored from the command line with ``config-help``, and the same content is published
in the :ref:`config_ref`.  If you want to know what else a task accepts, see
:ref:`exploring_config_schema`.

After the ``psfgen`` and ``validate`` tasks we declare an ``md`` task, with the
subdirective ``ensemble`` set to ``minimize`` and nothing else listed, so this task runs
an energy minimization under ``namd3`` with defaults for everything else.

After the first ``md`` task is the ``solvate`` task.  Notice that it has no subdirectives; in this case default values are used for any subdirectives.  After this task has finished, we have a run-ready *nonequilibrated* system.  Equilibration here is three tasks: another minimization, a short NVT ``md`` task to bring the system to temperature, and then a single :ref:`density_equilibrate <subs_buildtasks_density_equilibrate>` task that settles the box density under NPT.

That last task replaces what used to be written by hand as a ladder of progressively longer NPT ``md`` tasks (``nsteps: 200``, then 400, then 800, and so on).  Why such runs have to be broken into pieces at all is worth knowing, because it is a property of NAMD rather than of the chemistry: after solvation the initial water box is slightly too sparse, so under pressure control the cell shrinks.  NAMD fixes its spatial decomposition at the *start* of each run, and if a cell dimension shrinks past the patch margin partway through, the run dies with *"periodic cell has become too small for the patch grid"*.  Restarting resets that decomposition, so the equilibration has to proceed in restarts.  ``density_equilibrate`` does that for you, and sizes each chunk from the shrink rate it has just measured rather than from a guessed schedule.

**Deciding when the system is equilibrated is the task's job, not yours.**  It runs until the box density is *measurably* stationary and then stops on its own.  The test is not "the curve looks flat": density fluctuates intrinsically under NPT and those fluctuations are correlated in time, so the task estimates the integrated autocorrelation time of the series, uses it to compute an honest standard error of the mean, and requires that error to be small compared with the drift tolerance before it will even test for drift.  Only then does it ask whether the trend across the averaging window is below tolerance, and it must answer yes several times in a row.  A small, noisy box therefore runs longer than a large quiet one without anyone tuning a step count.

The task writes its own convergence plot, which shows the reasoning rather than only the trace:

.. figure:: density-convergence.png

    Box density vs. timestep for the BPTI system, as reported by ``density_equilibrate``.  The grey band at the left is the burn-in discarded as the barostat transient; the green band is the trailing window the convergence test is evaluated over; the red points are per-chunk means of that window; and the marker at the right is the step at which the run decided it was finished.  The final density is stated on the plot.

The ``mdplot`` task then plots the same quantity across *all* the post-solvation dynamics -- the minimization and the NVT warm-up as well as the NPT stage -- which is the wider view rather than the convergence argument:

.. figure:: solvated-density.png

    Density vs. timestep for the BPTI system post-solvation, across every simulation stage.  The dashed vertical rules mark the boundaries between stages, so the jump as the barostat engages reads as a change of ensemble rather than as a physical event.

Between the two, nothing about the equilibration has to be assumed: the first figure shows the criterion that ended the run, and the second shows that criterion in the context of everything that came before it.

Finally, we see a ``terminate`` task, whose main role is to generate some informative output and to provide a set of NAMD input files (PSF, PDB, xsc, coor, and vel) that all have a common base file name.  The ``package`` subdirective creates a tarball named for its own ``basename`` (here ``prod_6pti.tar.gz``) containing all input files necessary to execute a NAMD run, ready for transfer to the HPC resource of your choice.  Inside the tarball the files sit under a single directory also named for the package ``basename`` (``prod_6pti/``); the state files keep the terminate ``basename`` (``my_6pti``) while the generated NAMD config uses the package ``basename``.  The ``namd`` attribute of the ``package`` subdirective lets you specify any NAMD configuration options for the production config file; here we simply request the default parameters for an NPT run.

The ``terminate`` task also archives the build's other working files into a second tarball, named from the terminate ``basename`` (here ``my_6pti-artifacts.tar.gz``); the ``artifacts`` attribute can override that name.  Plot images are deliberately *not* swept into it -- they stay in the ``mdplots/`` subdirectory of the run, where they are easy to look at without unpacking anything.

Listing the contents of the package tarball:

.. code-block:: bash

   $ tar ztf prod_6pti.tar.gz
    prod_6pti/my_6pti.psf
    prod_6pti/my_6pti.pdb
    prod_6pti/my_6pti.coor
    prod_6pti/my_6pti.xsc
    prod_6pti/my_6pti.vel
    prod_6pti/prod_6pti.namd
    prod_6pti/my_6pti_minimal.prm

That is the whole production system: the five state files, the NAMD config, and a single parameter file.  ``my_6pti_minimal.prm`` is generated by the ``terminate`` task, which merges every parameter and stream file the build drew on and keeps only the records matching atom types actually present in the PSF.  You do not have to carry the CHARMM distribution to the machine that runs the job, and there is no search path to get right there -- the tarball is self-contained.

The archive tarball contains all intermediate files used in the build but which are not necessary for production MD runs.  These files can be useful for debugging or for understanding the build process.

Intermediate files are named by a fixed convention (``CC-MT-SSS_TASKNAME.ext``) that is
worth knowing when you read a build directory; see :ref:`file name conventions`.

.. raw:: html

        <div class="autogen-footer">
            <p>Example author: Cameron F. Abrams&nbsp;&nbsp;&nbsp;Contact: <a href="mailto:cfa22@drexel.edu">cfa22@drexel.edu</a></p>
        </div>
