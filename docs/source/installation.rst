Installation
============


Software Prequisites
--------------------

The commands ``vmd``, ``charmrun``, and ``namd2`` should be in your path.  If you choose to package the results using the ``topogromacs`` plugin, the ``gmx`` executable should also be in your path.  To build membrane systems using ``pestifer``'s ``packmol_memgen`` integration, you must have the ``ambertools`` package (v. 23.6 or better) installed (preferable via ``conda``).

Installation
------------

To use Pestifer, install it from PyPI:

.. code-block:: console

   $ pip install pestifer

Pestifer is under very active development.  To get the latest version from PyPI, just update:

.. code-block:: console

   $ pip install -U pestifer

.. To use packmol-memgen, you must have ambertools installed in a suitable conda environment.  You need not run pestifer in this specific environment, but pestifer will detect this environment and use it to run packmol-memgen.

.. Before you can use packmol-memgen, you must edit this file:
.. <conda-root>/envs/<env-name>/lib/python<python-ver>/site-packages/packmol_memgen/lib/pdbremix/v3numpy.py
.. and change instances of ``np.float`` to ``np.float64``
.. here, <conda-root> is your conda root directory (mine is ~/anaconda3)
.. <env-name> is the name of the environment in which you installed ambertools
.. <python-ver> is the python version in that environment
.. if you try to use packmol-memgen from pestifer and pestifer detects that the v3numpy.py file is not patched, it will exit with an error message.

.. If you use conda/anaconda, we recommended that you create a separate Python environment running ``pestifer``:

.. .. code-block:: console

..     $ conda create --name mol-env python
..     $ conda activate mol-env

.. Once this environment is created and activated, you can install both ``ambertools`` from ``conda-forge``:

.. .. code-block:: console

..     $ conda install -c conda-forge ambertools
..     $ conda install -c conda-forge htpolynet

.. If you are not a conda user, you can install ``pestifer`` from PyPI.

.. .. code-block:: console

..     $ pip install pestifer