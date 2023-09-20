Installation
============


Software Prequisites
--------------------

The commands ``vmd``, ``charmrun``, and ``namd2`` should be in your path.  If you choose to package the results using the ``topogromacs`` plugin, the ``gmx`` executable should also be in your path.

Installation
------------

To use Pestifer, install it from PyPI:

.. code-block:: console

   $ pip install pestifer

Pestifer is under very active development.  To get the latest version for PyPI, just update:

.. code-block:: console

   $ pip install -U pestifer


.. If you use conda/anaconda, we recommended that you create a separate Python environment running ``HTPolyNet``:

.. .. code-block:: console

..     $ conda create --name mol-env python
..     $ conda activate mol-env

.. Once this environment is created and activated, you can install both ``ambertools`` and ``HTPolyNet`` from ``conda-forge``:

.. .. code-block:: console

..     $ conda install -c conda-forge ambertools
..     $ conda install -c conda-forge htpolynet

.. If you are not a conda user, you can install ``HTPolyNet`` from PyPI.

.. .. code-block:: console

..     $ pip install htpolynet