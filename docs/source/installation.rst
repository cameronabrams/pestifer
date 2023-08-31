Installation
============


Software Prequisites
--------------------

The following commands should be in your path: ``vmd``, ``charmrun``, and ``namd2``.

Installation
------------

To use Pestifer, install it from a cloned Github repository:

.. code-block:: console

   $ git clone git@github.com:cameronabrams/pestifer.git
   $ cd pestifer
   $ pip install -e .

Pestifer is under very active development.

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