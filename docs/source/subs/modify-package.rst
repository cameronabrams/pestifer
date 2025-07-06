.. _subs modify-package:

modify-package
---------------

The subcommand is meant for development use only. It allows modifications to the **source package**, including adding and deleting examples and updating atomselect macros.

This subcommand will only work on a full source repository, so if you want to use it, you will need to clone the repository from GitHub:

.. code-block:: bash

    git clone git@github.com:cameronabrams/pestifer.git
    cd pestifer
    pip install -e . # so that you can tell `pestifer` to modify itself

I also recommend creating a Git branch for your modifications, so you can easily revert them if needed:

.. code-block:: bash

    git checkout -b my-modifications

Adding and Deleting Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

I developed most of the examples by iteration, so I would start with a simple YAML config file and then modify it to add new features or test new functionality.  