make-resi-database
------------------

With this subcommand, you can build PSF and PDB files of single molecules suitable for system-building with, for example, packmol.  An extensive library of lipids (all those defined in the main lipids topology file for CHARMM36FF) is already included with pestifer.  If you want to make more, for example, some that are only topologically defined in some of the CHARMMFF stream files, you can do this:

.. code-block:: console

    $ pestifer make-resi-database <stream-name> <resi-name>

where ``stream-name`` is the name of the CHARMM36 molecule type (e.g., lipid) and ``resi-name`` is a RESI defined in the appropriate stream file.

``make-resi-database`` has a lot of options; the defaults provided have been tested on all lipid RESIs in the CHARMM36FF and seem adequate.

``make-resi-database`` will generate a directory called ``data/<stream-name>/<resi-name>/``, under which you will find 10 unique PDB samples of the molecule in various configurations.

A ``bilayer`` task can be told which of these configs to use when packing a system.  To use your own locally generated PDBs from ``make-resi-database``, you need to declare the location of your local database in your config file.  To do this, you need to assign the value ``data`` to the ``pdb_depot`` variable in the ``paths`` section.  You build a database outside your build directory to collect molecules you want to use in multiple different builds, and just refer to that relative to your build directory.

A directory tree recognizable as a PDB depot has subdirectories that are known CHARMM36FF streams, each of which can contain just naked PDB files with names of the form ``<RESI>.pdb``, or subdirectories for each RESI with the configurational samples mentioned above. 

