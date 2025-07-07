.. _example env 4tvp:

Example 8: Closed, PGT122/35O22-Liganded HIV-1 BG505 Env SOSIP.664 Trimer (ligands removed)
-------------------------------------------------------------------------------------------

`PDB ID: 4tvp <https://www.rcsb.org/structure/4tvp>`_ is a structure of the HIV-1 Env ectodomain trimer in a closed conformation, with the PGT122 and 35O22 antibodies bound.  Here we prepare a system with just the ectodomain, omitting the ligands.

Many structures in the RCSB are only available in mmCIF format, rather than the older, outdated PDB format; `4tvp <https://www.rcsb.org/structure/4TVP>`_ is one example.  Fortunately, because ``pestifer`` uses the ``pidibble`` `package <https://pypi.org/project/pidibble/>`_, this is quite easy to handle in pestifer:

.. literalinclude:: ../../../pestifer/resources/examples/hiv-sosip-env-ectodomain2.yaml
    :language: yaml
