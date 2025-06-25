.. _subs_runtasks_pdb2pqr:

pdb2pqr
-------

A ``pdb2pqr`` task can be used to determine the protonation states of residues based on their calculated pKa values and a specified pH. The only subdirective that is required is ``pH``, which specifies the pH at which the pKa values are calculated. The task will use the PROPKA algorithm to calculate the pKa values and assign ionization states to residues in the task's incoming structure.  It is best to use this task after a ``psfgen`` task has processed the external structure file but before any solvation or membrane tasks.

Example 
+++++++

.. code-block:: yaml

    tasks:
      - psfgen:
          ...
      - pdb2pqr:
          pH: 7.0


How it works
++++++++++++

This task first takes the incoming structure and prepares an input for ``pdb2pqr`` by stripping all hydrogens, naming all histidine residues as ``HIS`` and all water residues as ``HOH``.  Since this happens after a ``psfgen`` task, the input structure's names are already in CHARMM format.  The task then runs ``pdb2pqr`` with the PROPKA algorithm to calculate pKa values and assign ionization states to residues based on the specified pH. The output of the command is a PDB file that indicates which residues had their protonation states updated.  This task also captures the summary table showing all residue pKa values.

The task then prepares for a new invocation of ``psfgen`` that will apply all the protonation states as patches.  This ``psfgen`` will use the incoming PSF/PDB combination, *NOT* the input or output of ``pdb2pqr``.  This task can handle the following patches:

- ``NNEU``: Neutral N-terminus
- ``CNEU``: Neutral C-terminus
- ``HS2``: Convert histidine protonated on ND1 (``HSD``) to histidine protonated on NE2 (``HSE``) (``HSD`` is the default assignment to all histidines)
- ``HSPP``: Convert histidine protonated on ND1 (``HSD``) to histidine protonated on both NE2 and ND1 (``HSP``)
- ``TYRO``: Deprotonated tyrosine
- ``SERD``: Deprotonated serine
- ``ASPP``: Protonate aspartic acid (``ASP``) on OD1 (``ASPP``)
- ``GLUP``: Protonate glutamic acid (``GLU``) on OE1 (``GLUP``)
- ``LSN`` : Deprotonated lysine (atom ``HZ3`` is removed)
- ``RN2`` : Deprotonated arginine (atom ``HH12`` is removed)

The task then runs ``psfgen`` to apply the patches to the incoming PSF/PDB combination, and it outputs a new PSF/PDB combination that has the protonation states applied.

References
++++++++++

* Very Fast Empirical Prediction and Rationalization of Protein pKa Values.
  Hui Li, Andrew D. Robertson and Jan H. Jensen. PROTEINS: Structure, Function,
  and Bioinformatics. 61:704-721 `DOI 10.1002/prot.20660 <https://doi.org/10.1002/prot.20660>`_ (2005)

* Very Fast Prediction and Rationalization of pKa Values for Protein-Ligand
  Complexes.  Delphine C. Bas, David M. Rogers and Jan H. Jensen.  PROTEINS:
  Structure, Function, and Bioinformatics 73:765-783 `DOI 10.1002/prot.22102 <https://doi.org/10.1002/prot.22102>`_ (2008)

*  PROPKA3: Consistent Treatment of Internal and Surface Residues in Empirical
   pKa predictions.  Mats H.M. Olsson, Chresten R. Sondergard, Michal Rostkowski,
   and Jan H. Jensen.  Journal of Chemical Theory and Computation, 7(2):525-537
   `DOI 10.1021/ct100578z <https://doi.org/10.1021/ct100578z>`_ (2011)

*  Improved Treatment of Ligands and Coupling Effects in Empirical Calculation
   and Rationalization of pKa Values.  Chresten R. Sondergaard, Mats H.M. Olsson,
   Michal Rostkowski, and Jan H. Jensen.  Journal of Chemical Theory and
   Computation, 7(7):2284-2295 `DOI 10.1021/ct200133y <https://doi.org/10.1021/ct200133y>`_ (2011)
