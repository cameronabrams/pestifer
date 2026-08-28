.. _subs_buildtasks_validate:

validate
--------

A ``validate`` task allows you to perform validation tests on the current state files (PDB/PSF/XSC) to ensure, for example, that modifications requested in previous tasks have been successfully applied.

Validation tests are performed in VMD sessions managed by the task.  There are three types of validation tests available:

1.  :ref:`Attribute test <subs_buildtasks_validate_attribute_test>`: Allows you to query atom attributes in an atom selection of the current state for expected values.
2.  :ref:`Connection test <subs_buildtasks_validate_connection_test>`: Allows you to check on existence of inter-residue connections in an atom selection of the current state.
3.  :ref:`Residue test <subs_buildtasks_validate_residue_test>`: Allows you to check a count of atoms and a count of residues in an atom selection of the current state.

After all tests are performed, the task parses all VMD logs for results and generates a summary report. By default, if any test fails, the task raises an exception and the run is aborted.

.. _subs_buildtasks_validate_letters:

Test the molecule, not the letter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Chain and segment letters in a *built* system are assigned by pestifer, not inherited from the
input file.  A test written against a letter is therefore testing an output of the build, and
the most natural-looking test of all -- "I excluded chain B, so ``chain B`` should be empty" --
fails on a perfectly correct build.

The reason is that pestifer splits each input chain into one segment per segment type, and every
segment beyond the first takes the next **unused** letter.  Excluding a chain returns its letter
to that pool, where the next segment needing one picks it up.  Building PDB 8DX0 (two protein
chains, each with waters and a magnesium tagged with the same chain id) shows it directly:

.. code-block:: text

   no exclusion                        exclude: [chainID == 'B']

   A  149 residues  protein            A  149 residues  protein
   B  149 residues  protein            B    1 residue   MG      <-- chain A's magnesium
   C    1 residue   MG                 C  123 residues  water
   D    2 residues  MG
   E  123 residues  water
   F  105 residues  water

After the exclusion, ``segname B`` and ``chain B`` both name the magnesium ion lifted out of
chain A.  A test like

.. code-block:: yaml

   - residue_test:                # DON'T: passes or fails for the wrong reason
       name: chain B excluded
       selection: segname B
       measure: residue_count
       value: 0

reports ``FAIL ... count 1 (expected 0)`` on a build that excluded chain B exactly as asked.

Note that whether you get away with it depends on *which* letter you excluded, which is why this
trap is easy to miss: the pool is consumed in order, so a late letter is rarely reused, while an
early one is taken almost immediately.  Excluding chain ``P`` and then testing ``chain P`` (as
:ref:`example 5 <example hivprotease>` does) is safe by accident, not by design.

Write the test against something the build cannot rename -- the molecule itself:

.. code-block:: yaml

   - residue_test:                # DO: the count of what you excluded
       name: only one protomer remains
       selection: protein and name CA
       measure: residue_count
       value: 149

or assert on a residue that is unique to the excluded molecule (a distinctive ``resname``, a
sequence position that exists in only one protomer).  If you must refer to a chain, take the
letter from the built system rather than from the input file: the mapping is recorded in the
run's diagnostics log and in the psfgen script the task writes.

.. _subs_buildtasks_validate_attribute_test:

Attribute test
~~~~~~~~~~~~~~

A common use of an attribute test is to check that a mutation requested in a previous task has been successfully applied.  An attribute test requires four fields to be specified, along with one optional field:

1. **name**: a descriptive name for the test -- this is up to you.
2. **selection**: a string that represents VMD atomselect logic
3. **attribute**: the name of the atom attribute to query
4. **value**: the expected value of the atom attribute
5. **value_count**: the expected number of atoms in the selection with this attribute (defaults to 1)

For example, if you mutated resid 11 on chain A from threonine to cysteine, you can check that the residue name of that residue is in fact CYS:

.. code-block:: yaml

   tasks:
     - ... (prior tasks here)
     - validate:
         tests:
           - attribute_test:
               name: point mutation
               selection: protein and chain A and resid 11 and name CA
               attribute: resname
               value: CYS
     - ... (subsequent tasks here)

This test will generate and execute a VMD input script that looks like this:

.. code-block:: tcl

    (... read in psf and pdb of current state ...)
    set test_selection [atomselect top "{self.selection}"]
    set result [$test_selection get {self.attribute}]
    set count [llength [lsearch -nocase -exact -all $result {target}]]
    if {$count != 1} {
        vmdcon "FAIL resname has unexpected count $count of value CYS (expected 1) in selection protein and chain A and resid 11 and name CA"
    } else {
        vmdcon "PASS resname has expected count $count of value CYS (expected 1) in selection protein and chain A and resid 11 and name CA"
    }

The validate tasks manages all VMD logs it generates, so it will extract pass and fail information once all tests are run.

.. _subs_buildtasks_validate_connection_test:

Connection test
~~~~~~~~~~~~~~~

A connection test can be used to check for the presence or absence of an inter-residue bond.  This test requires three fields to be specified, along with one optional field:

1. **name**: a descriptive name for the test -- this is up to you.
2. **selection**: a string that represents VMD atomselect logic
3. **connection_type**: the type of connection to check for (one of ``interresidue``, ``disulfide``, or ``glycosylation``)
4. **connection_count**: the expected number of connections (defaults to 1)

There are three types of connections that can be checked:

1. **interresidue**: Checks for bonds between atoms in different residues.  The ``selection`` field should be atomselect logic that selects atoms in at most two distinct residues.
2. **disulfide**: Checks for disulfide bonds between cysteine residues.  The ``selection`` field should be atomselect logic that selects two residues that participate in the disulfide bond.  This test will check that the bond is in fact between to CYS residues and two SG atoms.
3. **glycosylation**: Checks for glycosylation sites on asparagine residues.  The ``selection`` field should be atomselect logic that selects for one asparagine residue.  This test will validate that the ND2 atom of the ASN residue is bonded to the C1 atom of a glycan monomer.  For example:

.. code-block:: yaml

   tasks:
     - ... (prior tasks here)
     - validate:
         tests:
        - connection_test:
            name: graft 61 chain A
            selection: protein and chain A and resid 61 
            connection_type: glycosylation
            connection_count: 1
     - ... (subsequent tasks follow)


.. _subs_buildtasks_validate_residue_test:

Residue test
~~~~~~~~~~~~

This test allows you to check the number of atoms and/or residues in an atomselection.
A residue test requires four fields to be specified:

1. **name**: a descriptive name for the test -- this is up to you.
2. **selection**: a string that represents VMD atomselect logic
3. **measure**: the type of measurement to perform (one of ``atom_count`` or ``residue_count``)
4. **value**: the expected value of the measurement

For example, if you excluded a residue by residue name, you might want to test that no such residues remain in the system:

.. code-block:: yaml

   tasks:
     - ... (prior tasks here)
     - validate:
         tests:
           - residue_test:
               name: no POPC residues
               selection: resname POPC
               measure: residue_count
               value: 0
     - ... (subsequent tasks here)
