.. _example sars cov2 spike ba2:

Example 15: Fully Glycosylated, Closed SARS-CoV-2 Omicron BA.2 Variant Spike
----------------------------------------------------------------------------

This example highlights the use of Pestifer to build a fully glycosylated SARS-CoV2 Spike protein (BA.2 strain) using grafted glycans and cleaving at the furin cleavage sites.  This build is based on the `PDB entry 7xix <https://www.rcsb.org/structure/7XIX>`_, which contains a spike protein in the closed conformation.  The PDB file contains glycans, but they are not fully resolved, so we graft glycans from prototypical structures.

The glycans taken from prototypical structures are the following:

- `PDB ID 2byh <https://www.rcsb.org/structure/2WAH>`_ chain C is a "high-mannose" glycan with 9 mannoses; its full name is alpha-D-mannopyranose-(1-2)-alpha-D-mannopyranose-(1-6)-[alpha-D-mannopyranose-(1-3)]alpha-D-mannopyranose-(1-6)-[alpha-D-mannopyranose-(1-2)-alpha-D-mannopyranose-(1-3)]beta-D-mannopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose
- `PDB ID 4b7i <https://www.rcsb.org/structure/4B7I>`_ chain C is an "intermediate" glycan with 5 mannoses and a fucose; its full name is alpha-D-mannopyranose-(1-3)-[alpha-D-mannopyranose-(1-6)]alpha-D-mannopyranose-(1-6)-[alpha-D-mannopyranose-(1-3)]beta-D-mannopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose-(1-4)-[alpha-L-fucopyranose-(1-6)]2-acetamido-2-deoxy-beta-D-glucopyranose
- `PDB ID 4byh <https://www.rcsb.org/structure/4BYH>`_ chain C is a "complex" glycan; its full name is N-acetyl-alpha-neuraminic acid-(2-6)-beta-D-galactopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose-(1-2)-alpha-D-mannopyranose-(1-6)-[2-acetamido-2-deoxy-beta-D-glucopyranose-(1-2)-alpha-D-mannopyranose-(1-3)]beta-D-mannopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose-(1-4)-[alpha-L-fucopyranose-(1-6)]2-acetamido-2-deoxy-beta-D-glucopyranose

.. figure:: https://cdn.rcsb.org/images/carbohydrates/wa/2wah/2wah_SNFG_2.svg
            :alt: 2byh chain C glycan
            :width: 400px
            :align: center

            High-mannose glycan from PDB ID 2byh chain C.  Green circles denote mannoses, either α or β, and blue circles denote N-acetylglucosamines.

.. figure:: https://cdn.rcsb.org/images/carbohydrates/b7/4b7i/4b7i_SNFG_2.svg
            :alt: 4b7i chain C glycan
            :width: 400px
            :align: center

            Intermediate glycan from PDB ID 4b7i chain C.  Green circles denote mannoses, either α or β, blue circles denote N-acetylglucosamines, and the red triangle denotes fucose.

.. figure:: https://cdn.rcsb.org/images/carbohydrates/by/4byh/4byh_SNFG_2.svg
            :alt: 4byh chain C glycan
            :width: 400px
            :align: center

            Complex glycan from PDB ID 4byh chain C.  Green circles denote mannoses, either α or β, blue circles denote N-acetylglucosamines, red triangle denotes fucose, yellow circles denote galactose, and the purple diamond denotes sialic acid.

The script below shows the use of ``graft`` modifications to include the glycans.  The glycan assignments (i.e., which asparagines have high-mannose, intermediate, and complex glycans) are taken from `Watanabe et al. (2020) <https://doi.org/10.1126/science.abb9983>`_.  The commented-out integer labels on each graft directive indicate the residue numbers in the PDB file to which the glycans are grafted.  

The ``cleave`` task is used to cleave each protomer at its furin cleavage site (residue 685).

.. literalinclude:: ../../../pestifer/resources/examples/sars-cov2-S-BA2/sars-cov2-S-BA2.yaml
    :language: yaml


Note the various syntax used in the ``graft`` directives.  For example:

.. code-block:: yaml

    graft:
      - A_1304:4b7i,C_1-8 # 66

This indicates that the glycan from PDB ID 4b7i chain C, residues 1 to 8, is grafted onto resid 1304 of chain A on the spike.  That resid is not the asparagine at position 61; it is the primary NAG attached to Asn61.  Residue 1 of chain C of 4b71 is also a primary NAG, so the graft operation aligns the entire glycan such that its primary NAG aligns on the primary NAG already resolved in the spike's structure.  That NAG is deleted and then the glycan from 4b7i is attached directly from the C1 atom of the primary NAG to the ND2 atom of Asn61.

.. code-block:: yaml
    
    graft:
      - D_1-2:4byh,C_1#2-10 # 616

In contrast, this indicates that the glycan from PDB ID 4byh chain C, residues 1 and 2, is grafted onto resid 1 and 2 of chain D on the spike. Chain D happens to be just the two NAGs at Asn 616 on one protomer. The ``1#2`` notation means to take resid 1 and 2 from chain C of 4byh and use them together as an alignment basis before grafting.

A future release of pestifer will allow for more transparent specification of glycans.

The prototypical glycans in the same PDB structure but not used here are:

- PDB ID 2wah chain D; beta-D-mannopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose
- PDB ID 4byh chain D; beta-D-galactopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose-(1-2)-alpha-D-mannopyranose-(1-3)-[beta-D-galactopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose-(1-2)-alpha-D-mannopyranose-(1-6)]beta-D-mannopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose-(1-4)-[alpha-L-fucopyranose-(1-6)]2-acetamido-2-deoxy-beta-D-glucopyranose

.. figure:: https://cdn.rcsb.org/images/carbohydrates/wa/2wah/2wah_SNFG_3.svg
            :alt: 2wah chain D glycan
            :width: 400px
            :align: center

            Prototypical glycan from PDB ID 2wah chain D.  Green circles denote mannoses, either α or β, and blue circles denote N-acetylglucosamines.

.. figure:: https://cdn.rcsb.org/images/carbohydrates/by/4byh/4byh_SNFG_3.svg
            :alt: 4byh chain D glycan
            :width: 400px
            :align: center

            Prototypical glycan from PDB ID 4byh chain D.  Green circles denote mannoses, either α or β, blue circles denote N-acetylglucosamines, red triangle denotes fucose, and yellow circles denote galactoses.

