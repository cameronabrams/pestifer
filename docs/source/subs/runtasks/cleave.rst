.. _subs_runtasks_cleave:

cleave 
------

A ``cleave`` task is used to cleave one or more chains in a system.  A ``cleave`` task requires a single ``sites`` subdirective which contains a list of cleavage site desigations.  The format for a cleavage site designation is:

``<chain>:<resi1>-<resi2>``

``<chain>`` is the chain ID, ``<resi1>`` and ``<resi2>`` are the residue numbers (including optionally appended insertion codes) of the two residues that form the peptide bond.  For example, ``A:12-15`` would create a cleavage between residue 12 and residue 15 in chain A. 

When a cleavage occurs, all residues in the parent chain N-terminal to the cleavage site retain their original chain ID, while those residues C-terminal to the cleavage site are assigned the first new available chain ID.

An example of a ``cleave`` task is shown in :ref:`example 13`.
