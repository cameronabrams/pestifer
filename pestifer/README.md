# pestifer

This python script automates the process of building an input script for `psfgen` using as much information as possible in a published PDB or mmCIF file.  

Cameron F Abrams -- cfa22@drexel.edu

2020-2023


### What it does

`pestifer.py` parses most of the PDB/mmCIF records necessary to build a complete ``psfgen`` input script.  It uses SEQRES records to build chains, the REMARK 465 records to indicate missing residues, and the LINK and SSBOND records to include the proper bond patches.  It will generate replica protomers using BIOMT records to generate multimers.  It can generate the proper commands to mutate residues and cleave chains.

`pestifer.py` is still a work in progress and has not been widely tested on random PDB entries. 