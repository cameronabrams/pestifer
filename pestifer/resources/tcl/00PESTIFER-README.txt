---------------------------------------------------------------------
TcL scripts and packages for pestifer

Base-level scripts:

- macros.tcl
    Defines several convenient VMD atomselect macros (created 
    automatically from the ycleptic base.yaml via pestifer inittcl)
- vmdrc.tcl
    A general VMD startup script that defines a lot of
    useful stuff

Packages:

- la
    The Linear Algebra package from Hume Integration Software
- orient
    The Orient package from the VMD folks
- pestifer
    Several custom packages for use in pestifer (Cameron F. Abrams)
        - PestiferAUTools
            Procedures for generating asymmetric units from full
            multimeric protein complexes
        - PestiferAxes
            Procedures for detecting axes of rotational symmetry
            in multimeric protein complexes
        - PestiferCRot
            Procedures for facilitating rotations around bonds
        - PestiferDeclash
            Procedures for rotating bonds to relieve steric clashes
        - PestiferGetLinks
            Procedures for automatically detecting LINKS from
            3D structures and atom types
        - PestiferMultimer
            Procedures for computing geometric features of 
            multimers
        - PestiferUtil
            General utility procedures

Scripts:

- domainswap.tcl
    A script that facilitates the generation of a NAMD input
    for an SMD run that conducts a domain-swap operation
- loop_closure.tcl
    Facilitates the ligation of fictitious C- and N-termini
    in protein loop-building
- memb.tcl
    Facilitates generation of PSF file for a PDB generated
    by packmol (usually for membrane-embedding)
- tg.tcl
    Facilitates the TopoGromacs package (experimental)
- tmd_prep.tcl (unused)

Cameron F. Abrams, <cfa22@drexel.edu>

---------------------------------------------------------------------