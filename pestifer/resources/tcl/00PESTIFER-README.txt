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

- pestifer
    Several custom packages for use in pestifer (Cameron F. Abrams).
    Those marked "in use" are loaded by a pestifer workflow; the rest
    ship and remain loadable by hand, but nothing in pestifer requires
    them and they are not maintained.
        - PestiferUtil        (in use)
            General utility procedures
        - PestiferCRot        (in use)
            Procedures for facilitating rotations around bonds
        - PestiferEnviron     (in use)
            Environment/setup procedures
        - PestiferIonize      (in use)
            Procedures for adding ions
        - PestiferAUTools
            Procedures for generating asymmetric units from full
            multimeric protein complexes
        - PestiferAxes
            Procedures for detecting axes of rotational symmetry
            in multimeric protein complexes
        - PestiferDeclash
            Procedures for rotating bonds to relieve steric clashes
        - PestiferGetLinks
            Procedures for automatically detecting LINKS from
            3D structures and atom types
        - PestiferMultimer
            Procedures for computing geometric features of
            multimers

Scripts:

- bilayer_embed.tcl
    Embeds a protein in a pre-built bilayer and generates the PSF
    for the combined system (sourced by make_membrane_system)
- bilayer_patch.tcl
    Builds a bilayer patch (sourced by make_membrane_system)
- domainswap.tcl
    A script that facilitates the generation of a NAMD input
    for an SMD run that conducts a domain-swap operation.  The
    domainswap task is no longer registered as a user-invocable
    task; the script is kept so the method is not lost.
- bilayer_orient.tcl
    Orients a transmembrane protein so that z is the membrane
    normal.  Superseded by the Python transform path; retained
    only as a test oracle (see below).

Third-party packages (not used by any pestifer workflow):

- la
    The Linear Algebra package from Hume Integration Software
- orient
    The Orient package from the VMD folks

Orientation moved into Python (pestifer.objs.rottrans, via
MakeMembraneSystemTask._orientation_align).  Their only consumer is
bilayer_orient.tcl, which is retained solely as the oracle for the
regression test that pins the Python path to the coordinates the Tcl
path produced.

Cameron F. Abrams, <cfa22@drexel.edu>

---------------------------------------------------------------------
