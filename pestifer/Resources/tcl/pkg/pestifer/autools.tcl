# Author: Cameron F. Abrams <cfa22@drexel.edu>
package provide PestiferAUTools 1.0

namespace eval ::PestiferAUTools:: {
    namespace export *
}

proc format_tmat { tmat biomtnum } {
    set result [list]
    set rownum 1
    foreach row [lrange $tmat 0 end-1] {
        set rowstr "REMARK 350   BIOMT[format %1d $rownum][format %4d $biomtnum]"
        foreach elem [lrange $row 0 end-1] {
            set rowstr "$rowstr[format %10.6f $elem]"
        }
        set le [lindex $row end]
        set rowstr "$rowstr[format %16.6f $le]"
        lappend result $rowstr
        incr rownum
    }
    return [join $result \n]
}

# make_au_from_aligning
#
# This procedure will generate a new PDB file containing a asymmetric
# unit and appropriate BIOMT operations to build a complex based
# on aligning the *base* chains onto congruent chains of any other
# subunits.  For example, if you have a molecule with chains A, B, and C
# where there is quasi symmetry among the three, you can build a
# a new pdb containing just A and the BIOMT transformations that will
# allow for replicating A and putting it in the place of B and C via
# alignment.
#
# Arguments:
#  - molid : molecule ID
#  - protomer_chains: list of chains by subunit; first entry signifies
#    the base chain(s).
#  - resids: string of resid range description for an atomselection
#    to use as an alignment basis
#  - outpub: string name of output pdb file to be generated
#
# Example:
#   Suppose you have a system (molid 0) with three protomers: 
#   chains A, B (protomer 1); chains C, D (2); and chains E, F (3).
#   Suppose there is *nearly* three-fold symmetry, and you want to 
#   generate a *strictly* three-fold symmetric structure by aligning 
#   replicas of chains A and B respectively onto chains C and D and 
#   chains E and F.  Suppose further that the resids you want to use to 
#   do the alignment are 1-100.  In a VMD session where the molecule
#   is loaded, issue this command
#
#   make_au_from_aligning 0 {{A B} {C D} {E F}} "1-100" "my_au.pdb"
#
#   and the file "my_au.pdb" will contain the BIOMT transformations
#   in REMARK 350 records and all atoms for chain A, B.
#
proc PestiferAUTools::make_au_from_aligning { molid protomer_chains resids outpdb } {
    set tmats [list]
    set c1 [lindex $protomer_chains 0]
    set savepdb [atomselect $molid "(protein or glycan) and chain $c1"]
    set r_p1 [atomselect $molid "(protein or glycan) and chain $c1 and resid $resids"]
    for {set i 1} { $i < [llength $protomer_chains] } { incr i } {
        set c [lindex $protomer_chains $i]
        set p [atomselect $molid "(protein or glycan) and chain $c and resid $resids"]
        set tmat [measure fit $r_p1 $p]
        lappend tmats $tmat
    }
    set fp [open "tmp.pdb" "w"]
    puts $fp "REMARK 350"
    puts $fp "REMARK 350 BIOMOLECULE: 1"
    set chains [join $c1 ", "]
    puts $fp "REMARK 350 APPLY THE FOLLOWING TO CHAINS: $chains"
    puts $fp "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000"
    puts $fp "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000"
    puts $fp "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000"
    set i 2
    foreach tmat $tmats {
        puts $fp "[format_tmat $tmat $i]"
        incr i 
    }
    puts $fp "REMARK 350"
    close $fp
    $savepdb writepdb "autmp.pdb"
    set out [open $outpdb "w"]
    fconfigure $out -translation binary
    foreach fname {"tmp.pdb" "autmp.pdb"} {
        set in [open $fname]
        fcopy $in $out
        close $in
    }
    close $out
    file delete tmp.pdb autmp.pdb
}

proc PestiferAUTools::overlay_protomers { molid protomer_chains resids } {
    set c1 [lindex $protomer_chains 0]
    set r_p1 [atomselect $molid "(protein or glycan) and chain $c1 and resid $resids"]
    set rmsd [list]
    for {set i 1} { $i < [llength $protomer_chains] } { incr i } {
        set c [lindex $protomer_chains $i]
        set p [atomselect $molid "(protein or glycan) and chain $c"]
        set pga [atomselect $molid "(protein or glycan) and chain $c and resid $resids"]
        $p move [measure fit $pga $r_p1]
        lappend rmsd [measure rmsd $pga $r_pi]
    }
    return $rmsd
}