# Author: Cameron F. Abrams, <cfa22@drexel.edu>

# A VMD/TcL script to perform a domain-swapping TMD simulation using NAMD colvars

proc repl_w_avg { L } {
    set sum 0.0
    foreach l $L {
        set sum [expr $sum + $l]
    }
    set avg [expr $sum/[llength $L]]
    set nL [list]
    foreach l $L {
        lappend nL $avg
    }
    return $nL
}

set swap_domain_def ""
set anchor_domain_def ""
set chain_swap_pairs [list]
set total_complex_def "protein or glycan"
set psf "UNSET"
set pdb "UNSET"
set refpdb "ref.pdb"
set cvfile "cv.inp"
set force_constant 200.0
set target_numsteps 10000
# puts "$argv"

for { set i 0 } { $i < [llength $argv] } { incr i } {
    if { [lindex $argv $i] == "-pdb"} {
       incr i
       set pdb [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-psf"} {
       incr i
       set psf [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-refpdb"} {
       incr i
       set refpdb [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-cv"} {
       incr i
       set cvfile [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-anchor_domain_def"} {
       incr i
       set anchor_domain_def [join [split [lindex $argv $i] ","] " "]
    }
    if { [lindex $argv $i] == "-total_complex_def"} {
       incr i
       set total_complex_def [join [split [lindex $argv $i] ","] " "]
    }
    if { [lindex $argv $i] == "-swap_domain_def"} {
       incr i
       set swap_domain_def [join [split [lindex $argv $i] ","] " "]
    }
    if { [lindex $argv $i] == "-chain_swap_pairs"} {
        incr i
        set otok [split [lindex $argv $i] ":"]
        foreach o $otok {
            set itok [split $o ","]
            lappend chain_swap_pairs $itok
        }
    }
    if { [lindex $argv $i] == "-force_constant"} {
       incr i
       set force_constant [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-target_numsteps"} {
       incr i
       set target_numsteps [lindex $argv $i]
    }
}

mol new $psf
mol addfile $pdb

set ref_coms [list]
[atomselect top all] set occupancy 0
set com [measure center [atomselect top "$total_complex_def"]]
set idx 1

set lefts [list]
set rights [list]
foreach elem $chain_swap_pairs {
    lappend lefts [lindex $elem 0]
    lappend rights [lindex $elem 1]
}

foreach c $lefts {
    set swap_domain_sel [atomselect top "chain $c and $swap_domain_def"]
    $swap_domain_sel set occupancy $idx
    set swap_sel($c) $swap_domain_sel
    set swap_idx($c) $idx
    set idx [expr $idx + 1]

    set anchor_domain_sel [atomselect top "chain $c and $anchor_domain_def"]
    $anchor_domain_sel set occupancy $idx
    set anchor_idx($c) $idx
    set anchor_sel($c) $anchor_domain_sel
    set idx [expr $idx + 1]
}

[atomselect top all] writepdb $refpdb

set cvfp [open $cvfile "w"]
# colvars for anchors
set anchor_cvs [list]
set anchor_centers [list]
foreach c $lefts {
    set anchor_xyz [measure center $anchor_sel($c)]
    lappend anchor_cvs "anchor${c}"
    lappend anchor_centers 0.0
    puts $cvfp "colvar {"
    puts $cvfp "    name anchor${c}"
    puts $cvfp "    distance {"
    puts $cvfp "        group1 {"
    puts $cvfp "            atomsFile $refpdb"
    puts $cvfp "            atomsCol O"
    puts $cvfp "            atomsColValue $anchor_idx($c)"
    puts $cvfp "        }"
    puts $cvfp "        group2 {"
    puts $cvfp "            dummyAtom ([lindex $anchor_xyz 0],[lindex $anchor_xyz 1],[lindex $anchor_xyz 2])"
    puts $cvfp "        }"
    puts $cvfp "    }"
    puts $cvfp "}"
}

# colvars for pivot arms
set pivotarm_cvs [list]
set pivotarm_centers [list]
foreach c $lefts {
    lappend pivotarm_cvs "pivotarm${c}"
    set sloc [measure center $swap_sel($c)]
    lappend pivotarm_centers [veclength [vecsub $sloc $com]]
    puts $cvfp "colvar {"
    puts $cvfp "    name pivotarm${c}"
    puts $cvfp "    distance {"
    puts $cvfp "        group1 {"
    puts $cvfp "            atomsFile $refpdb"
    puts $cvfp "            atomsCol O"
    puts $cvfp "            atomsColValue $swap_idx($c)"
    puts $cvfp "        }"
    puts $cvfp "        group2 {"
    puts $cvfp "            dummyAtom ([lindex $com 0],[lindex $com 1],[lindex $com 2])"
    puts $cvfp "        }"
    puts $cvfp "    }"
    puts $cvfp "}"
}
set pivotarm_centers [repl_w_avg $pivotarm_centers]

# colvars for spacers
set spacer_cvs [list]
set spacer_centers [list]
foreach c $lefts d $rights {
    lappend spacer_cvs "spacer${c}to${d}"
    set cloc [measure center $swap_sel($c)]
    set dloc [measure center $swap_sel($d)]
    lappend spacer_centers [veclength [vecsub $cloc $dloc]]
    puts $cvfp "colvar {"
    puts $cvfp "    name spacer${c}to${d}"
    puts $cvfp "    distance {"
    puts $cvfp "        group1 {"
    puts $cvfp "            atomsFile $refpdb"
    puts $cvfp "            atomsCol O"
    puts $cvfp "            atomsColValue $swap_idx($c)"
    puts $cvfp "        }"
    puts $cvfp "        group2 {"
    puts $cvfp "            atomsFile $refpdb"
    puts $cvfp "            atomsCol O"
    puts $cvfp "            atomsColValue $swap_idx($d)"
    puts $cvfp "        }"
    puts $cvfp "    }"
    puts $cvfp "}"
}
set spacer_centers [repl_w_avg $spacer_centers]

# colvars for alignments
set align_cvs [list]
set align_centers [list]
set align_targets [list]
foreach c $lefts d $rights {
    lappend align_cvs "align${c}to${d}"
    lappend align_centers [measure rmsd $swap_sel($c) $swap_sel($d)]
    lappend align_targets 0.0
    puts $cvfp "colvar {"
    puts $cvfp "    name align${c}to${d}"
    puts $cvfp "    rmsd {"
    puts $cvfp "        atoms {"
    puts $cvfp "            atomsFile $refpdb"
    puts $cvfp "            atomsCol O"
    puts $cvfp "            atomsColValue $swap_idx($c)"
    puts $cvfp "            rotateReference off"
    puts $cvfp "            centerReference off"
    puts $cvfp "        }"
    puts $cvfp "        refPositionsFile $refpdb"
    puts $cvfp "        refPositionsCol O"
    puts $cvfp "        refPositionsColValue $swap_idx($d)"
    puts $cvfp "    }"
    puts $cvfp "}"
}
set align_centers [repl_w_avg $align_centers]

# biases for anchors
puts $cvfp "harmonic {"
puts $cvfp "    colvars $anchor_cvs"
puts $cvfp "    centers $anchor_centers"
puts $cvfp "    forceConstant 200"
puts $cvfp "}"
# biases for pivotarms
puts $cvfp "harmonic {"
puts $cvfp "    colvars $pivotarm_cvs"
puts $cvfp "    centers $pivotarm_centers"
puts $cvfp "    forceConstant 200"
puts $cvfp "}"
# biases for spacers
puts $cvfp "harmonic {"
puts $cvfp "    colvars $spacer_cvs"
puts $cvfp "    centers $spacer_centers"
puts $cvfp "    forceConstant 200"
puts $cvfp "}"
# biases for alignments
puts $cvfp "harmonic {"
puts $cvfp "    colvars $align_cvs"
puts $cvfp "    centers $align_centers"
puts $cvfp "    forceConstant 200"
puts $cvfp "    targetCenters $align_targets"
puts $cvfp "    targetNumSteps $target_numsteps"
puts $cvfp "}"
close $cvfp

exit
