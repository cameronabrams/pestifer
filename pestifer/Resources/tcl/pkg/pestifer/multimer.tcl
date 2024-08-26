# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# procedures for computing geometrical features of multimers
#
#
package provide PestiferMultimer 1.0

namespace eval ::PestiferMultimer:: {
    namespace export *
}



# Compute and return geometric center of a collection of points
proc centroid { points } {
    set C {0 0 0}
    foreach p $points {
        set C [vecadd $p $C]
    }
    set C [vecscale $C [expr 1.0/[llength $points]]]
    return $C
}

# Given the list of vertices of a polygon (assumed),
# compute the sector angles formed by adjacent pairs of
# vertices subtended by the centroid
#
# Arguments:
# vertex_list : list of points (VMD vectors)
# Returns:
# vector of angles in radians
proc PestiferMultimer::sector_angles { vertex_list } {
    set C [centroid $vertex_list]
    set Cx [list]
    foreach p $vertex_list {
        set this_pC [vecsub $p $C]
        lappend Cx [vecnorm $this_pC]
    }
    set ang_list [list]
    set N [llength $vertex_list]
    for {set i 0} {$i < [expr $N-1]} {incr i} {
        set Cahat [lindex $Cx $i]
        set Cbhat [lindex $Cx [expr $i + 1]]
        set cos_ab [vecdot $Cahat $Cbhat]
        set ang_ab [expr acos($cos_ab)]
        lappend ang_list $ang_ab
    }
    set Cahat [lindex $Cx [expr $N - 1]]
    set Cbhat [lindex $Cx 0]
    set cos_ab [vecdot $Cahat $Cbhat]
    set ang_ab [expr acos($cos_ab)]
    lappend ang_list $ang_ab
    return $ang_list
}

proc PestiferMultimer::opening_angles { sel_list } {
    set pi [expr acos(-1.0)]
    set points [list]
    foreach sel $sel_list {
        lappend points [measure center $sel]
    }
    return [vecscale [PestiferMultimer::sector_angles $points] [expr 180.0/$pi]]
}

proc PestiferMultimer::opening_angles_trace { molid selstrings {ignoreH 1} } {
    set sel_list [list]
    foreach selstr $selstrings {
        if { $ignoreH == 1 } {
            set selstr "($selstr) and noh"
        }
        lappend sel_list [ atomselect $molid $selstr ]
    }
    set numframes [molinfo $molid get numframes]
    vmdcon -info "molid $molid has $numframes frames"
    set results [list]
    for {set i 0} {$i<$numframes} {incr i} {
        foreach s $sel_list {
            $s frame $i
        }
        set angles [PestiferMultimer::opening_angles $sel_list]
        lappend results $angles
    }
    return $results
}

proc PestiferMultimer::rotation_matrices_trace { molid selstrings {ignoreH 1} } {
    set sel_list [list]
    foreach selstr $selstrings {
        if { $ignoreH == 1 } {
            set selstr "($selstr) and noh"
        }
        lappend sel_list [ atomselect $molid $selstr ]
    }
    set numframes [molinfo $molid get numframes]
    vmdcon -info "molid $molid has $numframes frames"
    set results [list]
    for {set i 0} {$i<$numframes} {incr i} {
        foreach s $sel_list {
            $s frame $i
        }
        set rotmats [PestiferMultimer::rotation_matrices $sel_list]
        lappend results [regsub -all "\{|\}" $rotmats ""]
    }
    return $results
}

proc PestiferMultimer::rotation_matrices { sel_list } {
    set pi [expr acos(-1.0)]
    set rotmats [list]
    for { set i 0 } { $i < [llength $sel_list] } { incr i } {
        if { [expr $i + 1] == [llength $sel_list] } {
            set j 0
        } else {
            set j [expr $i + 1]
        }
        set selA [lindex $sel_list $i]
        set selB [lindex $sel_list $j]
        if { [$selA num] != [$selB num] } {
            vmdcon -error "[$selA text] ([$selA num] atoms) not conrugent with [$selB text] ([$selB num] atoms)"
        } else {
            set tmat [measure fit $selA $selB]
            set cosTheta [lindex [lindex $tmat 0] 0]
            set r33 [list [lrange [lindex $tmat 0] 0 2] [lrange [lindex $tmat 1] 0 2] [lrange [lindex $tmat 2] 0 2] ]
            set fl_r33 [regsub -all "\{|\}" $r33 ""]
            lappend rotmats $fl_r33
        }
    }
    return $rotmats
}

