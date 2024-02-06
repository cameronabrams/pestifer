# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# procedures for computing geometrical features of multimers
#
#
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
proc sector_angles { vertex_list } {
    set C [centroid $vertex_list]
    set Cx [list]
    foreach p $vertex_list {
        set this_pC [vecsub $p $C]
        lappend Cx [vecscale $this_pC [expr 1.0/[veclength $this_pC]]]
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

proc opening_angles { sel_list } {
    set pi [expr acos(-1.0)]
    set points [list]
    foreach sel $sel_list {
        lappend points [measure center $sel]
    }
    return [vecscale [sector_angles $points] [expr 180.0/$pi]]
}

proc opening_angles_trace { molid refchains refresids } {
    set sel_list [list]
    foreach c $refchains {
        lappend sel_list [atomselect $molid "chain $c and resid $refresids"]
    }
    set numframes [molinfo $molid get numframes]
    puts "molid $molid has $numframes frames"
    set results [list]
    for {set i 0} {$i<$numframes} {incr i} {
        foreach s $sel_list {
            $s frame $i
        }
        set angles [opening_angles $sel_list]
        lappend results $angles
    }
    return $results
}
