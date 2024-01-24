# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Use VMD to detect glycan-based linkages
#

proc getlinks { molid } {
    set gc1 [atomselect $molid "glycan and name C1"]
    set gc1_ndx [$gc1 get index]
    set gc1_ch [$gc1 get chain]
    set gc1_resid [$gc1 get resid]
    set gc1_insertion [$gc1 get insertion]
    set gc1bl [$gc1 getbonds]
    # C1_RRR1_A1-C2_RRR2_A2
    set links [list]
    foreach c1 $gc1_ndx ch $gc1_ch resid $gc1_resid ins $gc1_insertion neighbors $gc1bl {
        if {[llength $ins] > 0} {
            puts "C1 at residue $resid in chain $ch has insertion $ins"
            set at1 "${ch}_${resid}${ins}_C1"
        } else {
            set at1 "${ch}_${resid}_C1"
        }

        foreach ndx $neighbors {
            set nid [lindex [[atomselect $molid "index $ndx"] get {chain resid insertion name}] 0]
            puts "$nid"
            set nch [lindex $nid 0]
            set nre [lindex $nid 1]
            set nins [lindex $nid 2]
            set nan [lindex $nid 3]
            if { $resid != $nre } {
                if {[llength $nins] > 0} {
                    set at2 "${nch}_${nre}${nins}_${nan}"
                } else {
                    set at2 "${nch}_${nre}_${nan}"
                }
                set lstr "${at2}-$at1"
                lappend links $lstr
            }
        }
    }
    return $links
}