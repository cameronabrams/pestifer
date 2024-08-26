# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Use VMD to detect glycan-based linkages and SS bonds; output them in YAML format for a pestifer input file
#
package provide PestiferGetLinks 1.0

namespace eval ::PestiferGetLinks:: {
    namespace export *
}

proc PestiferGetLinks::getssbonds { molid } {
    set sg [atomselect $molid "name SG"]
    set sg_ndx [$sg get index]
    set sg_ch [$sg get chain]
    set sg_resid [$sg get resid]
    set sg_insertion [$sg get insertion]
    set sg_bondlist [$sg getbonds]
    set ssbonds [list]
    foreach s $sg_ndx ch $sg_ch resid $sg_resid ins $sg_insertion neighbors $sg_bondlist {
        if {[llength $ins] > 0} {
            vmdcon -info "C1 at residue $resid in chain $ch has insertion $ins"
            set at1 "${ch}_${resid}${ins}"
        } else {
            set at1 "${ch}_${resid}"
        }
        if {[llength $neighbors] == 0} {
            vmdcon -info "Resid $resid of chain $ch has no ligands on its SG"
        }
        foreach ndx $neighbors {
            set nid [lindex [[atomselect $molid "index $ndx"] get {chain resid insertion name}] 0]
            set nch [lindex $nid 0]
            set nre [lindex $nid 1]
            set nins [lindex $nid 2]
            set nan [lindex $nid 3]
            if {[llength $nins] > 0} {
                set at2 "${nch}_${nre}${nins}"
            } else {
                set at2 "${nch}_${nre}"
            }
            if { $nan == "SG" } {
                if { $ch == $nch } {
                    if { $resid < $nre } {
                        lappend ssbonds "${at1}-${at2}"
                    }
                } elseif { $ch < $nch } {
                    lappend ssbonds "${at1}-${at2}"
                }
            }
        }
    }
    return $ssbonds
}

proc PestiferGetLinks::getlinks { molid } {
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
            vmdcon -info "C1 at residue $resid in chain $ch has insertion $ins"
            set at1 "${ch}_${resid}${ins}_C1"
        } else {
            set at1 "${ch}_${resid}_C1"
        }

        if {[llength $neighbors] == 0} {
            vmdcon -info "Resid $resid of chain $ch has no ligands on its C1"
        }

        foreach ndx $neighbors {
            set nid [lindex [[atomselect $molid "index $ndx"] get {chain resid insertion name}] 0]
            # puts "$nid"
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