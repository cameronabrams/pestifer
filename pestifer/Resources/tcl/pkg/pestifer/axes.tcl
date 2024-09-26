# Author: Cameron F. Abrams <cfa22@drexel.edu>
package provide PestiferAxes 1.0

namespace eval ::PestiferAxes:: {
    namespace export *
}

proc PestiferAxes::axis { rlist } {
    # Given a list whose elements are parallel arrays of atom coordinates, one per protomer and therefore
    # assumed to be symmetry related, computes selected cross products of displacement vectors connecting
    # parallel coordinates and averages the resulting vectors to generate an approximation to the axis
    # around which the protomers are arranged.
    #
    #   2
    #   |
    #   |
    #   0-----1
    #
    #   v1=r0-r1
    #   v2=r0-r2
    #   c=v1 x v2
    #
    set npoints [llength $rlist]
    set accum [veczero]
    set count 0
    for {set i 0} {$i < $npoints} {incr i 3} {
        set r0 [lindex $rlist $i]
        set r1 [lindex $rlist [expr $i + 1]]
        set r2 [lindex $rlist [expr $i + 2]]
        set v1 [vecsub $r1 $r0]
        set v2 [vecsub $r2 $r0]
        set c [vecnorm [veccross $v1 $v2]]
        set accum [vecadd $accum $c]
        incr count
    }
    return [vecscale $accum [expr 1.0/$count]]

}

proc PestiferAxes::get_axis {mol frame chains selstr} {
    # using vector arithmetic, computes the approximate symmetry axis displayed by same atom selection in each
    # protomer listed in chains
    set acc [veczero]
    set count 0
    set R [list]
    # elements of R should be parallel arrays of symmetry-related atoms coordinates per chain/protomer
    foreach c $chains {
        set aa [atomselect $mol "chain $c and $selstr"]
        $aa frame $frame
        lappend R [$aa get {x y z}]
    }
    set natoms [llength [lindex $R 0]]
    for { set r 0 } { $r < $natoms } { incr r } {
        set chain [list]
        for { set c 0 } { $c < [llength $chains] } { incr c } {
            lappend chain [lindex [lindex $R $c] $r]
        }
        set this_c [axis $chain]
        set acc [vecadd $acc $this_c]
        incr count
    }

    return [vecscale $acc [expr 1.0/$count]]
}

proc PestiferAxes::get_center {mol frame chains selstr} {
    set acc [veczero]
    set count 0
    set R [list]
    # elements of R should be parallel arrays of symmetry-related atoms coordinates per chain/protomer
    foreach c $chains {
        set aa [atomselect $mol "chain $c and $selstr"]
        $aa frame $frame
        lappend R [$aa get {x y z}]
    }
    set natoms [llength [lindex $R 0]]
    for { set r 0 } { $r < $natoms } { incr r } {
        set chain [list]
        for { set c 0 } { $c < [llength $chains] } { incr c } {
            lappend chain [lindex [lindex $R $c] $r]
        }
        set this_c [vecadd {*}$chain]
        set acc [vecadd $acc $this_c]
        incr count [llength $chains]
    }

    return [vecscale $acc [expr 1.0/$count]]
}

proc PestiferAxes::get_tilt {molid frame chains1 selstr1 chains2 selstr2} {
    set trimer_axis [get_axis $molid $frame $chains1 $selstr1]
    set base_axis [get_axis $molid $frame $chains2 $selstr2]
    set tilt_dot [vecdot $trimer_axis $base_axis]
    set tilt_angle [expr acos($tilt_dot)]
    set pi [expr acos(-1)]
    return [expr $tilt_angle * 180/$pi]
}

proc PestiferAxes::get_tilt_z {molid frame chains1 selstr1} {
    set trimer_axis [get_axis $molid $frame $chains1 $selstr1]
    set tilt_dot [vecdot $trimer_axis { 0 0 -1 }]
    set tilt_angle [expr acos($tilt_dot)]
    set pi [expr acos(-1)]
    return [expr $tilt_angle * 180/$pi]
}

# proc PestiferAxes::draw_axis {molid v} {
# }