# Author: Cameron F. Abrams, <cfa22@drexel.edu>

for { set i 0 } { $i < [llength $argv] } { incr i } {
    if { [lindex $argv $i] == "-pdb"} {
       incr i
       set pdb [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-psf"} {
       incr i
       set psf [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-o"} {
       incr i
       set outbasename [lindex $argv $i]
    }
}

mol new $psf
mol addfile $pdb waitfor all
check_pierced_rings 0 6 1.5
check_pierced_rings 0 5 1.5
exit
