# Author: Cameron F. Abrams, <cfa22@drexel.edu>

# This VMD script takes the PSF and namdbin files and
# generates a PDB format file from the coordinates
for { set i 0 } { $i < [llength $argv] } { incr i } {
    if { [lindex $argv $i] == "-pdb"} {
       incr i
       set pdb [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-psf"} {
       incr i
       set psf [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-coor"} {
       incr i
       set namdbin [lindex $argv $i]
    }
}
mol new $psf
mol addfile $pdb waitfor all
set a [atomselect top all]
$a writenamdbin $namdbin
exit
