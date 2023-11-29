# Author: Cameron F. Abrams, <cfa22@drexel.edu>
# coords.tcl 
# Generic VMD script for altering coordinates
#
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
    if { [lindex $argv $i] == "-c"} {
       incr i
       set coordfile [lindex $argv $i]
    }
}

set outpdb ${outbasename}.pdb

mol load psf $psf pdb $pdb

source $coordfile

[atomselect top all] writepdb $outpdb
exit
