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
    if { [lindex $argv $i] == "-p"} {
       incr i
       set patchfile [lindex $argv $i]
    }
}

set outpsf ${outbasename}.psf
set outpdb ${outbasename}.pdb

readpsf $psf pdb $pdb

source $patchfile

regenerate angles dihedrals
guesscoord

writepsf $outpsf
writepdb $outpdb
exit
