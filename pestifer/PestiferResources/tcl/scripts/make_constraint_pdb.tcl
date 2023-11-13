# Author: Cameron F. Abrams, <cfa22@drexel.edu>

set constrained_atoms_def ""
set pdb "UNSET"
set refpdb "ref.pdb"
set force_constant 200.0

for { set i 0 } { $i < [llength $argv] } { incr i } {
    if { [lindex $argv $i] == "-pdb"} {
       incr i
       set pdb [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-refpdb"} {
       incr i
       set refpdb [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-constrained_atoms_def"} {
       incr i
       set constrained_atoms_def [join [split [lindex $argv $i] ","] " "]
    }
    if { [lindex $argv $i] == "-force_constant"} {
       incr i
       set force_constant [lindex $argv $i]
    }
}

mol new $pdb

set a [atomselect top all]
$a set occupancy 0.0
set c [atomselect top "$constrained_atoms_def"]
$c set occupancy $force_constant
$a writepdb $refpdb

exit
