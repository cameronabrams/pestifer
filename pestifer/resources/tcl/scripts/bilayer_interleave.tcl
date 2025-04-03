# Author: Cameron F. Abrams, <cfa22@drexel.edu>

# VMD script for creating a new psf/pdb pair for a complete, 
# membrane-embedded protein pdb and the original protein psf/pdb

set scriptname bilayer_interleave

set pdbA ""
set pdbB ""
set outbasename "bilayer-interleaved"

for { set i 0 } { $i < [llength $argv] } { incr i } {
    if { [lindex $argv $i] == "-pdbA"} {
       incr i
       set pdbA [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-pdbB"} {
       incr i
       set pdbB [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-o"} {
       incr i
       set outbasename [lindex $argv $i]
    }
}
mol new $pdbA waitfor all
set upper_source [molinfo top get id]
mol new $pdbB waitfor all
set lower_source [molinfo top get id]

foreach source { $upper_source $lower_source } {

}