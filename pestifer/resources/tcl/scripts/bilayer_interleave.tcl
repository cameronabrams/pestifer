# Author: Cameron F. Abrams, <cfa22@drexel.edu>

# VMD/psfgen script for creating a new psf/pdb pair for a
# bilayer created using the upper leaflet of one
# input bilayer and the lower leaflet of another

package require PestiferEnviron 1.0
namespace import ::PestiferEnviron::*

package require pbctools

set scriptname bilayer_interleave

set pdbA ""
set xscA ""
set pdbB ""
set xscB ""
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
   if { [lindex $argv $i] == "-xscA"} {
      incr i
      set xscA [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-xscB"} {
      incr i
      set xscB [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-o"} {
      incr i
      set outbasename [lindex $argv $i]
   }
}
mol new $pdbA waitfor all
set upper_source [molinfo top get id]
pbc readxst $xscA -molid $upper_source
mol new $pdbB waitfor all
set lower_source [molinfo top get id]
pbc readxst $xscB -molid $lower_source

foreach source { $upper_source $lower_source } criterion { "z>0" "z<0" } {
   set whole [atomselect $source all]
   set bilayer [atomselect $source "lipid"]
   set com [measure center $bilayer weight mass]
   set com_z [lindex $com 2]
   $whole moveby [list 0 0 -$com_z]
   set leaflet [atomselect $source "same molecule as ($criterion)"]
   
}