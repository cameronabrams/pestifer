# Author: Cameron F. Abrams, <cfa22@drexel.edu>

# VMD/psfgen script for creating a new psf/pdb pair for a
# bilayer created using the upper leaflet of one
# input bilayer and the lower leaflet of another

proc leaflet_apportionment { molid } {
   set whole [atomselect $molid all]
   set bilayer [atomselect $molid "lipid"]
   set com [measure center $bilayer weight mass]
   set com_z [lindex $com 2]
   $whole moveby [list 0 0 -$com_z]
   set residue_list [lsort -unique [atomselect $molid "residue"]]
   set residues(upper) [list]
   set residues(lower) [list]
   foreach residue $residue_list {
      set ressel [atomselect $molid "residue $residue"]
      set com [measure center $ressel weight mass]
      set com_z [lindex $com 2]
      if { $com_z > 0 } {
         lappend residues(upper) $residue
      } else {
         lappend residues(lower) $residue
      }
   }
   set residues(upper) [lsort -unique $residues(upper)]
   set residues(lower) [lsort -unique $residues(lower)]
   return $residues
}

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
set source(upper) [molinfo top get id]
# pbc readxst $xscA -molid $source(upper)
mol new $pdbB waitfor all
set source(lower) [molinfo top get id]
# pbc readxst $xscB -molid $source(lower)

foreach leaflet { upper lower } {
   set residues [leaflet_apportionment $source($leaflet)]
   set leaflet_sel [atomselect $source "residue $residues($leaflet)"]
   $leaflet_sel writepdb "${outbasename}-${leaflet}.pdb" 
}
exit
