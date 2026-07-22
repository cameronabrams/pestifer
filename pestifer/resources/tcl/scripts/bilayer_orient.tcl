# Author: Cameron F. Abrams, <cfa22@drexel.edu>

## VMD for orienting a transmembrant protein in a coordinate system in which z is membrane normal


package require Orient
namespace import ::Orient::*

set scriptname bilayer_orient

set pdb "";  # psf of protein-only system
set psf "";  # pdb of protein-only system
set z_head_group ""; # atomselection whose center of mass defines one end of the protein axis
set z_tail_group ""; # atomselection whose center of mass defines the other end of the protein axis
set outbasename "oriented"

for { set i 0 } { $i < [llength $argv] } { incr i } {
   if { [lindex $argv $i] == "-psf"} {
      incr i
      set psf [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-pdb"} {
      incr i
      set pdb [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-z_head_group"} {
      incr i
      set z_head_group [deprotect_str_arg [lindex $argv $i]]
   }
   if { [lindex $argv $i] == "-z_tail_group"} {
      incr i
      set z_tail_group [deprotect_str_arg [lindex $argv $i]]
   }
   if { [lindex $argv $i] == "-o"} {
      incr i
      set outbasename [lindex $argv $i]
   }
}

# check for required arguments
if { $psf == "" } {
   puts "Error: -psf argument required"
   exit
}
if { $pdb == "" } {
   puts "Error: -pdb argument required"
   exit
}
if { $z_head_group == "" && !$no_orient } {
   puts "Error: -z_head_group argument required if -no_orient is not indicated"
   exit
}
if { $z_tail_group == "" && !$no_orient } {
   puts "Error: -z_tail_group argument required if -no_orient is not indicated"
   exit
}

# load the protein system
mol new $psf
mol addfile $pdb waitfor all
set protein [molinfo top get id]
set pro_sel [atomselect $protein "all"]

vmdcon -info "orienting protein axis to bilayer normal"
set head [atomselect $protein "$z_head_group"]
set tail [atomselect $protein "$z_tail_group"]
set head_com [measure center $head weight mass]
set tail_com [measure center $tail weight mass]
set pro_axis [vecnorm [vecsub $head_com $tail_com]]
set A [orient $pro_sel $pro_axis {0 0 1}]
$pro_sel move $A

# write output pdb
$pro_sel writepdb "${outbasename}.pdb"
