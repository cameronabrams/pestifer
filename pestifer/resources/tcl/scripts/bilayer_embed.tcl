# Author: Cameron F. Abrams, <cfa22@drexel.edu>

# VMD/psfgen script for creating a new psf/pdb pair for a complete, 
# membrane-embedded protein pdb from an intact membrane system
# and the original protein psf/pdb

# if referenced using the Psfgen scriptwriter, all common psfgen
# pre-build commands are invoked automatically
# and the script is run in the context of the psfgen package

# package require PestiferEnviron 1.0
# namespace import ::PestiferEnviron::*

package require Orient
namespace import ::Orient::*

set scriptname bilayer_embed

set pdb "";  # psf of protein-only system
set psf "";  # pdb of protein-only system
set bilayer_pdb ""; # the membrane-only system
set bilayer_psf ""; # the membrane-only system
set bilayer_xsc ""; # the membrane-only system
set z_head_group ""
set z_tail_group ""
set z_ref_group ""
set z_value 0.0
set outbasename "embedded"
set no_orient 0

for { set i 0 } { $i < [llength $argv] } { incr i } {
   if { [lindex $argv $i] == "-psf"} {
      incr i
      set psf [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-pdb"} {
      incr i
      set pdb [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-bilayer_pdb"} {
      incr i
      set bilayer_pdb [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-bilayer_psf"} {
      incr i
      set bilayer_psf [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-bilayer_xsc"} {
      incr i
      set bilayer_xsc [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-no_orient"} {
      incr i
      set no_orient 1
   }
   if { [lindex $argv $i] == "-z_head_group"} {
      incr i
      set z_head_group [deprotect_str_arg [lindex $argv $i]]
   }
   if { [lindex $argv $i] == "-z_tail_group"} {
      incr i
      set z_tail_group [deprotect_str_arg [lindex $argv $i]]
   }
   if { [lindex $argv $i] == "-z_ref_group"} {
      incr i
      set z_ref_group [deprotect_str_arg [lindex $argv $i]]
   }
   if { [lindex $argv $i] == "-z_value"} {
      incr i
      set z_value [lindex $argv $i]
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
if { $bilayer_pdb == "" } {
   puts "Error: -bilayer_pdb argument required"
   exit
}
if { $bilayer_psf == "" } {
   puts "Error: -bilayer_psf argument required"
   exit
}
if { $bilayer_xsc == "" } {
   puts "Error: -bilayer_xsc argument required"
   exit
}
# check for optional arguments
if { $z_head_group == "" && !$no_orient } {
   puts "Error: -z_head_group argument required if -no_orient is not indicated"
   exit
}
if { $z_tail_group == "" && !$no_orient } {
   puts "Error: -z_tail_group argument required if -no_orient is not indicated"
   exit
}

mol new $psf
mol addfile $pdb waitfor all
set protein [molinfo top get id]
vmdcon -info "protein system from $psf and $pdb into molecule $protein"

mol new $bilayer_psf
mol addfile $bilayer_pdb waitfor all
set bilayer [molinfo top get id]
vmdcon -info "bilayer system from $bilayer_psf and $bilayer_pdb into molecule $bilayer"
set bilayer_sel [atomselect $bilayer "all"]

vmdcon -info "selecting \"$z_head_group\" and \"$z_tail_group\" from protein (molecule $protein)"
set head [atomselect $protein "$z_head_group"]
set tail [atomselect $protein "$z_tail_group"]
vmdcon -info "head: [$head num] atoms"
vmdcon -info "tail: [$tail num] atoms"
set head_com [measure center $head weight mass]
set tail_com [measure center $tail weight mass]
set bilayer_com [measure center $bilayer_sel weight mass]
set bilayer_x [lindex $bilayer_com 0]
set bilayer_y [lindex $bilayer_com 1]
set bilayer_z [lindex $bilayer_com 2]
set pro_sel [atomselect $protein "all"]
if { !$no_orient } {
   vmdcon -info "orienting protein axis to bilayer normal"
   set pro_axis [vecnorm [vecsub $head_com $tail_com]]
   set A [orient $pro_sel $pro_axis {0 0 1}]
   $pro_sel move $A
}
set pro_mid_z_ref_sel [atomselect $protein "$z_ref_group"]
set pro_embed_mid_z [lindex [measure center $pro_mid_z_ref_sel weight mass] 2]
set pro_com [measure center $pro_sel weight mass]
set pro_x [lindex $pro_com 0]
set pro_y [lindex $pro_com 1]
set pro_x_shift [expr $bilayer_x - $pro_x]
set pro_y_shift [expr $bilayer_y - $pro_y]
set pro_z_shift [expr $bilayer_z - $pro_embed_mid_z - $z_value]

$pro_sel moveby [list $pro_x_shift $pro_y_shift $pro_z_shift]
$pro_sel writepdb "${outbasename}_embedded.pdb"

set bad_atoms [measure contacts 2.4 $bilayer_sel $pro_sel]
set bad_membrane_idx [lindex $bad_atoms 0]
set bad_membrane_sel [atomselect $bilayer "index $bad_membrane_idx"]

readpsf $psf pdb ${outbasename}_embedded.pdb
readpsf $bilayer_psf pdb $bilayer_pdb

foreach seg [$bad_membrane_sel get segname] resid [$bad_membrane_sel get resid] {
   delatom $seg $resid
}
regenerate angles dihedrals

writepsf cmap ${outbasename}.psf
writepdb ${outbasename}.pdb

resetpsf
mol delete $protein
mol delete $bilayer
mol new ${outbasename}.psf
mol addfile ${outbasename}.pdb waitfor all
set allatoms [atomselect top "all"]
set net_charge [vecsum [$allatoms get charge]]
vmdcon -info "net charge of embedded system before filling: $net_charge"

set embedded_system [molinfo top get id]
set embedded_protein_sel [atomselect ${embedded_system} "protein or glycan"]
set embedded_protein_minmax [measure minmax $embedded_protein_sel]
set embedded_protein_min [lindex $embedded_protein_minmax 0]
set embedded_protein_max [lindex $embedded_protein_minmax 1]
set embedded_protein_min_z [lindex $embedded_protein_min 2]
set embedded_protein_max_z [lindex $embedded_protein_max 2]
set fp [open "${outbasename}_embed_prefilling.yaml" "w"]
puts $fp "embedded_protein_min_z: $embedded_protein_min_z"
puts $fp "embedded_protein_max_z: $embedded_protein_max_z"
puts $fp "net_charge: $net_charge"
puts $fp "prefill_xsc: $bilayer_xsc"
close $fp
