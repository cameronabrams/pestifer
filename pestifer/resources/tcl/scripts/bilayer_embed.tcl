# Author: Cameron F. Abrams, <cfa22@drexel.edu>

# VMD/psfgen script for creating a new psf/pdb pair for a complete, 
# membrane-embedded protein pdb from an intact membrane system
# and the original protein psf/pdb

# if referenced using the Psfgen scriptwriter, all common psfgen
# pre-build commands are invoked automatically
# and the script is run in the context of the psfgen package

package require Orient
namespace import ::Orient::*
package require pbctools
package require solvate
package require autoionize

set scriptname bilayer_embed

set pdb "";  # psf of protein-only system
set psf "";  # pdb of protein-only system
set bilayer_pdb ""; # the membrane-only system
set bilayer_psf ""; # the membrane-only system
set bilayer_xsc ""; # the membrane-only system
set z_head_group ""; # atomselection whose center of mass defines one end of the protein axis
set z_tail_group ""; # atomselection whose center of mass defines the other end of the protein axis
set z_ref_group ""; # atomselection whose center of mass lies at the membrane midplane + z_value
set z_value 0.0 ; # offset of the protein center of mass from the bilayer midplane
set outbasename "embedded"
set no_orient 0; # override the orientation of the protein axis to align with the bilayer normal
set zdist 10.0; # distance between z-extremal protein atoms and z-boundaries of box
set sc 0.0
set cation POT
set anion CLA

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
      set no_orient_val [lindex $argv $i]
      if { $no_orient_val == "1" || $no_orient_val == "True" || $no_orient_val == "true" || $no_orient_val == "yes" || $no_orient_val == "Yes" } {
         set no_orient 1
      } else {
         set no_orient 0
      }
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
   if { [lindex $argv $i] == "-zdist"} {
      incr i
      set zdist [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-sc"} {
      incr i
      set sc [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-cation"} {
      incr i
      set cation [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-anion"} {
      incr i
      set anion [lindex $argv $i]
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

# load the protein system
mol new $psf
mol addfile $pdb waitfor all
set protein [molinfo top get id]
set pro_sel [atomselect $protein "all"]

# load the bilayer system, take measurements
mol new $bilayer_psf
mol addfile $bilayer_pdb waitfor all
set bilayer [molinfo top get id]
pbc readxst $bilayer_xsc -molid $bilayer
set bilayer_box [lindex [pbc get -molid $bilayer] 0]
set bilayer_box_Lx [lindex $bilayer_box 0]
set bilayer_box_Ly [lindex $bilayer_box 1]
set bilayer_sel [atomselect $bilayer "all"]

set bilayer_minmax [measure minmax $bilayer_sel]
set bilayer_min [lindex $bilayer_minmax 0]
set bilayer_max [lindex $bilayer_minmax 1]
set bilayer_min_x [lindex $bilayer_min 0]
set bilayer_max_x [lindex $bilayer_max 0]
set bilayer_mid_x [expr ($bilayer_min_x + $bilayer_max_x) / 2]
set bilayer_min_y [lindex $bilayer_min 1]
set bilayer_max_y [lindex $bilayer_max 1]
set bilayer_mid_y [expr ($bilayer_min_y + $bilayer_max_y) / 2]
set bilayer_min_z [lindex $bilayer_min 2]
set bilayer_max_z [lindex $bilayer_max 2]
set bilayer_mid_z [expr ($bilayer_min_z + $bilayer_max_z) / 2]

vmdcon -info "bilayer z-span [lindex $bilayer_box 2]; by apparent atom measurements, min-z $bilayer_min_z max-z $bilayer_max_z"

set bilayer_com [measure center $bilayer_sel weight mass]
set bilayer_com_z [lindex $bilayer_com 2]

if { !$no_orient } {
   vmdcon -info "orienting protein axis to bilayer normal"
   set head [atomselect $protein "$z_head_group"]
   set tail [atomselect $protein "$z_tail_group"]
   set head_com [measure center $head weight mass]
   set tail_com [measure center $tail weight mass]
   set pro_axis [vecnorm [vecsub $head_com $tail_com]]
   set A [orient $pro_sel $pro_axis {0 0 1}]
   $pro_sel move $A
} else {
   vmdcon -info "not orienting protein axis to bilayer normal"
}

# perform a translation of the protein to the middle of the bilayer
set pro_mid_z_ref_sel [atomselect $protein "$z_ref_group"]
set pro_embed_mid_z [lindex [measure center $pro_mid_z_ref_sel weight mass] 2]
set pro_com [measure center $pro_sel weight mass]
set pro_x [lindex $pro_com 0]
set pro_y [lindex $pro_com 1]
set pro_x_shift [expr $bilayer_mid_x - $pro_x]
set pro_y_shift [expr $bilayer_mid_y - $pro_y]
set pro_z_shift [expr $bilayer_com_z - $pro_embed_mid_z - $z_value]
$pro_sel moveby [list $pro_x_shift $pro_y_shift $pro_z_shift]
$pro_sel writepdb "${outbasename}_embedded.pdb"

set pro_minmax [measure minmax $pro_sel]
set pro_min [lindex $pro_minmax 0]
set pro_max [lindex $pro_minmax 1]
set pro_min_x [lindex $pro_min 0]
set pro_max_x [lindex $pro_max 0]
set pro_min_y [lindex $pro_min 1]
set pro_max_y [lindex $pro_max 1]
set pro_min_z [lindex $pro_min 2]
set pro_max_z [lindex $pro_max 2]
set box_min_z [expr $pro_min_z - $zdist]
set box_max_z [expr $pro_max_z + $zdist]

# delete atoms that are in conflict with the protein
set bad_atoms [measure contacts 2.4 $bilayer_sel $pro_sel]
set bad_membrane_idx [lindex $bad_atoms 0]
set bad_membrane_sel [atomselect $bilayer "index $bad_membrane_idx"]

readpsf $psf pdb ${outbasename}_embedded.pdb
readpsf $bilayer_psf pdb $bilayer_pdb

catch {foreach seg [$bad_membrane_sel get segname] resid [$bad_membrane_sel get resid] {
   delatom $seg $resid
}} result
regenerate angles dihedrals

writepsf cmap ${outbasename}_prefill.psf
writepdb ${outbasename}_prefill.pdb

# add water slabs to fill gaps above and below the protein
resetpsf
mol delete $protein
mol delete $bilayer
mol new ${outbasename}_prefill.psf
mol addfile ${outbasename}_prefill.pdb waitfor all
set raw_embedded_system [molinfo top get id]
set allatoms [atomselect $raw_embedded_system "all"]
$allatoms writenamdbin ${outbasename}.coor
set segnames [lsort -unique [$allatoms get segname]]
set solsegnums [list]
foreach sg $segnames {
   if { [string match "WT*" $sg] } {
      if {[regexp {[0-9]+} $sg match number]} {
         regsub -all {[^0-9]} $sg "" digits
         lappend solsegnums $digits
      }
   }
}
set solsegnums [lsort $solsegnums]
vmdcon -info "solvent (WT) segment numbers in raw-embedded bilayer: $solsegnums"
if { [llength $solsegnums] == 0 } {
   set nextsolsegnum 1
} else {
   set nextsolsegnum [expr {[lindex $solsegnums end] + 1}]
}
vmdcon -info "next solvent segname: WT$nextsolsegnum"

set net_charge [vecsum [$allatoms get charge]]
mol delete $raw_embedded_system

set addl_water [list]

if { $box_min_z < $bilayer_min_z } {
   # make a slab of water thick enough to fill this gap
   set gapsize [expr $bilayer_min_z - $box_min_z]
   if {$gapsize < 3.0} {
      set extragap [expr 3.0 - $gapsize]
      set $box_min_z [expr $box_min_z - $extragap]
   }
   vmdcon -info "solvating into {{0 0 $box_min_z} {$bilayer_box_Lx $bilayer_box_Ly $bilayer_min_z}}"
   vmdcon -info "Running solvate [list [list 0 0 $box_min_z] [list $bilayer_box_Lx $bilayer_box_Ly $bilayer_min_z]] -o ${outbasename}_water_lower"
   solvate -minmax [list [list 0 0 $box_min_z] [list $bilayer_box_Lx $bilayer_box_Ly $bilayer_min_z]] -o ${outbasename}_water_lower
   if { $sc > 0.0 } {
      vmdcon -info "Adding $sc M salt (cation $cation, anion $anion) to lower water slab"
      autoionize -psf ${outbasename}_water_lower.psf -pdb ${outbasename}_water_lower.pdb -o ${outbasename}_solution_lower -sc $sc -cation $cation -anion $anion
   } {
      file copy ${outbasename}_water_lower.psf ${outbasename}_solution_lower.psf
      file copy ${outbasename}_water_lower.pdb ${outbasename}_solution_lower.pdb
   }
   lappend addl_solution ${outbasename}_solution_lower
} else {
   set box_min_z $bilayer_min_z
}

if { $box_max_z > $bilayer_max_z } {
   # make a slab of water thick enough to fill this gap
   set gapsize [expr $box_max_z - $bilayer_max_z]
   if {$gapsize < 3.0} {
      set extragap [expr 3.0 - $gapsize]
      set $box_max_z [expr $box_max_z + $extragap]
   }
   vmdcon -info "Running solvate -minmax [list [list 0 0 $bilayer_max_z] [list $bilayer_box_Lx $bilayer_box_Ly $box_max_z]] -o ${outbasename}_water_upper"
   solvate -minmax [list [list 0 0 $bilayer_max_z] [list $bilayer_box_Lx $bilayer_box_Ly $box_max_z]] -o ${outbasename}_water_upper
   if { $sc > 0.0 } {
      vmdcon -info "Adding $sc M salt (cation $cation, anion $anion) to upper water slab"
      autoionize -psf ${outbasename}_water_upper.psf -pdb ${outbasename}_water_upper.pdb -o ${outbasename}_solution_upper -sc $sc -cation $cation -anion $anion
   } else {
      file copy ${outbasename}_water_upper.psf ${outbasename}_solution_upper.psf
      file copy ${outbasename}_water_upper.pdb ${outbasename}_solution_upper.pdb
   }
   lappend addl_solution ${outbasename}_solution_upper
} else {
   set box_max_z $bilayer_max_z
}

set newsegids [list]
resetpsf
readpsf ${outbasename}_prefill.psf pdb ${outbasename}_prefill.pdb
foreach aw $addl_solution {
   set apdb ${aw}.pdb
   set apsf ${aw}.psf
   mol new $apsf
   mol addfile $apdb waitfor all
   set solution [molinfo top get id]
   set water_sel [atomselect $solution "water"]
   set new_segids [lsort -unique [$water_sel get segname]]
   foreach ns $new_segids {
      set pdbfrag [atomselect $solution "segname $ns"]
      $pdbfrag writepdb ${aw}_${ns}.pdb
      segment WT${nextsolsegnum} {
         pdb ${aw}_${ns}.pdb
         first none
         last none
         auto none
      }
      coordpdb ${aw}_${ns}.pdb WT${nextsolsegnum}
      lappend newsegids WT${nextsolsegnum}
      set nextsolsegnum [expr $nextsolsegnum + 1]
   }
   set ions [atomselect $solution "ion"]
   set ion_segs [lsort -unique [$ions get segname]]
   foreach is $ion_segs {
      set ion_frag [atomselect $solution "segname $is"]
      $ion_frag writepdb ${aw}_${is}.pdb
      segment $is {
         pdb ${aw}_${is}.pdb
         first none
         last none
         auto none
      }
      coordpdb ${aw}_${is}.pdb $is
      lappend newsegids $is
   }
   mol delete $solution
   file delete $apdb
   file delete $apsf
   foreach ns $new_segids {
      file delete ${aw}_${ns}.pdb
   }
}

writepsf cmap ${outbasename}_solvent_appended.psf
writepdb ${outbasename}_solvent_appended.pdb

# delete any waters that conflict with the protein after the slab addition
mol new ${outbasename}_solvent_appended.psf
mol addfile ${outbasename}_solvent_appended.pdb waitfor all
set embedded_system [molinfo top get id]
set all [atomselect $embedded_system "all"]
$all moveby [list 0 0 expr (-1*$box_min_z)]
set box_Lz [expr $box_max_z - $box_min_z]
set pro_sel [atomselect $embedded_system "protein or glycan"]
set segs_to_search [join $newsegids]
set water_sel [atomselect $embedded_system "segname $segs_to_search"]
set bad_atoms [measure contacts 2.4 $water_sel $pro_sel]
set bad_water_idx [lindex $bad_atoms 0]
set bad_water_sel [atomselect $embedded_system "index $bad_water_idx"]
catch {foreach seg [$bad_water_sel get segname] resid [$bad_water_sel get resid] {
   delatom $seg $resid
}} result
regenerate angles dihedrals
writepsf cmap ${outbasename}_filled.psf
writepdb ${outbasename}_filled.pdb

# write the resulting box dimensions to an xsc file
set fp [open "${outbasename}.xsc" "w"]
puts $fp "# PESTIFER generated xst file"
puts $fp "#\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w"
puts $fp "0 $bilayer_box_Lx 0 0  0 $bilayer_box_Ly 0  0 0 $box_Lz  [expr ($bilayer_box_Lx / 2)] [expr ($bilayer_box_Ly / 2)] [expr ($box_Lz / 2)] 0 0 0 0 0 0"
close $fp

resetpsf
mol delete $embedded_system

# neutralize the system if necessary

vmdcon -info "net charge of embedded system before filling: $net_charge"

if {[expr abs($net_charge)] > 0.0001} {
   autoionize -psf ${outbasename}_filled.psf -pdb ${outbasename}_filled.pdb -o ${outbasename} -neutralize
} else {
   vmdcon -info "no autoionization required; net charge is $net_charge"
   file copy ${outbasename}_filled.psf ${outbasename}.psf
   file copy ${outbasename}_filled.pdb ${outbasename}.pdb
}
mol new ${outbasename}.psf
mol addfile ${outbasename}.pdb waitfor all
[atomselect top "all"] writenamdbin ${outbasename}.coor

vmdcon -info "Final embedded system: ${outbasename}.psf ${outbasename}.pdb ${outbasename}.coor ${outbasename}.xsc"
