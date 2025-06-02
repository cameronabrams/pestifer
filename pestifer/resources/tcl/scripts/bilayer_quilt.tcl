# Author: Cameron F. Abrams, <cfa22@drexel.edu>

# VMD/psfgen script for creating a new psf/pdb pair for a
# bilayer created using the upper leaflet of one
# input bilayer and the lower leaflet of another, and
# replicated in x and y to make a "quilt"

package require PestiferEnviron 1.0
namespace import ::PestiferEnviron::*

package require pbctools

set scriptname bilayer_quilt

set propdb ""; # pdb of protein that will eventually be embedded
set margin 10.0; # Angstroms; min distance from any protein atom to box edge in x and y
# A = system defining upper leaflet
set psfA ""
set pdbA ""
set xscA ""
# B = system defining lower leaflet
set psfB ""
set pdbB ""
set xscB ""

set outbasename "quilt"
set npatchx 1
set npatchy 1
set dimx 0
set dimy 0
for { set i 0 } { $i < [llength $argv] } { incr i } {
   if { [lindex $argv $i] == "-propdb"} {
      incr i
      set propdb [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-margin"} {
      incr i
      set margin [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-psfA"} {
      incr i
      set psfA [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-psfB"} {
      incr i
      set psfB [lindex $argv $i]
   }
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
   if { [lindex $argv $i] == "-nx"} {
      incr i
      set npatchx [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-ny"} {
      incr i
      set npatchy [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-dimx"} {
      incr i
      set dimx [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-dimy"} {
      incr i
      set dimy [lindex $argv $i]
   }
   if { [lindex $argv $i] == "-o"} {
      incr i
      set outbasename [lindex $argv $i]
   }
}

# determine the lateral dimensions of the system
set min_sys_Lx 0
set min_sys_Ly 0
if { $propdb != "" } {
   # there is a protein pdb that should determine the box size
   mol new $propdb waitfor all
   set prop_molid [molinfo top get id]
   set pro [atomselect $prop_molid all]
   set mm [measure minmax $pro]
   set minx [lindex $mm 0 0]
   set miny [lindex $mm 0 1]
   set maxx [lindex $mm 1 0]
   set maxy [lindex $mm 1 1]
   set pro_Lx [expr $maxx - $minx]
   set pro_Ly [expr $maxy - $miny]
   set min_sys_Lx [expr $pro_Lx + 2*$margin]
   set min_sys_Ly [expr $pro_Ly + 2*$margin]
   vmdcon -info "minimum system (Lx,Ly): $min_sys_Lx, $min_sys_Ly"
   mol delete $prop_molid
} elseif { $dimx > 0 && $dimy > 0 } {
   # explicit dimensions provided
   set min_sys_Lx $dimx
   set min_sys_Ly $dimy
} elseif { $npatchx > 0 && $npatchy > 0 } {
   # explicit dimensions provided in number of patchlengths
   vmdcon -info "x and y dimensions provided in number of patchlengths; x $npatchx, y $npatchy"
} else {
   vmdcon -err "No protein pdb + margin or patch dimensions provided"
   exit
}

# read in the two source bilayer systems, measure lateral areas and determine the common patch area
mol new $psfA
mol addfile $pdbA waitfor all
set source_upper [molinfo top get id]
pbc readxst $xscA -molid $source_upper
set boxA [molinfo $source_upper get {a b c}]
mol new $psfB
mol addfile $pdbB waitfor all
set source_lower [molinfo top get id]
pbc readxst $xscB -molid $source_lower
set boxB [molinfo $source_lower get {a b c}]

set patch_areaA [expr {[lindex $boxA 0] * [lindex $boxA 1]}]
set patch_areaB [expr {[lindex $boxB 0] * [lindex $boxB 1]}]
set patchx [expr max([lindex $boxA 0], [lindex $boxB 0])]
set patchy [expr max([lindex $boxA 1], [lindex $boxB 1])]

vmdcon -info "patch areaA: $patch_areaA"
vmdcon -info "patch areaB: $patch_areaB"
vmdcon -info "patchx: $patchx"
vmdcon -info "patchy: $patchy"

if {$min_sys_Lx > 0 && $min_sys_Ly > 0} {
   # this will be true if a protein is used to size the box or if explicit dimensions
   # are given
   set npatchx [expr max(1, int($min_sys_Lx / $patchx)+1)]
   set npatchy [expr max(1, int($min_sys_Ly / $patchy)+1)]
   vmdcon -info "calculated npatchx x npatchy: $npatchx x $npatchy"
}

# determine the implied quilt lateral dimensions for each pure patch
set quilt_areaA [expr $patch_areaA * $npatchx * $npatchy]
set quilt_areaB [expr $patch_areaB * $npatchx * $npatchy]
vmdcon -info "quilt areaA: $quilt_areaA"
vmdcon -info "quilt areaB: $quilt_areaB"

# extract upper leaflet of system A and lower leaflet of system B
catch {set sliced_upper [leaflet_apportionment $source_upper]}
split_psf $psfA $pdbA $sliced_upper "sliceA"
catch {set sliced_lower [leaflet_apportionment $source_lower]}
split_psf $psfB $pdbB $sliced_lower "sliceB"

set upper_patch "sliceA1"
set lower_patch "sliceB2"

mol delete $source_upper
mol delete $source_lower

# count lipid moleculs in each leaflet
mol new $upper_patch.psf 
mol addfile $upper_patch.pdb waitfor all
set sliceAmolid [molinfo top get id]
set nlipidsA [llength [lsort -unique [[atomselect $sliceAmolid "lipid"] get residue]]]
mol delete $sliceAmolid
mol new $lower_patch.psf
mol addfile $lower_patch.pdb waitfor all
set sliceBmolid [molinfo top get id]
set nlipidsB [llength [lsort -unique [[atomselect $sliceBmolid "lipid"] get residue]]]
mol delete $sliceBmolid

vmdcon -info "upper-patch nlipids: $nlipidsA"
vmdcon -info "lower-patch nlipids: $nlipidsB"
# there should be no molecules loaded at this point
if { [verify_no_mols] != 1} {
   vmdcon -err "There are still molecules loaded"
   exit
}

# reconcile the patch area and number of lipids
set saplA [expr $patch_areaA / $nlipidsA]
set saplB [expr $patch_areaB / $nlipidsB]
vmdcon -info "upper-patch SAPL: $saplA"
vmdcon -info "lower-patch SAPL: $saplB"

# determine the difference in area between the two leaflets
# and determine how many lipids to delete from the larger leaflet
set quilt_area_AB [expr $quilt_areaA - $quilt_areaB]
vmdcon -info "Quilt area difference (upper - lower): $quilt_area_AB"
set nlipids_deleteA 0
set nlipids_deleteB 0
if { $quilt_area_AB > 0 } {
   vmdcon -info "Upper leaflet quilt is larger by $quilt_area_AB"
   set nlipids_deleteA [expr int($quilt_area_AB / $saplA)]
   vmdcon -info "Deleting $nlipids_deleteA lipids from upper leaflet quilt"
   set Lx [lindex $boxA 0]
   set Ly [lindex $boxA 1]
} elseif { $quilt_area_AB < 0 } {
   vmdcon -info "Lower leaflet quilt is larger by $quilt_area_AB"
   set nlipids_deleteB [expr int(-$quilt_area_AB / $saplB)]
   vmdcon -info "Deleting $nlipids_deleteB lipids from lower leaflet quilt"
   set Lx [lindex $boxB 0]
   set Ly [lindex $boxB 1]
} else {
   vmdcon -info "Upper and lower leaflet quilts are equal in area"
   set Lx [lindex $boxA 0]
   set Ly [lindex $boxA 1]
}

# build the quilt by replication of the the hybrid patch
mol new "${upper_patch}.psf" 
mol addfile "${upper_patch}.pdb" waitfor all
set umolid [molinfo top get id]
set upper_patch_sel [atomselect $umolid "all"]
mol new "${lower_patch}.psf" 
mol addfile "${lower_patch}.pdb" waitfor all
set lmolid [molinfo top get id]
set lower_patch_sel [atomselect $lmolid "all"]

# set next_available_chain A
set segtypes {lipid ion water}
set seglabels {L I WT}
set segidx {1 1 1}
set maxres_per_seg 9999
for {set nx 0} { $nx < $npatchx } { incr nx } {
   set xoffset [expr $nx * $Lx]
   for {set ny 0} { $ny < $npatchy } { incr ny } {
      set yoffset [expr $ny * $Ly]
      vmdcon -info "patch ($nx,$ny) offset ($xoffset,$yoffset)"
      set movevec [list $xoffset $yoffset 0]
      set unmovevec [vecscale -1 $movevec]
      vmdcon -info "moving upper patch by $movevec"
      $upper_patch_sel moveby $movevec
      set segidx [write_psfgen $umolid $segtypes $seglabels $segidx $maxres_per_seg]
      vmdcon -info "segidx: $segidx"
      vmdcon -info "moving upper patch by $unmovevec"
      $upper_patch_sel moveby $unmovevec
      vmdcon -info "moving lower patch by $movevec"
      $lower_patch_sel moveby $movevec
      set segidx [write_psfgen $lmolid $segtypes $seglabels $segidx $maxres_per_seg]
      vmdcon -info "moving lower patch by $unmovevec"
      $lower_patch_sel moveby $unmovevec
   }
}
mol delete $umolid
mol delete $lmolid
# there should be no molecules loaded at this point
if { [verify_no_mols] != 1} {
   vmdcon -err "There are still molecules loaded"
   exit
}

# write the first quilt psf and pdb files
regenerate angles dihedrals
set firstname  "${outbasename}_first_quilting"
writepsf "${firstname}.psf"
writepdb "${firstname}.pdb"

# determine final box size and set the origin to be the center of the quilt in all three dimensions
set quilt_boxX [expr $Lx * $npatchx]
set quilt_boxY [expr $Ly * $npatchy]
set quilt_boxZ [expr {max([lindex $boxA 2],[lindex $boxB 2])}]

set origin_x [expr $quilt_boxX / 2]
set origin_y [expr $quilt_boxY / 2]
set origin_z [expr $quilt_boxZ / 2]

# write the box dimensions and origin to an xsc file
set fp [open "${outbasename}.xsc" "w"]
puts $fp "# PESTIFER generated xsc"
puts $fp "#\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w"
puts $fp "0 $quilt_boxX 0 0  0 $quilt_boxY 0  0 0 $quilt_boxZ  $origin_x $origin_y $origin_z 0 0 0 0 0 0"
close $fp

# read in the first quilt psf and pdb files
mol new ${firstname}.psf
mol addfile ${firstname}.pdb waitfor all
set molid [molinfo top get id]

# delete any lipids necessary to equalize the two leaflet areas
set slicedquilt [leaflet_apportionment $molid]

if { $nlipids_deleteA > 0 } {
   vmdcon -info "Deleting $nlipids_deleteA lipids from upper leaflet quilt"
   set u_allres [lindex $slicedquilt 0]
   set u_l [atomselect $molid "lipid and residue $u_allres"]
   set u_lres [lsort -unique [$u_l get residue]]
   set nres [llength $u_lres]
   set shuffled_resnums [shuffle_list $u_lres]
   set remove_resnums [lrange $shuffled_resnums 0 [expr {$nlipids_deleteA - 1}]]
   set remove_sel [atomselect $molid "residue $remove_resnums"]
} elseif { $nlipids_deleteB > 0} {
   vmdcon -info "Deleting $nlipids_deleteB lipids from lower leaflet quilt"
   set l_allres [lindex $slicedquilt 1]
   set l_l [atomselect $molid "lipid and residue $l_allres"]
   set l_lres [lsort -unique [$l_l get residue]]
   set nres [llength $l_lres]
   set shuffled_resnums [shuffle_list $l_lres]
   set remove_resnums [lrange $shuffled_resnums 0 [expr {$nlipids_deleteB - 1}]]
   set remove_sel [atomselect $molid "residue $remove_resnums"]
} else {
   vmdcon -info "No lipids to delete"
   set remove_sel [atomselect $molid "name A_FAKE_NAME"]
}
set badseg [$remove_sel get segname]
set badres [$remove_sel get resid]
set badname [$remove_sel get name]

mol delete $molid
if { [verify_no_mols] != 1} {
   vmdcon -err "There are still molecules loaded"
   exit
}
resetpsf
readpsf ${firstname}.psf pdb ${firstname}.pdb
if { [llength $badseg] > 0 } {
   foreach seg $badseg res $badres name $badname {
      delatom $seg $res $name
   }
   regenerate angles dihedrals
}

writepsf "${outbasename}_prelabel.psf"
writepdb "${outbasename}_prelabel.pdb"
resetpsf

# relabel all segments
mol new ${outbasename}_prelabel.psf
mol addfile ${outbasename}_prelabel.pdb waitfor all
set molid [molinfo top get id]

set segidx {1 1 1}
write_psfgen $molid $segtypes $seglabels $segidx $maxres_per_seg
regenerate angles dihedrals
writepsf "${outbasename}.psf"
writepdb "${outbasename}.pdb"
