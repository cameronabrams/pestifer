# Author: Cameron F. Abrams, <cfa22@drexel.edu>

# VMD script for creating a new psf/pdb pair for a complete, 
# membrane-embedded protein pdb and the original protein psf/pdb

set scriptname memb

set pdb ""
set psf ""
set addpdb ""
set outbasename "memb-parameterized"

for { set i 0 } { $i < [llength $argv] } { incr i } {
    if { [lindex $argv $i] == "-psf"} {
       incr i
       set psf [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-pdb"} {
       incr i
       set pdb [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-addpdb"} {
       incr i
       set addpdb [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-o"} {
       incr i
       set outbasename [lindex $argv $i]
    }
}

if { ! [file exists $addpdb] } {
   vmdcon -err "${scriptname}: $addpdb not found."
   exit
}

# if no pdb and/or no psf, we must be processing a membrane-only addpdb
set init_model ""
set init_chains [list]
set next_available_chain A
if { $pdb != "" && $psf != "" } {
   readpsf $psf pdb $pdb
   mol new $psf
   mol addfile $pdb waitfor all
   set init_model [atomselect top all]
   set init_chains [lsort -unique [$init_model get chain]]
   set next_available_chain [letter_up [lindex $init_chains end]]
   set init_molid [molinfo top get id]
}

mol new $addpdb waitfor all
set add_molid [molinfo top get id]
set waters [atomselect $add_molid "water"]
set ions [atomselect $add_molid "ion"]
set lipids [atomselect $add_molid "lipid"]
set apparent_insert [atomselect $add_molid "not water and not ion and not lipid"]

if { $init_model != "" } {
   vmdcon -info "init_model: [$init_model num] atoms"
   vmdcon -info "apparent insert: [$apparent_insert num] atoms"
   vmdcon -info "added lipids: [$lipids num] atoms"
   vmdcon -info "added waters: [$waters num] atoms"
   vmdcon -info "added ions: [$ions num] atoms"
}

set maxr_per_seg 1000

foreach segtype [list lipid ion water] S [list L I W] {
   set a [atomselect $add_molid "$segtype"]
   $a set chain $next_available_chain

   set ridx [lsort -integer -unique [$a get residue]]
   set nres [llength $ridx]
   set nseg [expr $nres / $maxr_per_seg + 1]
   set mm [expr $nres % $nseg]
   vmdcon -info "[$a num] $segtype atoms (chain $next_available_chain) in $nres residues divide into $nseg segments"

   for { set seg 1 } { $seg <= $nseg } { incr seg } {
      set segname "${next_available_chain}${seg}"
      set left [expr (${seg}-1)*$maxr_per_seg]
      set right [expr $left + $maxr_per_seg - 1]
      if { $right >= [llength $ridx] } {
         set right end
      }
      vmdcon -info "Segment $segname is residues [lindex $ridx $left] to [lindex $ridx $right]"
      set res [lrange $ridx $left $right]
      set segsel [atomselect top "$segtype and residue $res"]
      set sridx [lsort -integer -unique [$segsel get residue]]
      set rser 1
      foreach x $sridx {
         set sridx_sermap($x) $rser
         incr rser
      }
      set rser [list]
      foreach x [$segsel get residue] {
         lappend rser $sridx_sermap($x)
      }
      $segsel set resid $rser
      vmdcon -info "Segment $segname has [$segsel num] atoms"
      $segsel writepdb "${segname}_tmp.pdb"
      segment $segname {
         auto none
         first none
         last none
         pdb ${segname}_tmp.pdb
      }
      coordpdb ${segname}_tmp.pdb $segname
   }
   set next_available_chain [letter_up $next_available_chain]
}
