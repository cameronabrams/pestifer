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

foreach segtype [list lipid ion water] {

   set a [atomselect $add_molid "$segtype"]
   set input_chains [lsort -unique [$a get chain]]
   foreach ic $input_chains {
      set new_chain($ic) $next_available_chain
      set next_available_chain [letter_up $next_available_chain]
   }

   foreach chain $input_chains {
      set tsel [atomselect $add_molid "chain $chain and $segtype"]
      set nc $new_chain($chain)
      $tsel set chain $nc
      $tsel writepdb "${nc}_tmp.pdb"
      segment $nc {
         pdb ${nc}_tmp.pdb
      }
      coordpdb ${nc}_tmp.pdb $nc
   }
}