# Author: Cameron F. Abrams, <cfa22@drexel.edu>

# VMD/psfgen script for creating a new psf/pdb pair for a complete, 
# membrane-embedded protein pdb from packmol and the original protein psf/pdb
# if a protein psf/pdb is not provided, the script will assume that the script will assume that the
# input pdb is a membrane-only system

# if referenced using the Psfgen scriptwriter, all common psfgen
# pre-build commands are invoked automatically
# and the script is run in the context of the psfgen package

package require PestiferEnviron 1.0
namespace import ::PestiferEnviron::*

set scriptname memb

set pdb "";  # psf of protein-only system
set psf "";  # pdb of protein-only system
set addpdb ""; # pdb output of packmol that contains protein and membrane
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
   vmdcon -err "${scriptname}: packmol output $addpdb not found."
   exit
}

# if no pdb and/or no psf, we must be processing a membrane-only addpdb
set init_model ""
set init_chains [list]
set next_available_chain A
set init_molid -1
if { $pdb != "" && $psf != "" } {
   mol new $psf
   mol addfile $pdb waitfor all
   set init_model [atomselect top all]
   set init_chains [lsort -unique [$init_model get chain]]
   set next_available_chain [letter_up [lindex $init_chains end]]
   set init_molid [molinfo top get id]
}

mol new $addpdb waitfor all
set add_molid [molinfo top get id]
if { $init_molid == -1 } {
   set environ_molid $add_molid
} else {
   set insert_init [atomselect $init_molid all]
   set insert_ser [$insert_init get serial]
   set insert_min_ser [lindex $insert_ser 0]
   set insert_max_ser [lindex $insert_ser end]
   set insert_from_environ [atomselect $add_molid "serial $insert_min_ser to $insert_max_ser"]
   # select atoms from the packmol output that are not in the original protein
   set environ [atomselect $add_molid "serial > $insert_max_ser"]
   # set the coordinates of the protein atoms read in from the psf/pdb
   # to the coordinates of the atoms in the packmol output
   # this gets the protein atoms and the membrane atoms in the same
   # coordinate system, since packmol moves the protein atoms
   $insert_init set {x y z} [$insert_from_environ get {x y z}]
   $insert_init writepdb "insert.pdb"
   readpsf $psf pdb insert.pdb
   $environ writepdb "environ.pdb"
   mol new environ.pdb waitfor all
   set environ_molid [molinfo top get id]
}

set waters [atomselect $environ_molid "water"]
set ions [atomselect $environ_molid "ion"]
set lipids [atomselect $environ_molid "lipid"]
set other [atomselect $environ_molid "not water and not ion and not lipid"]

if { $init_model != "" } {
   vmdcon -info "insert: [$insert_init num] atoms"
}
vmdcon -info "lipids: [$lipids num] atoms"
vmdcon -info "waters: [$waters num] atoms"
vmdcon -info "ions: [$ions num] atoms"
vmdcon -info "other: [$other num] atoms"

set maxr_per_seg 1000

write_psfgen $environ_molid $next_available_chain {lipid water ion} {L I W} $maxr_per_seg

regenerate angles dihedrals

# pestifer takes over to write the psf and pdb