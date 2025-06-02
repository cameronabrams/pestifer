# Author: Cameron F. Abrams, <cfa22@drexel.edu>

# VMD/psfgen script for creating a new psf/pdb pair for a bilayer patch
# from a packmol pdb output

# if referenced using the Psfgen scriptwriter, all common psfgen
# pre-build commands are invoked automatically
# and the script is run in the context of the psfgen package

package require PestiferEnviron 1.0
namespace import ::PestiferEnviron::*

set scriptname bilayer_patch

set pdb "";  # output of packmol
set outbasename "bilayer_patch-parameterized"

for { set i 0 } { $i < [llength $argv] } { incr i } {
    if { [lindex $argv $i] == "-pdb"} {
       incr i
       set pdb [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-o"} {
       incr i
       set outbasename [lindex $argv $i]
    }
}

mol new $pdb waitfor all
set environ_molid [molinfo top get id]

set waters [atomselect $environ_molid "water"]
set ions [atomselect $environ_molid "ion"]
set lipids [atomselect $environ_molid "lipid"]
set other [atomselect $environ_molid "not water and not ion and not lipid"]

vmdcon -info "lipids: [$lipids num] atoms"
vmdcon -info "waters: [$waters num] atoms"
vmdcon -info "ions: [$ions num] atoms"
vmdcon -info "other: [$other num] atoms"

set segtypes {lipid ion water}
set seglabels {L I WT}
set segidx {1 1 1}
set maxres_per_seg 9999

write_psfgen $environ_molid $segtypes $seglabels $segidx $maxres_per_seg

mol delete $environ_molid
regenerate angles dihedrals

# pestifer takes over to write the psf and pdb