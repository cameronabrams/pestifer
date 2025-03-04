# Author: Cameron F. Abrams, <cfa22@drexel.edu>

# VMD script for generating solvated/neutralized system from the protein psf/pdb

set scriptname solv

set pad 10; # default pad in angstroms
set pdb "empty.pdb"
set psf "empty.psf"
set outpre "ionized"
set cubic 0
set cellfile "cell.inp"

for { set i 0 } { $i < [llength $argv] } { incr i } {
    if { [lindex $argv $i] == "-pad" } {
       incr i
       set pad [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-cubic" } {
       set cubic 1
    }
    if { [lindex $argv $i] == "-pdb"} {
       incr i
       set pdb [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-psf"} {
       incr i
       set psf [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-o"} {
       incr i
       set outbasename [lindex $argv $i]
    }
}

if { ! [file exists $psf] } {
   vmdcon -err "${scriptname}: $psf not found."
   exit
}
if { ! [file exists $pdb] } {
   vmdcon -err "${scriptname}: $pdb not found."
   exit
}
mol new $psf
mol addfile $pdb

set a [atomselect top all]

set box { { ? ? ? } { ? ? ? } }
set basisvec { ? ? ? }
set origin { ? ? ? }

set minmax [measure minmax $a]
set maxspan [expr -9999]
foreach d {0 1 2} {
   set thisspan [expr [lindex $minmax 1 $d ] - [lindex $minmax 0 $d]]
   if { $thisspan > $maxspan } {
      set maxspan $thisspan
   }
}
vmdcon -info "${scriptname}: Maximum span [format %.2f $maxspan] Angstroms"

set sympad [list 0 0 0]
if { $cubic == 1 } {
   vmdcon -info "${scriptname}: Enforcing cubic box"
   foreach d {0 1 2} {
      set thisspan [expr [lindex $minmax 1 $d ] - [lindex $minmax 0 $d]]
      lset sympad $d [expr 0.5*($maxspan-$thisspan)]
   }
}

# a note on origin from NAMD docs:
# cellOrigin $ <$ center of periodic cell (Å) $ >$
# Acceptable Values: position
# Default Value: 0 0 0
# Description: When position rescaling is used to control pressure, this location will remain constant. Also used as the center of the cell for wrapped output coordinates.

foreach d {0 1 2} {
  lset box 0 $d [expr [lindex $minmax 0 $d] - $pad - [lindex $sympad $d]]
  lset box 1 $d [expr [lindex $minmax 1 $d] + $pad + [lindex $sympad $d]]
  lset basisvec $d [expr [lindex $box 1 $d ] - [lindex $box 0 $d]] 
  lset origin $d [expr 0.5*([lindex $box 1 $d ] + [lindex $box 0 $d])] 
}

package require solvate
package require autoionize
psfcontext mixedcase

solvate $psf $pdb -minmax $box -o ${outbasename}_water
autoionize -psf ${outbasename}_water.psf -pdb ${outbasename}_water.pdb -neutralize -o ${outbasename}

# generate an input file for the first solvated MD simulation
# namd config file
set fp [open "${outbasename}_cell.tcl" "w"]
puts $fp "cellbasisvector1 [lindex $basisvec 0] 0 0"
puts $fp "cellbasisvector2 0 [lindex $basisvec 1] 0"
puts $fp "cellbasisvector3 0 0 [lindex $basisvec 2]"
puts $fp "cellorigin $origin"
close $fp
vmdcon -info "${scriptname}: Generated ${outbasename}_cell.tcl."

set fp [open "${outbasename}.xsc" "w"]
puts $fp "#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z"
puts $fp "0 [lindex $basisvec 0] 0 0 0 [lindex $basisvec 1] 0 0 0 [lindex $basisvec 2] $origin"
close $fp
vmdcon -info "${scriptname}: Generated ${outbasename}.xsc."


quit
