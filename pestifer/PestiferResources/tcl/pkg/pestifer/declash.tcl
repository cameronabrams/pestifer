# Author: Cameron F. Abrams, <cfa22@drexel.edu>

package provide PestiferDeclash 1.0

namespace eval ::PestiferDeclash:: {
    namespace export *
}
package require PestiferCRot
namespace import PestiferCRot::*
#
# A very simple Metropolis algorithm to alter conformation of raw model-built
# protein gap loops to minimize clashes
proc PestiferDeclash::declash_loop { molid segname loop maxcycles } {
  set nr [llength $loop]
  set loopsel [atomselect $molid "segname $segname and resid [join $loop]"]
  set residue_numbers [[atomselect $molid "[$loopsel text] and name CA"] get residue]
  set env [atomselect $molid "same residue as exwithin 4.0 of (segname $segname and resid [join $loop])"]
  set residuenum_end [lindex $residue_numbers end]
  vmdcon -info "DECLASH_LOOP) molid $molid segname $segname loop $loop maxcycles $maxcycles"
  vmdcon -info "DECLASH_LOOP) loopsel [$loopsel num] atoms; residue_numbers $residue_numbers"
  vmdcon -info "DECLASH_LOOP) env [$env num] atoms"
  for { set i 0 } { $i < $nr } { incr i } {
    # rotate phi angle and psi angle to minimize number of contacts between residue and 
    # its environment
    set rsel [atomselect $molid "segname $segname and resid [lindex $loop $i] to [lindex $loop end]"]
    set residuenum1 [lindex $residue_numbers $i]
    set CON_STRUCT [measure contacts 2.0 $rsel $env]
    set CON [llength [lindex $CON_STRUCT 0]]
    vmdcon -info "DECLASH_LOOP) [lindex $loop $i] INIT $CON"
    for { set t 0 } { $t < $maxcycles && $CON > 0 } { incr t } {
      set SAVEPOS [$loopsel get {x y z}]
      set rphi [expr (1-2*rand())*120.0]
      #set rpsi [expr (1-2*rand())*60.0]
      Crot_phi_toCterm $residuenum1 $residuenum_end $segname $molid $rphi
      if { $i > [expr $nr - 1] } {
        Crot_psi_toCterm $residuenum1 $residuenum_end $segname $molid $rphi
      }
      $env update
      set TRICON_STRUCT [measure contacts 2.0 $rsel $env]
      set TRICON  [llength [lindex $TRICON_STRUCT 0]]
      if { [expr $TRICON < $CON] } {
        # accept this move
        set CON $TRICON
        vmdcon -info "DECLASH_LOOP) ${segname}:[lindex $loop 0]-[lindex $loop $i] $t $CON"
      } else {
        # reject this move
        $loopsel set {x y z} $SAVEPOS
      }
    }
    $rsel delete
  }
  $env delete
  $loopsel delete
}

# General MC for reducing steric clashes between a pendant group and the
# rest of a molecule.  
#
# Parameters
# ----------
# molid: the molecule id
# indices: list of atom indices defining the pendant group
# bonds: list of rotatable bonds in the pendant group
# movers: list of lists of movable atoms per rotatable bond
# maxcycles: maximum number of MC cycles (one bond rotation per)
# clashdist: distance threshold for designating whether two atoms "clash" (A)
#
# indices, bonds, and movers must be determined outside this procedure
# The MC loop terminates either when the maximum number of iterations is
# reacted or the number of interatomic clashes drops to zero.
#
# If used correctly, this procedure can supersede the declash_loop procedure
#
proc PestiferDeclash::declash_pendant { molid indices bonds movers maxcycles clashdist } {
   set degs {-120 -60 60 120 180}
   set atomsel [atomselect $molid "noh and index $indices"]
   set environ [atomselect $molid "noh"]
   set nbonds [llength $bonds]
   vmdcon -info "Declash environ has [$environ num] atoms"
   set ncontacts [llength [lindex [measure contacts $clashdist $atomsel $environ] 0]]
   vmdcon -info "Declash pendant: Pendant has $ncontacts initial atomic clashes"
   if {$ncontacts > 0} {
      vmdcon -info "Declashing via maximally $maxcycles cycles"
      for {set i 0} {$i < $maxcycles} {incr i} {
         set ridx [expr {int(rand()*$nbonds)}]
         set bo [lindex $bonds $ridx]
         set b [[atomselect $molid "index $bo"] get {x y z}]
         set this_mover [atomselect $molid "index [lindex $movers $ridx]"]
         set didx [expr {int(rand()*[llength $degs])}]
         set deg [lindex $degs $didx]
         set tmat [trans bond [lindex $b 0] [lindex $b 1] $deg degrees]
         $this_mover move $tmat
         set newcontacts [llength [lindex [measure contacts $clashdist $atomsel $environ] 0]]
         if { $newcontacts >= $ncontacts } {
            set deg [expr -1 * ($deg)]
            set tmat [trans bond [lindex $b 0] [lindex $b 1] $deg degrees]
            $this_mover move $tmat
         } else {
            set ncontacts $newcontacts
         }
         vmdcon -info "  cycle $i bond $bo deg $deg ncontacts $ncontacts"
         $this_mover delete
         if { $ncontacts == 0 } {
            break
         }
      }
   }
   $environ delete
   $atomsel delete
}
