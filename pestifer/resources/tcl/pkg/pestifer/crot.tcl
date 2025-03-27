# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# TcL procedures for facilitating rotations around bonds
package provide PestiferCRot 1.0

namespace eval ::PestiferCRot:: {
    namespace export *
}


# Determine whether or not Sresidues q and r are in the same chain
# Parameters:
# -----------
# - q, r: absolute residue numbers
# - molid: molecule id
#
# Returns:
# --------
# - 0 if q and r are not in the same chain, or if either q or r is less than zero
# - 1 if q and r are in the same chain
proc PestiferCRot::residues_in_same_chain { molid q r } {
   if { $q < 0 || $r < 0 } {
      return 0
   }
   set Q [atomselect $molid "residue $q"]
   set R [atomselect $molid "residue $r"]
   set result 1
   if { [$Q num]==0 || [$R num]==0 } {
      set result 0
   } else {
      set CQ [lindex [$Q get chain] 0]
      set CR [lindex [$R get chain] 0]
      if { $CQ != $CR } {
         set result 0
      }
   }
   $Q delete
   $R delete
   return $result
}

# Measure phi, psi, and omega at residue r
# Parameters:
# -----------
# - r: absolute residue number
# - molid: molecule id
#
# Returns:
# --------
# - list of phi, psi, and omega values, in that order
proc PestiferCRot::get_phi_psi_omega { r molid } {
   set protein_check [expr [[atomselect $molid "residue $r and protein"] num]]
   if { $protein_check == 0 } {
      vmdcon -info "Residue $r of mol $molid is not protein"
      return [list "NaN" "NaN" "NaN"]
   }
   set r_prev [expr $r-1]
   set r_next [expr $r+1]
   set phi "NaN"
   set psi "NaN"
   set omega "NaN"
   if { [residues_in_same_chain $molid $r_prev $r] > 0} {
      set nnc [[atomselect $molid "residue $r_prev and name C"] get index]
      set n [[atomselect $molid "residue $r and name N"] get index]
      set ca [[atomselect $molid "residue $r and name CA"] get index]
      set c  [[atomselect $molid "residue $r and name C"] get index]
      set phi [measure dihed [list $nnc $n $ca $c]]
   }
   if { [residues_in_same_chain $molid $r $r_next] > 0 } {
      set n [[atomselect $molid "residue $r and name N"] get index]
      set ca [[atomselect $molid "residue $r and name CA"] get index]
      set c  [[atomselect $molid "residue $r and name C"] get index]
      set ccn [[atomselect $molid "residue $r_next and name N"] get index]
      set ccca [[atomselect $molid "residue $r_next and name CA"] get index]
      set psi [measure dihed [list $n $ca $c $ccn]]
      set omega [measure dihed [list $ca $c $ccn $ccca]]
   }
   return [list $phi $psi $omega]
}

# alpha-helix folder
# Parameters:
# -----------
# - rbegin: absolute residue number at beginning of section to be folded
# - rend: absolute residue number at and of section to be folded
# - rterm: absolute residue number at end of fragment to which rbegin and rend belong;
#          residues greater than rend but less than or equal to rterm are not folded
# - molid: molecule id
proc PestiferCRot::fold_alpha { rbegin rend rterm molid } {
   if { $rbegin < $rend } { # folding from N->C along section
      vmdcon -info "fold_alpha from $rbegin to $rend including fragment up to $rterm on mol $molid"
      set rotators_side C
      set increment 1
      proc finished { r rend } {
         return [expr $r > $rend]
      }
   } else { # folding from C->N along section
      vmdcon -info "fold_alpha from $rbegin to $rend including fragment up to $rterm on mol $molid"
      set rotators_side N
      set increment -1
      proc finished { r rend } {
         return [expr $r < $rend]
      }
   }
   for { set r $rbegin } { [finished $r $rend] == 0 } { incr r $increment } {
      set pp [get_phi_psi_omega $r $molid]
      set phi [lindex $pp 0]
      set psi [lindex $pp 1]
      set omega [lindex $pp 2]
      vmdcon -info "fold_alpha PRE residue $r ($phi,$psi,$omega)"
      if { [string compare $phi "NaN"] != 0 } {
         set dphi [expr (-57.0 - ($phi))]
         if { $rotators_side == "N" } {
            vmdcon -info "dphi negated"
            set dphi [expr -1*$dphi]
         }
         vmdcon -info "-> dphi $dphi"
         brot $molid $r $rterm phi $rotators_side $dphi
      }
      if {  [string compare $psi "NaN"] != 0 } {
         set dpsi [expr (-47.0 - ($psi))]
         if { $rotators_side == "N" } {
            set dpsi [expr -1*$dpsi]
         }
         vmdcon -info "-> dpsi $dpsi"
         brot $molid $r $rterm psi $rotators_side $dpsi
      }
      if { [string compare $omega "NaN"] != 0} {
         set domega [expr (180.0 - ($omega))]
         if { $rotators_side == "N" } {
            set domega [expr -1*$domega]
         }
         vmdcon -info "-> domega $domega"
         brot $molid $r $rterm omega $rotators_side $domega
      }
      set pp [get_phi_psi_omega $r $molid]
      set phi [lindex $pp 0]
      set psi [lindex $pp 1]
      set omega [lindex $pp 2]
      vmdcon -info "fold_alpha POST residue $r ($phi,$psi,$omega)"
   }
}

# Generalized bond rotation for proteins
# Parameters
# ----------
# - molid: molecule id
# - r0: absolute residue number of residue in which angle is located
# - r1: absolute residue number marking the terminus of the rotatable section
# - angle_name: one of phi, psi, or omega
# - rot: one of N (rotatable section is N-terminal to bond) or
#               C (rotatable section is C-terminal to bond)
# - deg: degrees of the rotation
#
# Cases:
#  r0<r1: rotatable section is C-terminal to bond; by necessity, rot==C
#  r0>r1: rotatable section is N-terminal to bond; by necessity, rot==N
#  r0==r1: rotatable section determined by rot
proc PestiferCRot::brot { molid r0 r1 angle_name rot deg} {
   set r_prev [expr $r0-1]
   set r_next [expr $r0+1]
   set c_prev -1
   set n_next -1
   set ca_next -1
   set sels [list]
   if {[residues_in_same_chain $molid $r0 $r_prev]} {
      set c_prev [atomselect $molid "residue $r_prev and name C"]
      lappend sels $c_prev
   }
   set n [atomselect $molid "residue $r0 and name N"]
   lappend sels $n
   set ca [atomselect $molid "residue $r0 and name CA"]
   lappend sels $ca
   set c [atomselect $molid "residue $r0 and name C"]
   lappend sels $c
   if {[residues_in_same_chain $molid $r0 $r_next]} {
      set n_next [atomselect $molid "residue $r_next and name N"]
      lappend sels $n_next
      set ca_next [atomselect $molid "residue $r_next and name CA"]
      lappend sels $ca_next
   }
   if { $angle_name == "phi" } {
      set pn [lindex [$n get {x y z}] 0]
      set pc [lindex [$ca get {x y z}] 0]
      if { $rot == "C" } { # rotators are c-terminal to the N-CA bond
         set rotators [atomselect $molid "(residue $r0 and not (name N HN HT1 HT2 HT3) ) or (residue > $r0 and residue <= $r1)"]
      } elseif { $rot == "N" } { # rotators are n-terminal to the N-CA bond
         set rotators [atomselect $molid "(residue $r0 and (name N HN HT1 HT2 HT3)) or (residue < $r0 and residue >= $r1)"]
         vmdcon -info "N-term phi brot $r0 to $r1 includes [$rotators num] atoms"
      }
   } elseif { $angle_name == "psi" } {
      set pn [lindex [$ca get {x y z}] 0]
      set pc [lindex [$c get {x y z}] 0]
      if { $rot == "C" } { # rotators are c-terminal to the CA-C bond
         set rotators [atomselect $molid "(residue $r0 and (name C O OT1 OT2 OXT) ) or (residue > $r0 and residue <= $r1)"]
      } elseif { $rot == "N" } { # rotators are n-terminal the CA-C bond
         set rotators [atomselect $molid "(residue $r0 and not (name C O OT1 OT2 OXT)) or (residue < $r0 and residue >= $r1)"]
         vmdcon -info "N-term psi brot $r0 to $r1 includes [$rotators num] atoms"
      }
   } elseif { $angle_name == "omega" } {
      set pn [lindex [$c get {x y z}] 0]
      set pc [lindex [$n_next get {x y z}] 0]
      if { $rot == "C" } { # rotators are c-terminal to C=N bond
         set rotators [atomselect $molid "residue > $r0 and residue <= $r1"]
      } elseif { $rot == "N" } { # rotators are n-terminal the C=N bond
         set rotators [atomselect $molid "residue <= $r0 and residue >= $r1"]
      }
   } elseif { $angle_name == "chi" } {
      if { $rot == 1 } {
         set pn [lindex [$ca get {x y z}] 0]
         set pc [lindex [[atomselect $molid "residue $r0 and name CB"] get {x y z}] 0]
         set rotators [atomselect $molid "residue $r0 and not name N HN CA C O"]
      } elseif { $rot == 2 } {
         set pn [lindex [[atomselect $molid "residue $r0 and name CB"] get {x y z}] 0]
         set pc [lindex [[atomselect $molid "residue $r0 and name CG"] get {x y z}] 0]
         set rotators [atomselect $molid "residue $r0 and not name N HN CA CB C O"]
      }
   } else {
      vmdcon -error "angle name $angle_name not recognized"
   }
   lappend sels $rotators

   $rotators move [trans bond $pn $pc $deg degrees]

   foreach s $sels {
      $s delete
   }
}

# Chain rotation procedures




# rotates all atoms in chain c-terminal to residue r up to and 
# including residue rend in chain c around residue r's phi angle 
# by deg degrees in molecule with id molid.  Recall phi:  (-C)-N--CA-C
# Note:  the key in this procedure is the definition of the atomselection
# containing the atoms that will rotate.  We do NOT include
# in this set atoms N or HN in resid "r" nor atoms "C" and "O" in 
# resid "rend".  The former is because they are N-terminal to the
# phi bond, and the latter is because psfgen seems to build 
# the C and O atom positions based on the positions of N, CA, and CB
# in the NEXT residue, rather than on the internal coordinates
# of the built residue.  That means, when we use psfgen to build
# the raw loop from residue i to j inclusive, it is the CA-C bond
# residue j that is non-natural in length.
proc PestiferCRot::Crot_phi { r rend segname molid deg } {
   set rot [atomselect $molid "((residue $r and not name N and not name HN) or (residue > $r and residue < $rend) or (residue $rend and not name C and not name O)) and segname $segname"]; checknum [$rot num] "no rot in Crot_phi";
   set n [atomselect $molid "residue $r and name N"] ; checknum [$n num] "No N in Crot_phi";
   set ca [atomselect $molid "residue $r and name CA"]; checknum [$ca num] "No CA in Crot_phi";
   set p1 [lindex [$n get {x y z}] 0]
   set p2 [lindex [$ca get {x y z}] 0]
   set ax [vecsub $p1 $p2]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $ca delete
   $n delete
}

# same as above, but treats last residue as C-terminus
proc PestiferCRot::Crot_phi_toCterm { r rend segname molid deg } {
   set rot [atomselect $molid "((residue $r and not name N and not name HN) or (residue > $r and residue <= $rend)) and segname $segname"]; checknum [$rot num] "no rot in Crot_phi";
   set n [atomselect $molid "residue $r and name N"] ; checknum [$n num] "No N in Crot_phi";
   set ca [atomselect $molid "residue $r and name CA"]; checknum [$ca num] "No CA in Crot_phi";
   set p1 [lindex [$n get {x y z}] 0]
   set p2 [lindex [$ca get {x y z}] 0]
   set ax [vecsub $p2 $p1]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $ca delete
   $n delete
}

proc PestiferCRot::Crot_psi { r rend segname molid deg } {
   set rot [atomselect $molid "((residue $r and name O) or (residue > $r and residue < $rend) or (residue $rend and not name C and not name O)) and segname $segname"]; checknum [$rot num] "no rot in Crot_psi";
   set ca [atomselect $molid "residue $r and name CA"]; checknum [$ca num] "No CA in Crot_psi";
   set co [atomselect $molid "residue $r and name C"] ; checknum [$co num] "No C in Crot_psi";
   set p1 [lindex [$ca get {x y z}] 0]
   set p2 [lindex [$co get {x y z}] 0]
   set ax [vecsub $p1 $p2]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $ca delete
   $co delete
}

proc PestiferCRot::Crot_psi_toCterm { r rend segname molid deg } {
   set rot [atomselect $molid "((residue $r and name O) or (residue > $r and residue <= $rend)) and segname $segname"]; checknum [$rot num] "no rot in Crot_psi";
   set ca [atomselect $molid "residue $r and name CA"]; checknum [$ca num] "No CA in Crot_psi";
   set co [atomselect $molid "residue $r and name C"] ; checknum [$co num] "No C in Crot_psi";
   set p1 [lindex [$ca get {x y z}] 0]
   set p2 [lindex [$co get {x y z}] 0]
   set ax [vecsub $p2 $p1]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $ca delete
   $co delete
}

# rotate the side chain of global residue r in mol molid around
# chi1 by deg degrees
proc PestiferCRot::SCrot_chi1 { r molid deg } {
   set rot [atomselect $molid "residue $r and not name N HN CA CB HA C O"]; checknum [$rot num] "no rot in SCrot_chi1";
   set ca [atomselect $molid "residue $r and name CA"]; checknum [$ca num] "no CA in SCrot_chi1";
   set cb [atomselect $molid "residue $r and name CB"]; checknum [$cb num] "no CB in SCrot_chi1";
   set p1 [lindex [$ca get {x y z}] 0]
   set p2 [lindex [$cb get {x y z}] 0]
   set ax [vecsub $p1 $p2]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $ca delete
   $cb delete
}

# rotate the side chain of residue r of chain c in mol molid around
# chi2 by deg degrees
proc PestiferCRot::SCrot_chi2 { r molid deg } {
   set rot [atomselect $molid "residue $r and not name N HN CA CB HA C O HB1 HB2"]; checknum [$rot num] "no rot in SCrot_chi1";
   set cb [atomselect $molid "residue $r and name CB"]; checknum [$cb num] "no CB in SCrot_chi1";
   set cg [atomselect $molid "residue $r and name CG"]; checknum [$cg num] "no CG in SCrot_chi1";
   set p1 [lindex [$cb get {x y z}] 0]
   set p2 [lindex [$cg get] {x y z}] 0]
   set ax [vecsub $p1 $p2]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $cb delete
   $cg delete
}

proc PestiferCRot::Crot_psi_special { r rend segname molid deg } {
   set rot [atomselect $molid "((residue > $r and residue < $rend) or (residue $rend and not name C and not name O)) and segname $segname"]; checknum [$rot num] "no rot in Crot_psi";
   set ca [atomselect $molid "residue $r and name CA"]; checknum [$ca num] "No CA in Crot_psi";
   set co [atomselect $molid "residue $r and name C"] ; checknum [$co num] "No C in Crot_psi";
   set p1 [lindex [$ca get {x y z}] 0]
   set p2 [lindex [$co get {x y z}] 0]
   set ax [vecsub $p1 $p2]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $ca delete
   $co delete
}

# rotates all residues C-terminal to "r" up to and including "rend" 
# (EXCEPT NOT C or O in rend) around the peptide bond between "r" and
# "r"+1.
proc PestiferCRot::Crot_omega { r rend segname molid deg } {
   set rot [atomselect $molid "((residue > $r and residue < $rend) or (residue $rend and not name C and not name O)) and segname $segname"]; checknum [$rot num] "no rot in Crot_omega";
   set co [atomselect $molid "residue $r and name C"] ; checknum [$co num] "No C in Crot_omega";
   set n [atomselect $molid "residue [expr $r + 1] and name N"] ; checknum [$n num] "No N in Crot_omega";
   set p1 [lindex [$co get {x y z}] 0]
   set p2 [lindex [$n get {x y z}] 0]
   set ax [vecsub $p1 $p2]
   $rot move [trans center $p1 axis $ax $deg degrees]
   $rot delete
   $n delete
   $co delete
}

# fold the run of new residues from r to rend into an alpha helix
# each phi,psi is 180 deg
# for alpha helix, phi ~ -60 deg, psi ~ -50
# Crot angles for phi are -120 and for psi -130
proc PestiferCRot::new_alpha { rbegin rend segname molid } {
   for { set r $rbegin } { $r <= $rend } { incr r } {
      # puts "phi/psi of $segname $r to alpha helix"
      Crot_phi_toCterm $r $rend $segname $molid -120
      if { $r < $rend } {
         Crot_psi_toCterm $r $rend $segname $molid -130
      }
   }
}

