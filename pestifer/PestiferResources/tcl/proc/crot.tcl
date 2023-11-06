# Chain rotation procedures
proc checknum { num msg } {
  if { $num == 0 } { 
    puts "$msg"
    #exit
  }
}

proc get_phi_psi { r segname molid } {
   set nn [atomselect $molid "segname $segname and residue [expr $r - 1]"]
   set phi "NaN"
   if { [$nn num] > 0 } {
      set nnc [[atomselect $molid "segname $segname and residue [expr $r - 1] and name C"] get index]
      set n [[atomselect $molid "segname $segname and residue $r and name N"] get index]
      set ca [[atomselect $molid "segname $segname and residue $r and name CA"] get index]
      set c  [[atomselect $molid "segname $segname and residue $r and name C"] get index]
      set phi [measure dihed [list $nnc $n $ca $c]]
   }
   set cc [atomselect $molid "segname $segname and residue [expr $r + 1]"]
   set psi "NaN"
   if { [$nn num] > 0 } {
      set n [[atomselect $molid "segname $segname and residue $r and name N"] get index]
      set ca [[atomselect $molid "segname $segname and residue $r and name CA"] get index]
      set c  [[atomselect $molid "segname $segname and residue $r and name C"] get index]
      set ccn [[atomselect $molid "segname $segname and residue [expr $r + 1] and name N"] get index]
      set psi [measure dihed [list $n $ca $c $ccn]]
   }
   return [list $phi $psi]
}

# r0 and r1 are ABSOLUTE RESIDUE NUMBERS
proc brot { molid r0 r1 angle_name rot deg} {

   set cn -1
   if { $r0 > 0 } {
      set cn [atomselect $molid "residue [expr $r0-1] and name C"]; checknum [$cn num] "No CN in brot";
   }
   set n [atomselect $molid "residue $r0 and name N"] ; checknum [$n num] "No N in brot";
   set ca [atomselect $molid "residue $r0 and name CA"]; checknum [$ca num] "No CA in brot";
   set c [atomselect $molid "residue $r0 and name C"] ; checknum [$n num] "No C in brot";
   set nc [atomselect $molid "residue [expr $r0+1] and name N"] ; checknum [$nc num] "No NC in brot";

   set center -1
   set ax -1

   if { $angle_name == "phi" } {
      set pn [lindex [$n get {x y z}] 0]
      set pc [lindex [$ca get {x y z}] 0]
      if { $rot == "C" } {
         set ax [vecsub $pn $pc]
         set center $pn
         if { $r0 < $r1 } {
            set rotators [atomselect $molid "(residue $r0 and not (name N HN HT1 HT2 HT3) ) or (residue > $r0 and residue <= $r1)"]
         } elseif { $r0 == $r1 } {
            set rotators [atomselect $molid "residue $r0 and not (name N HN HT1 HT2 HT3)"]
         } else {
            puts "Error: You are asking for a C-rot at residue $r0 up to residue $r1 but $r1 is N-terminal to $r0"
         }
      } elseif { $rot == "N" } {
         set ax [vecsub $pc $pn]
         set center $pc
         if { $r0 > $r1 } {
            set rotators [atomselect $molid "(residue $r0 and (name N HN HT1 HT2 HT3)) or (residue < $r0 and residue >= $r1)"]
         } elseif { $r0 == $r1 } {
            set rotators [atomselect $molid "residue $r0 and (name N HN HT1 HT2 HT3)"]
         } else {
            puts "Error: You are asking for a N-rot at residue $r0 back to residue $r1 but $r1 is C-terminal to $r0"
         }
      }
   } elseif { $angle_name == "psi" } {
      set pn [lindex [$ca get {x y z}] 0]
      set pc [lindex [$c get {x y z}] 0]
      if { $rot == "C" } {
         set ax [vecsub $pn $pc]
         set center $pn
         if { $r0 < $r1 } {
            set rotators [atomselect $molid "(residue $r0 and (name C O OT1 OT2 OXT) ) or (residue > $r0 and residue <= $r1)"]
         } elseif { $r0 == $r1 } {
            set rotators [atomselect $molid "residue $r0 and (name C O OT1 OT2 OXT)"]
         } else {
            puts "Error: You are asking for a C-rot at residue $r0 up to residue $r1 but $r1 is N-terminal to $r0"
         }
      } elseif { $rot == "N" } {
         set ax [vecsub $pc $pn]
         set center $pc
         if { $r0 > $r1 } {
            set rotators [atomselect $molid "(residue $r0 and not (name C O OT1 OT2 OXT)) or (residue < $r0 and residue >= $r1)"]
         } elseif { $r0 == $r1 } {
            set rotators [atomselect $molid "residue $r0 and not (name C O OT1 OT2 OXT)"]
         } else {
            puts "Error: You are asking for a N-rot at residue $r0 back to residue $r1 but $r1 is C-terminal to $r0"
         }
      }
   } elseif { $angle_name == "omega" } {
      puts "Error: omega rotations not yet implemented."
   } else {
      puts "Error: angle name $angle_name not recognized"
   }
   if { $center != -1 && $ax != -1 } {
      $rotators move [trans center $center axis $ax $deg degrees]
   } else {
      puts "Error: no rotation performed"
   }
   $rotators delete
   if { $cn != -1 } {
      $cn delete
   }
   $n delete
   $ca delete
   $c delete
   $nc delete
}

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
proc Crot_phi { r rend segname molid deg } {
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
proc Crot_phi_toCterm { r rend segname molid deg } {
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

proc Crot_psi { r rend segname molid deg } {
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

proc Crot_psi_toCterm { r rend segname molid deg } {
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
proc SCrot_chi1 { r molid deg } {
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
proc SCrot_chi2 { r molid deg } {
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

proc Crot_psi_special { r rend segname molid deg } {
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
proc Crot_omega { r rend segname molid deg } {
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
proc new_alpha { rbegin rend segname molid } {
   for { set r $rbegin } { $r <= $rend } { incr r } {
      # puts "phi/psi of $segname $r to alpha helix"
      Crot_phi_toCterm $r $rend $segname $molid -120
      if { $r < $rend } {
         Crot_psi_toCterm $r $rend $segname $molid -130
      }
   }
}

proc wrap_domain { x xlo xhi } {
   set dsize [expr $xhi - $xlo]
   if { [expr $x < $xlo] } {
      set y [expr $x + $dsize]
   } elseif { [expr $x > $xhi] } {
      set y [expr $x - $dsize]
   } else {
      set y $x
   }
   return $y
}

proc fold_alpha { rbegin rend segname molid } {
   for { set r $rbegin } { $r <= $rend } { incr r } {
      set pp [get_phi_psi $r $segname $molid]
      set phi [lindex $pp 0]
      set psi [lindex $pp 1]
      # puts "$r $phi $psi"
      # flush stdout
      if { $phi != "NaN "} {
         set dphi [wrap_domain [expr (-60 - ($phi))] -360.0 360.0]
         # puts "dphi $dphi"
         # flush stdout
         Crot_phi_toCterm $r $rend $segname $molid $dphi
      }
      if { $r < $rend } {
         if { $psi != "NaN"} {
            set dpsi [wrap_domain [expr (-50 - ($psi))] -360.0 360.0]
            # puts "dpsi $dpsi"
            # flush stdout
            Crot_psi_toCterm $r $rend $segname $molid $dpsi
         }
      }
   }
}