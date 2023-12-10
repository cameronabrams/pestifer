# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# TcL procedures for facilitating rotations around bonds

# Wrap x into periodic domain [xlo,xhi]
# proc wrap_domain { x xlo xhi } {
#    set dsize [expr ($xhi) - ($xlo)]
#    if { [expr $x < $xlo] } {
#       set y [expr $x + $dsize]
#    } elseif { [expr $x > $xhi] } {
#       set y [expr $x - $dsize]
#    } else {
#       set y $x
#    }
#    return $y
# }

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
proc residues_in_same_chain { molid q r } {
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
proc get_phi_psi_omega { r molid } {
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
proc fold_alpha { rbegin rend rterm molid } {
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
         # set dphi [wrap_domain [expr (-57 - ($phi))] -180.0 180.0]
         vmdcon -info "-> dphi $dphi"
         brot $molid $r $rterm phi $rotators_side $dphi
      }
      if {  [string compare $psi "NaN"] != 0 } {
         # set dpsi [wrap_domain [expr (-47 - ($psi))] -180.0 180.0]
         set dpsi [expr (-47.0 - ($psi))]
         vmdcon -info "-> dpsi $dpsi"
         brot $molid $r $rterm psi $rotators_side $dpsi
      }
      if { [string compare $omega "NaN"] != 0} {
         # set domega [wrap_domain [expr (180 - ($omega))] -180.0 180.0]
         set domega [expr (180.0 - ($omega))]
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


# Generalized bond rotation for terminal groups
# Parameters
# ----------
# - molid : molecule id
# - atomsel : atomselection of the group containing the bond to rotate
# - a[0,1]idx: atom indices of the two bonded atoms defining the bond to rotate
# - deg: degrees of rotation around bond [0,360]
# proc trot { molid atomsel a0idx a1idx deg} {
#    set all_ids [$atomsel get index]
#    set bonds [$atomsel getbonds]
#    # TODO: identify subset of serials including a1ser and all atoms "downstream"
#    # from it in the terminal group
#    # to do this, tag all atoms whose bonded neighbors are *within* the selection
#    set all_internal_ids [list]
#    foreach ao $all_ids bo $bonds {
#       set is_internal 1
#       foreach ba $bo {
#          set se [lsearch $all_ids $ba]
#          if { se == -1 } {
#             set is_internal 0
#          }
#       }
#       if { $is_internal == 1} {
#          lappend all_internal_ids $ao
#       }
#    }

# }

# Align an atomselection by overlapping one triplet of atoms in that
# atomselection with a given triplet of atoms
#
# Parameters
# ----------
#
# - to_sel: atomselection
#      The atomselection containing the ordered triplet of atoms to which
#      the alignment is made
# - from_sel: atomselection
#      The atomselection containing the ordered triplet of atoms from
#      which the alignment is made
#  - mover_sel: atomselection
#      The atomselection of all atoms that are moved
#
proc hatme { a_vec } {
   set l [veclength $a_vec]
   return [vecscale $a_vec [expr 1.0/$l]]
}
proc align_triples {to_sel from_sel mover_sel} {
   set to_   [$to_sel   get {x y z}]
   set from_ [$from_sel get {x y z}]
   # translate mover to overlap root atoms
   set dr [vecsub [lindex $to_ 0] [lindex $from_ 0]]
   $mover_sel moveby $dr
   set from_ [$from_sel get {x y z}]
   # calculate vectors from root to first in the align_to and the mover
   set to_vec     [vecsub [lindex $to_   0] [lindex $to_   1]]
   set from_vec   [vecsub [lindex $from_ 0] [lindex $from_ 1]]
   set to_uvec    [hatme $to_vec  ]
   set from_uvec  [hatme $from_vec]
   set axis       [veccross $to_uvec $from_uvec]
   set ctheta     [vecdot $to_uvec $from_uvec]
   set theta      [expr -180.0/3.141593*acos($ctheta)]
   set trans      [trans bond [lindex $to_ 0] [vecadd [lindex $to_ 0] $axis] $theta deg]
   $mover_sel move $trans
   # project to two right points into plane orthogonal to 0->1 vector, compute angle between
   # projected bonds to compute dihedral angle required to rotate the mover right to
   # the to right
   set from_      [$from_sel get {x y z}]
   set from_vec   [vecsub [lindex $from_ 0] [lindex $from_ 1]]
   set from_uvec  [hatme $from_vec]
   set to_p       [lindex $to_   2]
   set from_p     [lindex $from_ 2]
   set to_pp      [vecsub $to_p   [vecscale [vecdot $to_p   $to_uvec]   $to_uvec  ]]
   set from_pp    [vecsub $from_p [vecscale [vecdot $from_p $from_uvec] $from_uvec]]
   set to_opp     [vecsub [lindex $to_ 0] $to_pp]
   set to_uopp    [hatme $to_opp]
   set from_opp   [vecsub [lindex $from_ 0] $from_pp]
   set from_uopp  [hatme $from_opp]
   set ctheta     [vecdot $to_uopp $from_uopp]
   set theta      [expr 180.0/3.141593*acos($ctheta)]
   set trans      [trans bond [lindex $to_ 0] [lindex $to_ 1] $theta deg]
   $mover_sel move $trans
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
proc brot { molid r0 r1 angle_name rot deg} {
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
      }
   } elseif { $angle_name == "psi" } {
      set pn [lindex [$ca get {x y z}] 0]
      set pc [lindex [$c get {x y z}] 0]
      if { $rot == "C" } { # rotators are c-terminal to the CA-C bond
         set rotators [atomselect $molid "(residue $r0 and (name C O OT1 OT2 OXT) ) or (residue > $r0 and residue <= $r1)"]
      } elseif { $rot == "N" } { # rotators are n-terminal the CA-C bond
         set rotators [atomselect $molid "(residue $r0 and not (name C O OT1 OT2 OXT)) or (residue < $r0 and residue >= $r1)"]
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
         set pc [lindex [[atomselect $molid "residue and name CB"] get {x y z}] 0]
         set rotators [atomselect $molid "residue $r0 and not name N HN CA C O"]
      } elseif { $rot == 2 } {
         set pn [lindex [[atomselect $molid "residue and name CB"] get {x y z}] 0]
         set pc [lindex [[atomselect $molid "residue and name CG"] get {x y z}] 0]
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
proc checknum { num msg } {
  if { $num == 0 } { 
    puts "$msg"
    #exit
  }
}

# Given atom indices a and b of a bonded pair of atoms, a list of atom indices 
# determining an atomselection and the associated bondlist, accumulate
# indices on the "a-side" of the a-b bond, assuming that cutting that bond
# would yield two disconnected fragments
proc half_traversal {a b iL bL aG} {
   set ai [lsearch $iL $a]
   if {$ai==-1} { return $aG }
   set ab [lindex $bL $ai]
   foreach bi_ $ab {
      if {[lsearch $iL $bi_]!=-1 && $bi_ != $b} {
         if {[lsearch $aG $bi_]==-1} {
            lappend aG $bi_
            set aG [half_traversal $bi_ $a $iL $bL $aG]
         }
      }
   }
   return $aG
}

proc sel_is_pendant { atomsel } {
   vmdcon -info "checking sel of [$atomsel num] atoms for pendancy"
   set indexlist [$atomsel get index]
   set bondlist  [$atomsel getbonds]
   set outers [list]
   foreach bonds $bondlist {
      foreach bond $bonds {
         if {[lsearch $indexlist $bond]==-1} {
            if {[lsearch $outers $bond]==-1} {
               lappend outers $bond
            }
         }
      }
   }
   vmdcon -info " -> [llength $outers] bond(s) to atoms outside selection"
   if {[llength $outers]==1} {
      return 1
   } else {
      return 0
   }
}
# Determine indices of atoms that move to execute
# rotation around bond within atomsel
proc determine_movers_on_rotation { bond atomsel } {
   set movers [list]
   set indexlist [$atomsel get index]
   set bondlist  [$atomsel getbonds]

   set a [lindex $bond 0]
   set b [lindex $bond 1]
   # traverse each atoms's half-tree gathering indices.  
   # only one traversal should yield an index to a bonded
   # neighbor that is NOT in the atomselection, meaning
   # this is the pendant atomselection's connection atom
   # in the rest of the molecule.  The atoms for which
   # the half-tree traversal yields the connector is 
   # roots the half that is non-rotating, meaning all
   # other atoms in the atomsel are movers.
   set a_gather $bond
   set a_gather [half_traversal $a $b $indexlist $bondlist $a_gather]
   set a_is_left 0
   foreach a_ $a_gather {
      set findme [lsearch $indexlist $a_]
      if {$findme == -1} {
         set a_is_left 1
      }
   }
   foreach i_ $indexlist {
      set findme [lsearch $a_gather $i_]
      if {$findme == -1} {
         lappend movers $i_
      }
   }
   return $movers
}

proc get_rotatable_bonds { atomsel molid } {
   set element [$atomsel get element]
   set name [$atomsel get name]
   set iL [$atomsel get index]
   set bL [$atomsel getbonds]
   set r5sel [atomselect $molid "ringsize 5 from ([$atomsel text])"]
   set r6sel [atomselect $molid "ringsize 6 from ([$atomsel text])"]
   set r5i [$r5sel get index]
   set r6i [$r6sel get index]
   set rbonds [list]
   for {set aidx 0} {$aidx < [llength $iL]} {incr aidx} {
      set ai [lindex $iL $aidx]
      set a_el [lindex $element $aidx]
      if { [string equal $a_el 'H'] } {
         continue
      }
      set a_nm [lindex $name $aidx]
      set B [lindex $bL $aidx]
      foreach bi $B {
         set bidx [lsearch $iL $bi]
         if { $bidx == -1 } { 
            # assume pendant connection is a rotatable bond
            set bo [lsort -integer [list $ai $bi]]
            # vmdcon -info "Adding $bo"
            lappend rbonds $bo
            set rbonds [lsort -unique $rbonds]
         } else {
            set b_el [lindex $element $bidx]
            if { [string equal $b_el 'H'] } {
               continue
            }
            set b_nm [lindex $name $bidx]
            # vmdcon -info "a $aidx $ai $a_el $a_nm -- b $bidx $bi $b_el $b_nm"
            set checkP1 [expr !(([string equal $a_nm 'N'])&&([string equal $b_nm 'C']))]
            set checkP2 [expr !(([string equal $a_nm 'C'])&&([string equal $b_nm 'N']))]
            set checkR5 [expr (([lsearch $r5i $ai]==-1)||([lsearch $r5i $bi]==-1))]
            set checkR6 [expr (([lsearch $r6i $ai]==-1)||([lsearch $r6i $bi]==-1))]
            # vmdcon -info " --> $checkP1 $checkP2 $checkR5 $checkR6"
            if {($checkP1)&&($checkP2)&&($checkR5)&&($checkR6)} {
               set bo [lsort -integer [list $ai $bi]]
               # vmdcon -info "Adding $bo"
               lappend rbonds $bo
               set rbonds [lsort -unique $rbonds]
               # vmdcon -info "There are now [llength $rbonds] rotatable bonds"
            }
         }
      }
   }
   return $rbonds
}
proc declash_pendant_sel { atomsel molid maxcycles } {
   # if { [sel_is_pendant $atomsel]==0 } {
   #    vmdcon -info "Selection from $molid via [$atomsel text] is not pendant"
   #    return
   # }
   set degs {-120, 120}
   set environ [atomselect $molid "not ([$atomsel text])"]
   set rbonds [get_rotatable_bonds $atomsel $molid]
   set nbonds [llength $rbonds]
   vmdcon -info "Declash pendant: Pendant has [$atomsel atoms] and $nbonds rotatable bonds"
   set ncontacts [llength [lindex [measure contacts 1.0 $atomsel $environ] 0]]
   vmdcon -info "Declash pendant: Pendant has $ncontacts initial atomic clashes"
   vmdcon -info "Declashing via maximally $maxcycles cycles"
   for {set i 0} {$i < $maxcycles} {incr i} {
      set ridx [expr {int(rand()*$nbonds)}]
      set bo [lindex $rbonds $ridx]
      set b [[atomselect $molid "index $bo"] get {x y z}]
      set moveridx [determine_movers_on_rotation $bo $atomsel]
      set movers [atomselect $molid "index $moveridx"]
      set posn [backup $movers {x y z}]
      set didx [expr {int(rand()*2)}]
      set deg [lindex $degx $didx]
      set tmat [trans center [lindex $b 0] bond [lindex $b 0] [lindex $b 1] $deg degrees]
      $movers move $tmat
      set newcontacts [llength [lindex [measure contacts 1.0 $atomsel $environ] 0]]
      if { $newcontacts >= $ncontacts } {
         restore $movers $posn
      } else {
         set ncontacts $newcontacts
      }
      vmdcon -info "  cycle $i bond $b deg $deg ncontacts $ncontacts"
   }
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

