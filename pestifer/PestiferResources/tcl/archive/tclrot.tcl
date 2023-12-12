# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Pestifer 1.2.9: unused TcL procs 
#
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
      # don't consider this atom if it only has one bonded neighbor
      if {[llength $B]<=1} { 
         continue
      }
      foreach bi $B {
         set bidx [lsearch $iL $bi]
         if { $bidx == -1 } { 
            # # assume pendant connection is a rotatable bond
            # set bo [lsort -integer [list $ai $bi]]
            # # vmdcon -info "Adding $bo"
            # lappend rbonds $bo
            # set rbonds [lsort -unique $rbonds]
         } elseif {[llength [lindex $bL $bidx]]>1} {
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

proc check_bonds {bonds} {
   foreach bo $bonds {
      if {[measure bond $bo]>2.0} {
         vmdcon -info "warning: bond $bo is too long"
      }
   }
}

proc declash_pendant_sel { atomsel molid maxcycles clashdist} {
   if { [sel_is_pendant $atomsel]==0 } {
      vmdcon -info "Selection from $molid via [$atomsel text] is not pendant"
      return
   }
   set degs {-120 120}
   set environ [atomselect $molid "all"]
   vmdcon -info "Declash environ has [$environ num] atoms"
   set rbonds [get_rotatable_bonds $atomsel $molid]
   set nbonds [llength $rbonds]
   vmdcon -info "Declash pendant: Pendant has [$atomsel num] atoms and $nbonds rotatable bonds"
   set ncontacts [llength [lindex [measure contacts $clashdist $atomsel $environ] 0]]
   vmdcon -info "Declash pendant: Pendant has $ncontacts initial atomic clashes"
   vmdcon -info "Declashing via maximally $maxcycles cycles"
   for {set i 0} {$i < $maxcycles} {incr i} {
      set ridx [expr {int(rand()*$nbonds)}]
      set bo [lindex $rbonds $ridx]
      set b [[atomselect $molid "index $bo"] get {x y z}]
      vmdcon -info " $i : rotating bond $bo ($ridx) $b"
      set moveridx [determine_movers_on_rotation $bo $atomsel]
      vmdcon -info "  movers $moveridx"
      set movers [atomselect $molid "index $moveridx"]
      # set posn [backup $movers {x y z}]
      set didx [expr {int(rand()*2)}]
      set deg [lindex $degs $didx]
      set tmat [trans center [lindex $b 0] bond [lindex $b 0] [lindex $b 1] $deg degrees]
      $movers move $tmat
      check_bonds $rbonds
      set newcontacts [llength [lindex [measure contacts $clashdist $atomsel $environ] 0]]
      if { $newcontacts >= $ncontacts } {
         set deg [expr -1 * ($deg)]
         set tmat [trans center [lindex $b 0] bond [lindex $b 0] [lindex $b 1] $deg degrees]
         $movers move $tmat
         # restore $movers {x y z} $posn
      } else {
         set ncontacts $newcontacts
      }
      vmdcon -info "  cycle $i bond $bo deg $deg ncontacts $ncontacts"
   }
}

# Align an atomselection by overlapping one triplet of atoms in that
# atomselection with a given triplet of atoms

# Parameters
# ----------

# - to_sel: atomselection
#      The atomselection containing the ordered triplet of atoms to which
#      the alignment is made
# - from_sel: atomselection
#      The atomselection containing the ordered triplet of atoms from
#      which the alignment is made
#  - mover_sel: atomselection
#      The atomselection of all atoms that are moved

proc align_triples {to_sel from_sel mover_sel} {
   set to_   [$to_sel   get {x y z}]
   set from_ [$from_sel get {x y z}]
   # translate mover to overlap root atoms
   set dr [vecsub [lindex $to_ 0] [lindex $from_ 0]]
   $mover_sel moveby $dr
   $from_sel moveby $dr
   set from_ [$from_sel get {x y z}]
   # calculate vectors from root to first in the align_to and the mover
   set to_vec     [vecsub [lindex $to_   0] [lindex $to_   1]]
   set from_vec   [vecsub [lindex $from_ 0] [lindex $from_ 1]]
   set to_uvec    [hatvec $to_vec  ]
   set from_uvec  [hatvec $from_vec]
   set axis       [veccross $to_uvec $from_uvec]
   set ctheta     [vecdot $to_uvec $from_uvec]
   set theta      [expr -180.0/3.141593*acos($ctheta)]
   set trans      [trans bond [lindex $to_ 0] [vecadd [lindex $to_ 0] $axis] $theta deg]
   $mover_sel move $trans
   $from_sel move $trans
   # project to two right points into plane orthogonal to 0->1 vector, compute angle between
   # projected bonds to compute dihedral angle required to rotate the mover right to
   # the to right
   set from_      [$from_sel get {x y z}]
   set from_vec   [vecsub [lindex $from_ 0] [lindex $from_ 1]]
   set from_uvec  [hatvec $from_vec]
   set to_p       [lindex $to_   2]
   set from_p     [lindex $from_ 2]
   set to_pp      [vecsub $to_p   [vecscale [vecdot $to_p   $to_uvec]   $to_uvec  ]]
   set from_pp    [vecsub $from_p [vecscale [vecdot $from_p $from_uvec] $from_uvec]]
   set to_opp     [vecsub [lindex $to_ 0] $to_pp]
   set to_uopp    [hatvec $to_opp]
   set from_opp   [vecsub [lindex $from_ 0] $from_pp]
   set from_uopp  [hatvec $from_opp]
   set ctheta     [vecdot $to_uopp $from_uopp]
   set theta      [expr 180.0/3.141593*acos($ctheta)]
   set trans      [trans bond [lindex $to_ 0] [lindex $to_ 1] $theta deg]
   $mover_sel move $trans
}
