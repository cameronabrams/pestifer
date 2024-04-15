# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
package provide PestiferPierce 1.0

namespace eval ::PestiferPierce:: {
    namespace export *
}
# Check for pierced rings
#
# Parameters
# ----------
# molid - the molecule id
# ringsize - desired ringsize for identifying rings
# TOL - tolerance in A
#
# Note: this procedure does not use PBCs; it is really meant for
# checking pierced rings in proteins in which the entire protein
# is centered in a periodic image.  TODO: fix this
# 
proc PestiferPierce::check_pierced_rings { molid ringsize TOL } {
  # select all atoms in rings
  set r6 [atomselect $molid "ringsize $ringsize from all"]
  set r6i [$r6 get index]
  # make a mapping of atom index to internal index
  set i 0
  foreach ii $r6i {
    set r6o($ii) $i
    incr i
  }
  # get all x, y, and z positions of ring atoms; we will assume that rings are not
  # broken across PBC
  set r6x [$r6 get x]
  set r6y [$r6 get y]
  set r6z [$r6 get z]

  for { set ri 0 } { $ri < [llength $r6i] } { incr ri $ringsize } {
    #puts "ring $ri"
    set this_ri {}
    set this_rx {}
    set this_ry {}
    set this_rz {}
    for { set t 0 } { $t < $ringsize } { incr t } {
      lappend this_ri [lindex $r6i [expr $ri + $t]]
      lappend this_rx [lindex $r6x [expr $ri + $t]]
      lappend this_ry [lindex $r6y [expr $ri + $t]]
      lappend this_rz [lindex $r6z [expr $ri + $t]]
    }
    #puts "this_ri $this_ri"
    set this_rr {}
    foreach x $this_rx y $this_ry z $this_rz {
      lappend this_rr [list $x $y $z]
    }
    #puts "this_rr $this_rr"
    set this_com [list [vecsum $this_rx] [vecsum $this_ry] [vecsum $this_rz]]
    set this_com [vecscale $this_com [expr 1./$ringsize]]
    set this_b12 [vecsub [lindex $this_rr 0] [lindex $this_rr 1]]
    set this_b23 [vecsub [lindex $this_rr 1] [lindex $this_rr 3]]
    set c123 [veccross $this_b12 $this_b23]
    set lc123 [veclength $c123]
    set chat123 [vecscale $c123 [expr 1.0/$lc123]]
    #puts "ring $this_ri : $this_com : $chat123"
    set neigh [atomselect $molid "same residue as (exwithin 4.0 of index $this_ri)"]
    set nb [$neigh getbonds]
    set na [$neigh get index]
    set nax [$neigh get x]
    set nay [$neigh get y]
    set naz [$neigh get z]
    set i 0
    foreach a $na {
      set ord($a) $i
      incr i
    }
    set ln 0
    set ndots 0
    # fix to exclude atoms in the ring from bond definition!
    foreach a $na bl $nb {
      if { [lsearch $this_ri $a] !=  -1 } {
        continue
      }
      if { [expr $ln%100 == 0] } {
        puts -nonewline "."
        incr ndots
        if { $ndots > 80 } {
          puts ""
          set ndots 0
        }
        flush stdout
      }
      incr ln
      set ai $ord($a)
      set apos [list [lindex $nax $ai] [lindex $nay $ai] [lindex $naz $ai]]
      foreach b $bl {
        if { [lsearch $this_ri $b] != -1 } {
          continue
        }
        if { [ expr $b < $a ] } {
          continue
        }
        if { [lsearch $na $b] != -1 } {
          set bi $ord($b)
          set bpos [list [lindex $nax $bi] [lindex $nay $bi] [lindex $naz $bi]]
          set avpos [vecscale [vecadd $apos $bpos] 0.5]
          set avec [vecsub $this_com $apos]
          set bvec [vecsub $this_com $bpos]
          set crit1 [expr [veclength [vecsub $avpos $this_com]] < $TOL]
          set d1 [vecdot $avec $chat123]
          set d2 [vecdot $bvec $chat123]
          set crit2 [expr ($d1*$d2)<0]
          if { $crit1 && $crit2 } {
            puts ""
            set piercer [atomselect $molid "index $a $b"]
            set b_resnames [$piercer get resname]
            set b_names [$piercer get name]
            set b_segnames [$piercer get segname]
            set b_resids [$piercer get resid]
            set b_chain [$piercer get chain]
            set piercee [atomselect $molid "index $this_ri"]
            set first "[lindex $b_chain 0]([lindex $b_segnames 0])_[lindex $b_resnames 0][lindex $b_resids 0][lindex $b_names 0]($a)"
            set second "[lindex $b_chain 1]([lindex $b_segnames 1])_[lindex $b_resnames 1][lindex $b_resids 1][lindex $b_names 1]($b)"
            set ring "[lindex [$piercee get chain] 0]([lindex [$piercee get segname] 0])_[lindex [$piercee get resname] 0][lindex [$piercee get resid] 0]"
            puts "Bond $first-$second pierces $ring"
          }
        }
      }
    }
    $neigh delete
  }
}