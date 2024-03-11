# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Utility procs for VMD/Psfgen TcL scripts
#

package provide PestiferUtil 1.0

namespace eval ::PestiferUtil:: {
    namespace export *
}

# Logs a message to both the VMD console and to an open
# log file
proc PestiferUtil::double_log { logf message } {
   vmdcon -info "$message"
   puts $logf "$message"
   flush $logf
}

# outputs the message if num is zero
proc PestiferUtil::checknum { num msg } {
  if { $num == 0 } { 
    vmdcon -info "$msg"
  }
}

# returns -1 if x is negative, 1 otherwise
proc PestiferUtil::sign { x } {
    if { $x < 0 } {
        return -1
    } else {
        return 1
    }
}

# returns the unit vector of an input vector
proc PestiferUtil::hatvec { a_vec } {
   return [vecnorm $a_vec]
}

# Wrap x into periodic domain [xlo,xhi]
proc PestiferUtil::wrap_domain { x xlo xhi } {
   set dsize [expr ($xhi) - ($xlo)]
   if { [expr $x < $xlo] } {
      set y [expr $x + $dsize]
   } elseif { [expr $x > $xhi] } {
      set y [expr $x - $dsize]
   } else {
      set y $x
   }
   return $y
}
proc PestiferUtil::getcellvec { pad } {
   if {![info exists pad]} {
      set pad 0.1
   }
   set sel [atomselect top all]
   $sel frame 0
   set mm [measure minmax $sel]
   set ll [lindex $mm 0]
   set ur [lindex $mm 1]
   set sp [vecsub $ur $ll]
   # puts "cellBasisVector1   [expr [lindex $sp 0] + $pad] 0.0 0.0"
   # puts "cellBasisVector2   0.0 [expr [lindex $sp 1] + $pad] 0.0"
   # puts "cellBasisVector3   0.0 0.0 [expr [lindex $sp 2] + $pad]"
   # puts "cellOrigin         [measure center $sel]"
   return "cellBasisVector1   [expr [lindex $sp 0] + $pad] 0.0 0.0\ncellBasisVector2   0.0 [expr [lindex $sp 1] + $pad] 0.0\ncellBasisVector3   0.0 0.0 [expr [lindex $sp 2] + $pad]\ncellOrigin         [measure center $sel]"
}

proc PestiferUtil::update_atomselect_macro { macro_name new_value overwrite } {
   catch { atomselect macro $macro_name } results options
   set code [dict get $options -code]
   set existing_macro ""
   if { $code == 0 } {
      set existing_macro $results
   }
   if {$overwrite == 1} {
      atomselect macro $macro_name $new_value
   } else {
      if {$existing_macro != ""} {
         atomselect macro $macro_name "($existing_macro) or ($new_value)"
      } else {
         atomselect macro $macro_name $new_value
      }
   }
}

proc PestiferUtil::resequence { badsel goodsel } {
   if {[$badsel num] != [$goodsel num]} {
      vmdcon -info "Cannot resequence; selections have differing numbers of atoms"
      return $badsel
   }
   set badresids [lsort -unique -integer [$badsel get resid]]
   set goodresids [lsort -unique -integer [$goodsel get resid]]
   foreach br $badresids gr $goodresids {
      if { $br != $gr } {
         vmdcon info "Cannot resequence; selections do not have the same set of resids"
         return $badsel
      }
   }
   set badidx [$badsel get index]
   set goods [$goodsel get {resname resid insertion name}]
   set goodpos [$goodsel get {x y z}]
   set badsav [$badsel get {resname resid insertion name}]
   set badpos [$badsel get {x y z}]

   set badovr [list]
   foreach g $goods gp $goodpos ba $badsav {
      set resid [lindex $g 1]
      set insertion [lindex $g 2]
      set name [lindex $g 3]
      for { set idx 0 } { $idx < [llength $badsav] } { incr idx } {
         set b [lindex $badsav $idx]
         if { $resid == [lindex $b 1] && $insertion == [lindex $b 2] && $name == [lindex $b 3] } {
            vmdcon -info "bad atom $ba overwritten by good atom $g"
            lappend badovr $gp
            break
         }
      }
   }
   if { [llength $badovr] != [llength $goods ] } {
      vmdcon -info "Error: bad selection could not inherit all good data"
      return $badsel
   }
   return $badpos
}

proc PestiferUtil::letter_up { c } {
   scan $c %c i
   return [binary format c* [expr $i + 1]]
}

proc PestiferUtil::renumber_as_insertion { atomsel root_resid root_insertion } {
    set orig_resid [$atomsel get resid]
    set orig_residue [$atomsel get residue]
    set residues [lsort -unique $orig_residue]
    set residue_insertion [list]
    set orig_insertion [$atomsel get insertion]
    set new_resid [list]
    set new_insertion [list]
    vmdcon -info "renumber_as_insertion on [$atomsel num] atoms"
    foreach ori $orig_resid {
        lappend new_resid $root_resid
    }
    # { } is blank insertion
    if { $root_insertion == { } } {
        set ins "A"
        set icode [scan $ins "%c"]
    } else {
        set code [scan $root_insertion "%c"]
        set icode [expr $code + 1]
        set ins [format "%c" $icode]
    }
    foreach residue $residues {
        lappend residue_insertion $ins
        incr icode
        set ins [format "%c" $icode]
    }
    vmdcon -info "start at $root_resid $ins"
    for {set i 0} {$i < [$atomsel num]} {incr i} {
        set this_residue [lindex $orig_residue $i]
        set ridx [lsearch $residues $this_residue]
        set this_insertion [lindex $residue_insertion $ridx]
        lappend new_insertion $this_insertion
    }
    
    $atomsel set resid $new_resid
    $atomsel set insertion $new_insertion
}