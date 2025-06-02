# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Utility procs for VMD/Psfgen TcL scripts
#

package provide PestiferUtil 1.0

namespace eval ::PestiferUtil:: {
    namespace export *
}

proc PestiferUtil::split_psf {psf pdb residues {prefix "section"}} {
   set nsections [llength $residues]
   for {set i 0} {$i < $nsections} {incr i} {
      vmdcon -info "split_psf: section $i has [llength [lindex $residues $i]] residues"
   }
   package require psfgen
   set section 1
   mol new $psf
   mol addfile $pdb waitfor all
   set tmpmolid [molinfo top get id]
   foreach x $residues {
      resetpsf
      readpsf $psf pdb $pdb
      set keepsel [atomselect $tmpmolid "residue $x"]
      set delsel [atomselect $tmpmolid "not residue $x"]
      catch {
      foreach seg [$delsel get segname] resid [$delsel get resid]  {
         delatom $seg $resid
      }} rv
      writepsf "${prefix}${section}.psf"
      writepdb "${prefix}${section}.pdb"
      vmdcon -info "split_psf: wrote ${prefix}${section}.psf and ${prefix}${section}.pdb"
      mol new "${prefix}${section}.psf"
      mol addfile "${prefix}${section}.pdb" waitfor all
      set chkmolid [molinfo top get id]
      set sel [atomselect $chkmolid "all"]
      set nres [llength [lsort -unique [$sel get residue]]]
      $sel delete
      vmdcon -info "split_psf: section $section has $nres residues"
      if { $nres != [llength $x]} {
         vmdcon -err "split_psf: section $section assigned wrong number of residues: $nres != [llength $x]"
         exit
      }
      mol delete $chkmolid
      set section [expr $section + 1]
   }
   resetpsf
   mol delete $tmpmolid
}

proc PestiferUtil::verify_no_mols {} {
   catch {mol list} rv 
   if { $rv == "ERROR) No molecules loaded." } {
      return 1
   }
   return 0
}

proc PestiferUtil::shuffle_list {list} {
    set n [llength $list]
    for {set i [expr {$n - 1}]} {$i > 0} {incr i -1} {
        set j [expr {int(rand() * ($i + 1))}]
        # Swap elements i and j
        set temp [lindex $list $i]
        lset list $i [lindex $list $j]
        lset list $j $temp
    }
    return $list
}

proc PestiferUtil::ess { count } {
   # returns the number of elements in a list
   # if count is zero, returns 1
   if { $count > 1 } {
      return "s"
   } else {
      return ""
   }
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
proc PestiferUtil::mysign { x } {
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

proc PestiferUtil::letter_up { c {cyclic 0}} {
   if { $c == "Z" } {
      return "a"
   }
   if { $c == "z" && $cyclic == 1 } {
      return "A"
   }
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

proc PestiferUtil::backup { atomsel attributes } {
    set data [list]
    foreach attr $attributes {
        lappend data [$atomsel get $attr]
    }
    return $data
}

proc PestiferUtil::restore { atomsel attributes data } {
    foreach attr $attributes values $data {
        $atomsel set $attr $values
    }
}

proc PestiferUtil::namdbin2pdb { psf coor pdb } {
   mol new $psf
   mol addfile $coor waitfor all
   set a [atomselect top all]
   $a writepdb $pdb
}

proc PestiferUtil::pdb2namdbin { pdb coor } {
   mol new $pdb waitfor all autobonds off
   set a [atomselect top all]
   $a writenamdbin $coor
}

proc PestiferUtil::log_addframe { molid logid } {
   if { $logid != "-1" } {
     animate dup $logid
     [atomselect $logid all] set x [[atomselect $molid all] get x]
     [atomselect $logid all] set y [[atomselect $molid all] get y]
     [atomselect $logid all] set z [[atomselect $molid all] get z]
#     [atomselect $molid all] writepdb "tmp.pdb"
#     animate read pdb tmp.pdb $logid
     vmdcon -info "Molid $molid - logging molecule $logid has [molinfo $logid get numframes] frames."
#     exec rm -f tmp.pdb
   }
}

proc PestiferUtil::measure_bonds { psf coor infile fixedpdb outfile rfzr } {
    mol new $psf
    mol addfile $coor
    set a [atomselect top "all"]
    $a set occupancy 0
    set ca [atomselect top "name CA or (water) or (glycan and element C O)"]
    $ca set occupancy 1
    set fp [open $infile "r"]
    set lines [split [read $fp] \n]
    close $fp
    set RES0 {}
    set RES1 {}
    set RES2 {}
    set SEG {}
    foreach l $lines {
        if { [llength $l] == 4 } {
        lappend SEG [lindex $l 0]
        lappend RES0 [lindex $l 1]
        lappend RES1 [lindex $l 2]
        lappend RES2 [lindex $l 3]
        }
    }
    vmdcon -info "SEG $SEG"
    vmdcon -info "RES0 $RES0"
    vmdcon -info "RES1 $RES1"
    vmdcon -info "RES2 $RES2"

    set bl {}
    set CC {}
    set NN {}
    foreach seg $SEG h $RES0 i $RES1 j $RES2 {
        set labelus [atomselect top "segname $seg and resid $h to $i"]
        $labelus set occupancy 0
        set ci [string index $i end]
        if {[string is alpha $ci]} {
            set resid [string range $i 0 end-1]
            set theC [atomselect top "segname $seg and resid $resid and insertion $ci and name C"]
        } else {
            set theC [atomselect top "segname $seg and resid $i and name C"]
        }
        set ni [string index $j end]
        if {[string is alpha $ni]} {
            set resid [string range $j 0 end-1]
            set theN [atomselect top "segname $seg and resid $resid and insertion $ni and name N"]
        } else {
            set theN [atomselect top "segname $seg and resid $j and name N"]
        }
        if {[expr $rfzr > 0]} {
            set theN_idx [$theN get index]
            set flex_zone [atomselect top "segname $seg and name CA and same residue as (exwithin $rfzr of index $theN_idx)"]
            $flex_zone set occupancy 0
        }
        vmdcon -info "resid $i on segname $seg has [$theC num] Cs with serial [$theC get serial]"
        vmdcon -info "resid $j on segname $seg has [$theN num] Ns with serial [$theN get serial]"
        set ii [$theC get index]
        set jj [$theN get index]
        lappend bl [measure bond [list $ii $jj]]
        lappend CC [$theC get serial]
        lappend NN [$theN get serial]
    }
    vmdcon -info "CC $CC"
    vmdcon -info "NN $NN"

    set fp [open $outfile "w"]
    puts $fp "# segname c-serial n-serial distance(A)"
    foreach s $SEG i $CC j $NN b $bl {
        puts $fp "$s $i $j [format %.4f $b]"
    }
    close $fp
    $a writepdb $fixedpdb
    vmdcon -info "measure_bonds: results written to $outfile"
    vmdcon -info "   and fixed pdb written to $fixedpdb"
}

proc PestiferUtil::MaxRad { atomsel } {
   set com [measure center $atomsel]
   set cx [lindex $com 0]
   set cy [lindex $com 1]
   set maxrad -99.99
   foreach x [$atomsel get x] y [$atomsel get y] {
      set dx [expr $x - $cx]
      set dy [expr $y - $cy]
      set thisr2 [expr $dx*$dx + $dy*$dy]
      if { $thisr2 > $maxrad } {
         set maxrad $thisr2
      }
   }
   return [expr sqrt($maxrad)]
}

proc PestiferUtil::IsIntersectionEmpty { list1 list2 } {

   foreach elem $list1 {
      if {[lsearch -exact $list2 $elem] != -1} {
         return 0
      }
   }
   return 1

}

proc PestiferUtil::deprotect_str_arg { arg } {
# replaces underscores with spaces
   set arg [string map { _ " " } $arg]
   return $arg
}