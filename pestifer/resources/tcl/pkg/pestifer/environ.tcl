# Author: Cameron F. Abrams <cfa22@drexel.edu>
package provide PestiferEnviron 1.0

namespace eval ::PestiferEnviron:: {
    namespace export *
}

proc PestiferEnviron::leaflet_apportionment { molid } {
   set whole [atomselect $molid all]
   set bilayer [atomselect $molid "lipid"]
   set com [measure center $bilayer weight mass]
   set com_z [lindex $com 2]
   $whole moveby [list 0 0 -$com_z]
   set residue_list [lsort -unique [atomselect $molid "residue"]]
   set residues(upper) [list]
   set residues(lower) [list]
   foreach residue $residue_list {
      set ressel [atomselect $molid "residue $residue"]
      set com [measure center $ressel weight mass]
      set com_z [lindex $com 2]
      if { $com_z > 0 } {
         lappend residues(upper) $residue
      } else {
         lappend residues(lower) $residue
      }
   }
   set residues(upper) [lsort -unique $residues(upper)]
   set residues(lower) [lsort -unique $residues(lower)]
   return $residues
}

proc PestiferEnviron::write_psfgen { molid {next_available_chain A} {segtypes {lipid water ion}} {Slet {L I W}} {maxr_per_seg 1000}} {
    # execute segment and coordpdb stanzas for all atoms in the molid
    foreach segtype $segtypes S $Slet {
        set a [atomselect $molid "$segtype"]
        if { [$a num] > 0 } {
            $a set chain $next_available_chain

            set ridx [lsort -integer -unique [$a get residue]]
            set nres [llength $ridx]
            set nseg [expr $nres / $maxr_per_seg + 1]
            set mm [expr $nres % $nseg]
            vmdcon -info "[$a num] $segtype atoms (chain $next_available_chain) in $nres residues divide into $nseg segments"

            for { set seg 1 } { $seg <= $nseg } { incr seg } {
                set segname "${next_available_chain}${seg}"
                set left [expr (${seg}-1)*$maxr_per_seg]
                set right [expr $left + $maxr_per_seg - 1]
                if { $right >= [llength $ridx] } {
                    set right end
                }
                vmdcon -info "Segment $segname is residues [lindex $ridx $left] to [lindex $ridx $right]"
                set res [lrange $ridx $left $right]
                set segsel [atomselect $molid "$segtype and residue $res"]
                set sridx [lsort -integer -unique [$segsel get residue]]
                set rser 1
                foreach x $sridx {
                    set sridx_sermap($x) $rser
                    incr rser
                }
                set rser [list]
                foreach x [$segsel get residue] {
                    lappend rser $sridx_sermap($x)
                }
                $segsel set resid $rser
                vmdcon -info "Segment $segname has [$segsel num] atoms"
                $segsel writepdb "${segname}_tmp.pdb"
                segment $segname {
                    auto none
                    first none
                    last none
                    pdb ${segname}_tmp.pdb
                }
                coordpdb ${segname}_tmp.pdb $segname
                $segsel delete
            }
            set next_available_chain [letter_up $next_available_chain]
        }
    }

}