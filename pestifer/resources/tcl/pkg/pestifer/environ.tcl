# Author: Cameron F. Abrams <cfa22@drexel.edu>
package provide PestiferEnviron 1.0

namespace eval ::PestiferEnviron:: {
    namespace export *
}

proc PestiferEnviron::leaflet_apportionment { molid } {
    vmdcon -info "Apportioning residues to upper and lower leaflets"
    set whole [atomselect $molid all]
    set bilayer [atomselect $molid "lipid"]
    set bilayer_com [measure center $bilayer weight mass]
    set bilayer_com_z [lindex $bilayer_com 2]
    vmdcon -info "Bilayer center of mass z coordinate: $bilayer_com_z"
    set residue_list [lsort -unique [$bilayer get residue]]
    vmdcon -info "[llength $residue_list] residues in the bilayer"
    set residues_upper [list]
    set residues_lower [list]
    foreach residue $residue_list {
        set ressel [atomselect $molid "residue $residue"]
        set residue_com [measure center $ressel weight mass]
        set residue_com_z [lindex $residue_com 2]
        if { $residue_com_z > $bilayer_com_z } {
            lappend residues_upper $residue
        } else {
            lappend residues_lower $residue
        }
        $ressel delete
    }
    set upperchamber [atomselect $molid "(water or ion) and z > $bilayer_com_z"]
    set urix [lsort -unique [$upperchamber get residue]]
    foreach residue $urix {
        lappend residues_upper $residue
    }
    set lowerchamber [atomselect $molid "(water or ion) and z < $bilayer_com_z"]
    set lrix [lsort -unique [$lowerchamber get residue]]
    foreach residue $lrix {
        lappend residues_lower $residue
    }
    $whole delete
    $bilayer delete
    $upperchamber delete
    $lowerchamber delete
    return [list $residues_upper $residues_lower]
}

proc PestiferEnviron::write_psfgen { molid {segtypes {lipid ion water}} 
                                    {seglabels {L I WT}} {segidx {1 1 1}} 
                                    {maxr_per_seg 9999} {sac_chain A}} {
    # execute segment and coordpdb stanzas for all atoms in the molid
    set new_segidx [list]
    vmdcon -info "Writing psfgen segments for molid $molid with chain $sac_chain, 
                 segment types $segtypes, labels $seglabels, segidx $segidx, and max residues per segment $maxr_per_seg"
    foreach segtype $segtypes seglabel $seglabels idx $segidx {
        vmdcon -info "Processing segment type $segtype with label $seglabel and index $idx"
        set a [atomselect $molid "$segtype"]
        if { [$a num] > 0 } {
            $a set chain $sac_chain
            set ridx [lsort -integer -unique [$a get residue]]
            set nres [llength $ridx]
            set nseg [expr $nres / $maxr_per_seg + 1]
            set mm [expr $nres % $nseg]
            vmdcon -info "[$a num] $segtype atoms $nres residues divide into $nseg segment[ess $nseg]"
            if { $mm > 0 } {
                vmdcon -info "the last of which has $mm residues"
            }
            for { set seg 1 } { $seg <= $nseg } { incr seg } {
                set this_segidx [expr $seg + $idx - 1]
                set segname "${seglabel}${this_segidx}"
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
            lappend new_segidx [expr $seg + $idx]
        } else {
            lappend new_segidx $idx
        }
    }
    vmdcon -info "All segments processed, new segment indices: $new_segidx"
    return $new_segidx
}