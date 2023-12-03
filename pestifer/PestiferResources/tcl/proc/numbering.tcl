# Author: Cameron F. Abrams, <cfa22@drexel.edu>

proc renumber_as_insertion { atomsel root_resid root_insertion } {
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