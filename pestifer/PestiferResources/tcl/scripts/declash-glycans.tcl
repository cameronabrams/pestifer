# Author: Cameron F. Abrams, <cfa22@drexel.edu>

# A VMD/TcL script to adjust all glycan dihedrals to eliminate steric clashes
set clashdist 2.0
set maxcycles 100
set logfilename "declash-glycans.log"
for { set i 0 } { $i < [llength $argv] } { incr i } {
    if { [lindex $argv $i] == "-pdb"} {
       incr i
       set pdb [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-psf"} {
       incr i
       set psf [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-d"} {
        incr i
        set datafile [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-o"} {
        incr i
        set outpdb [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-log"} {
        incr i
        set logfilename [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-clashdist"} {
        incr i
        set clashdist [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-maxcycles"} {
        incr i
        set maxcycles [lindex $argv $i]
    }
}

mol new $psf
mol addfile $pdb waitfor all
set a [atomselect top all]
set molid [molinfo top get id]
source $datafile
vmdcon -info "Declashing $nglycans glycans; clashdist $clashdist; maxcycles $maxcycles"
set logf [open $logfilename "w"]
for {set i 0} {$i<$nglycans} {incr i} {
    vmdcon -info "Glycan $i has [llength $glycan_idx($i)] atoms and [llength $rbonds($i)] rotatable bonds"
    declash_pendant $molid $glycan_idx($i) $rbonds($i) $movers($i) $maxcycles $clashdist $logf
}

$a writepdb ${outpdb}
close $logf
exit
