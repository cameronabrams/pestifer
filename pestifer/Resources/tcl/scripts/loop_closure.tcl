for { set i 0 } { $i < [llength $argv] } { incr i } {
    if { [lindex $argv $i] == "-pdb"} {
       incr i
       set coor [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-psf"} {
       incr i
       set psf [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-o"} {
       incr i
       set outbasename [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-p"} {
       incr i
       set patchfile [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-files"} {
       incr i
       set filelistfile [lindex $argv $i]
    }
}

set outpsf ${outbasename}.psf
set outpdb ${outbasename}.pdb

mol new $psf
mol addfile $coor
set m0 [molinfo top get id]

set segnames [lsort -unique [[atomselect $m0 "all"] get segname]]

set pdbnames [list]

foreach s $segnames {
    [atomselect top "segname $s"] writepdb SEG$s.pdb
    segment $s {
        pdb SEG$s.pdb
    }
    coordpdb SEG$s.pdb $s
    lappend pdbnames SEG$s.pdb
}

source $patchfile

regenerate angles dihedrals
guesscoord

writepsf $outpsf
writepdb $outpdb
set fp [open $filelistfile "w"]
foreach fn $pdbnames {
    puts $fp "$fn"
}
close $fp
exit
