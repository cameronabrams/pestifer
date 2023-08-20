
set psf [lindex $argv 0]
set coor [lindex $argv 1]
set patchfile [lindex $argv 2]
set outpsf [lindex $argv 3]
set outpdb [lindex $argv 4]

mol new $psf
mol addfile $coor
set m0 [molinfo top get id]

set segnames [lsort -unique [[atomselect $m0 "all"] get segname]]

foreach s $segnames {
    [atomselect top "segname $s"] writepdb SEG$s.pdb
    segment $s {
        pdb SEG$s.pdb
    }
    coordpdb SEG$s.pdb $s
}

source patches.inp

#### LIGATION LIST STARTS
#### LIGATION LIST ENDS

regenerate angles dihedrals
guesscoord

writepsf $outpsf
writepdb $outpdb
exit
