# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Aligns MPER-TM trimer (PDB 68EW) to C-termini of HIV-1 Env SOSIP open trimer (PDB 5VN3)
# in preparation to fuse to generate a trimer model with gp41's that go all the way to 710
#

# proc bump { molid sel } {
#     $sel writepdb tmp.pdb
#     mol addfile tmp.pdb waitfor all molid $molid
#     $sel frame last
#     file delete tmp.pdb
# }

proc sign { x } {
    if { $x < 0 } {
        return -1
    } else {
        return 1
    }
}

mol new 5vn3
set envid [molinfo top get id]
mol new 6e8w waitfor all
set mtmid [molinfo top get id]
mol top $envid

# keep only the last model in the series of NMR models
set nf [molinfo $mtmid get numframes]
animate delete beg 0 end [expr $nf - 2] $mtmid

set tri [atomselect $envid all]
set trictr [measure center $tri]
set mtm [atomselect $mtmid all]
set mtctr [measure center $mtm]

# orient 3-fold axis of MPER-TM trimer along z

set hd [atomselect $mtmid "resid < 664"]
set tl [atomselect $mtmid "resid > 690"]

set hp [measure center $hd]
set tp [measure center $tl]

set axv [vecsub $tp $hp]
set R [veclength $axv]
set axuv [vecscale $axv [expr 1.0/$R]]

set pi [expr 2*asin(1.0)]
set dPhi [expr -[sign [lindex $axuv 1]]*180.0/$pi*acos([lindex $axuv 0])]
set dTheta [expr -(90.0-180.0/$pi*asin([lindex $axuv 2]))]
if { $dTheta < 0.0 } { 
    set dTheta [expr $dTheta + 180.0]
}
puts "$pi $dPhi $dTheta"

$mtm move [trans center [measure center $mtm] axis z $dPhi]
# bump $mtmid $mtm
$mtm move [trans center [measure center $mtm] axis y $dTheta]
# bump $mtmid $mtm


# set cartesian positions of MPER-TM
set cshift [vecsub $trictr $mtctr]
$mtm moveby [list [lindex $cshift 0] [lindex $cshift 1] 63.0]
# bump $mtmid $mtm
rotate to align chains A-A, B-B, and D-C
$mtm move [trans center [measure center $mtm] axis z 132]
bump $mtmid $mtm
we could use a little fine-tuning
$mtm move [trans center [measure center $mtm] axis y -3]
# bump $mtmid $mtm
$mtm move [trans center [measure center $mtm] axis z 11]
# bump $mtmid $mtm
$mtm move [trans center [measure center $mtm] axis y -1]
# bump $mtmid $mtm
$mtm moveby {2 0 0}
# bump $mtmid $mtm
$mtm move [trans center [measure center $mtm] axis y -2]
# bump $mtmid $mtm

foreach ec {A B D} mc {A B C} {
    set h1r [[atomselect $mtmid "chain $mc and resid 673 and name CA"] get residue]
    set h2r [expr $h1r-13]
    brot $mtmid $h1r $h2r psi N 30
    # bump $mtmid $mtm
    set j1r [[atomselect $envid "chain $ec and resid 635 and name CA"] get residue]
    set j2r [[atomselect $envid "chain $ec and resid 664 and name CA"] get residue]
    brot $envid $j1r $j2r phi C 10
}
brot 1 115 107 phi N -14
# bump $mtmid $mtm
brot 1 115 107 psi N -13
# bump $mtmid $mtm

# animate write dcd "6e8w-star.dcd" sel $mtm $mtmid

$tri writepdb "5vn3-star.pdb"
$mtm writepdb "6e8w-star.pdb"
exit
