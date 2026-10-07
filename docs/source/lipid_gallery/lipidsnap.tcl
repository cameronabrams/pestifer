# Render one lipid conformer: heavy atoms, coloured by element, phosphorus marked.
# MUST be fed to VMD on stdin -- see scripts/pestifer_snapshot.tcl for why `-e` renders
# an empty scene.
set psf ""; set pdb ""; set out "conf.png"; set W 440; set H 760; set SCALE 0.030
for {set i 0} {$i < [llength $argv]} {incr i} {
    switch -- [lindex $argv $i] {
        -psf  {set psf [lindex $argv [incr i]]}
        -coor {set pdb [lindex $argv [incr i]]}
        -o    {set out [lindex $argv [incr i]]}
        -size {regexp {(\d+)x(\d+)} [lindex $argv [incr i]] -> W H}
        -scale {set SCALE [lindex $argv [incr i]]}
    }
}
mol new $psf type psf waitfor all
mol addfile $pdb type pdb waitfor all
mol delrep 0 top

proc addrep {sel style color {mat Opaque}} {
    mol representation $style
    mol color $color
    mol selection $sel
    mol material $mat
    mol addrep top
}
addrep "noh"                 "Licorice 0.28 24.0 24.0" "Name"
addrep "name P and not hydrogen" "VDW 0.55 24.0"       "Name"

display projection Orthographic
display depthcue off
display shadows off
axes location Off
color Display Background white
color Name C silver
light 0 on; light 1 on

display resetview
# Face the membrane normal: conformers are built with z along it, and VMD looks down -z.
rotate x by -90
rotate y by 15
# FIXED scale, not resetview's per-molecule fit: every conformer of a gallery entry must share
# one scale, or a folded conformer is silently zoomed up to fill the frame and reads as normal.
scale to $SCALE
display resize $W $H
display update
set tga "[file rootname $out].tga"
render TachyonInternal $tga
if {[file exists $tga]} {
    exec convert $tga $out
    file delete $tga
    puts "lipidsnap: wrote $out"
}
quit
