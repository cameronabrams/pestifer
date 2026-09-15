# Solvated-system figure for an example page, in the style of the BPTI series (13457d00) with the
# orange "what this example is about" rep of 0d6e72a9.  The script behind those figures was never
# committed; this one was rebuilt on 2026-09-15 and checked against 6pti-solvated.png (same water
# line color in the darkest decile, 149,178,178).
#
# Feed on STDIN, never -e (see scripts/pestifer-snapshot for why):
#   vmd -dispdev text -args PSF PDB OUT.tga "ORANGE SEL" ["AXES"] ["EXTRA CMD"] < render_solvated_figure.tcl
# then crop to the non-white bounding box and scale to 855 px wide, like the existing figures.
#   AXES       optional list of view directions to choose among, e.g. "{0 1 0} {0 -1 0}"
#   EXTRA CMD  optional VMD command applied after the view, e.g. "rotate z by 90"
# Used for the four figures added 2026-09-15 (examples/28-31), from the 3.22.1 sweep builds:
#   28 "resname SEP or (protein and resid 57 and resname SER and not backbone and not name HN HA)"
#   29 "resname GLYM and not backbone and not name HN HA1 HA2"
#   30 "resname MG"
#   31 "resname AFUC BGLC BXYL" "{0 1 0} {0 -1 0}" "rotate z by 90"
proc main {} {
    global argv
    lassign $argv psf pdb out hl axes
    # optional 5th argument: the view axes to choose among (default all six).  An elongated
    # system shot down its long axis is mostly water; restrict it to the short ones.
    display projection Orthographic
    display depthcue on
    display cuemode Exp2
    display cuestart 0.5
    display cueend 10.0
    display cuedensity 0.32
    display shadows off
    display ambientocclusion off
    axes location off
    color Display Background white
    display resize 1800 1800

    set m [mol new $psf type psf waitfor all]
    mol addfile $pdb type pdb waitfor all molid $m
    mol delrep 0 $m
    # water: the BPTI-series look, matched against 6pti-solvated.png (darkest-decile line color
    # 149,178,178 in both).  VMD 2.0's Glass1 renders paler than the reference, so a copy with more
    # opacity and a muted teal on a color id that no Structure color uses (cyan is turns)
    mol representation Lines 4.0
    color change rgb 17 0.40 0.56 0.56
    material add WaterGlass copy Glass1
    material change opacity WaterGlass 0.30
    mol color ColorID 17
    mol selection {water or ion}
    mol material WaterGlass
    mol addrep $m
    mol representation CPK 1.0 0.3 12.0 12.0
    mol color Name
    mol selection {protein}
    mol material Opaque
    mol addrep $m
    mol representation NewCartoon 0.3 10.0 4.1 0
    mol color Structure
    mol selection {protein}
    mol material Opaque
    mol addrep $m
    mol representation CPK 1.5 0.5 12.0 12.0
    mol color ColorID 3
    mol selection "($hl) and not hydrogen"
    mol material Opaque
    mol addrep $m

    # choose the axis-aligned view that puts the highlight nearest the camera, so the box stays
    # square to the frame (a tilted box frays at the corners)
    set p [measure center [atomselect $m "protein"]]
    set h [measure center [atomselect $m "$hl"]]
    set d [vecsub $h $p]
    # camera looks down -z of the rotated frame; each entry: rotation cmds, world direction toward camera
    set views {
        {{} {0 0 1}}
        {{rotate y by 180} {0 0 -1}}
        {{rotate x by -90} {0 -1 0}}
        {{rotate x by 90} {0 1 0}}
        {{rotate y by -90} {1 0 0}}
        {{rotate y by 90} {-1 0 0}}
    }
    set best ""; set bestv -1e9
    if {$axes ne ""} {
        set keep {}
        foreach v $views { if {[lsearch -exact $axes [lindex $v 1]] >= 0} { lappend keep $v } }
        set views $keep
    }
    foreach v $views {
        set s [vecdot $d [lindex $v 1]]
        if {$s > $bestv} { set bestv $s; set best [lindex $v 0] }
    }
    display resetview
    if {$best ne ""} { eval $best }
    if {[llength $argv] > 5} { eval [lindex $argv 5] }
    scale by 1.9
    display update
    render TachyonInternal $out
    puts "RENDERED $out view={$best} front=[format %.1f $bestv]"
}
main
quit
