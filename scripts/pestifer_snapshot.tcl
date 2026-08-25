# pestifer_snapshot.tcl -- headless VMD snapshots of a pestifer-built system.
#
#   vmd -dispdev text -args -psf X.psf -coor X.pdb -o out.png [options] < pestifer_snapshot.tcl
#
# NOTE THE REDIRECT.  The script must be fed to VMD on **stdin**, not with `-e`.  A `-e` script
# is executed before VMD's event loop, so the drawable geometry behind a representation is never
# generated: every render then silently writes an image containing nothing but the corner axes,
# and a Tachyon scene export is byte-identical with and without a molecule loaded.  Delivered on
# stdin the commands run inside the event loop and the geometry is there.  (Use the
# pestifer-snapshot wrapper, which gets this right.)
#
# It is specifically the *geometry*, not the molecule: under `-e`, numatoms, numreps, numframes
# and every atomselect count are correct, and STRIDE has already assigned secondary structure.
# So `waitfor all` is not the missing piece -- it makes the file load synchronous (which still
# matters for a .dcd, and is used below), not the representation drawable.  Nor is there an
# in-script escape hatch: `display update`, `display update ui`, `display update on`, Tcl's
# `update` and `update idletasks`, `after`, and `mol showrep` were all measured under `-e` and
# leave the scene at axes-only.  Hence the delivery mechanism is the fix.
#
# Options (only -psf and -coor are required):
#   -psf FILE        topology
#   -coor FILE       coordinates: .pdb, .coor (NAMD binary) or .dcd
#   -frame N         frame to render from a .dcd (default: last)
#   -o FILE          output image; .png .tga .bmp .ppm (default: snapshot.png)
#   -size WxH        image size in pixels (default: 1600x1200)
#   -view NAME       xy | xz | yz | auto  (default: auto -- xz for membranes, xy otherwise)
#   -style NAME      auto | protein | membrane | glycan  (default: auto)
#   -solvent MODE    off|0 (default), or how to draw water/bulk solvent: points|lines|licorice
#                    (1 == points).  Hydrogens are never drawn on solvent.  Ions come with it.
#   -highlight SELS  ';'-separated VMD selections drawn as thick licorice in contrasting
#                    colors on top of everything else, e.g. "resid 11 13; resid 5 55".
#                    Hydrogens are dropped unless a selection matches only hydrogens.
#   -ss 0|1          draw cysteine sidechains (no H) as one licorice rep (default: 0), so a
#                    disulfide reads as a bond rather than as two nearby spheres
#   -face SEL        rotate so this selection faces the camera, i.e. sits in the foreground
#                    rather than behind the bulk of the structure
#   -ghost 0|1       draw the bulk cartoon translucent (default: 0), so a feature enclosed by
#                    the structure -- DNA threaded through a polymerase, a buried cofactor --
#                    is visible through it.  Highlights stay opaque.
#   -side SEL        rotate so this selection sits to one side, in the image plane -- what you
#                    want for a two-domain construct, where -face would stack them front-to-back
#   -domains SELS    ';'-separated selections, each drawn as its own cartoon in a distinct color
#                    (replaces the single protein cartoon).  For fusions and multi-domain builds,
#                    where one color per chain says nothing because it is all one chain.
#   -bg NAME         background color (default: white)
#   -renderer NAME   TachyonInternal | TachyonLOptiXInternal | snapshot (default: TachyonInternal)
#   -fill F          fraction of the frame the geometry should fill (default: 0.92)
#   -zoom F          extra scale factor applied after the auto-fit (default: 1.0)
#   -rot "X Y Z"     extra rotation in degrees applied after -view
#
# Everything runs inside procs on purpose: VMD feeds an -e file to its Tcl console a line at a
# time, which echoes every top-level command's result and mishandles some multi-line blocks.
# Inside a proc the body is parsed as one unit and nothing is echoed.
#
# Exits nonzero on any error so a calling pipeline can detect failure.

# CHARMM carbohydrate resnames psfgen writes (not the PDB names pestifer aliases away)
set PS_GLYCAN_RESNAMES {AGLC BGLC AGAL BGAL AMAN BMAN AFUC BFUC AXYL BXYL ANE5AC
                        BGLCNA AGLCNA BGALNA AGALNA ARHM IDOA BGLCA}
# `lipid` is a VMD macro but misses the CHARMM36 names pestifer builds with
set PS_LIPID_EXTRA {CHL1 CHOL POPC POPE POPS POPI POPA PLPC DMPC DPPC DOPC DLPC DSPC
                    PSM SSM DPSM OSM SAPI TOCL1 TOCL2 LSM}
# Counterions added by solvation: they belong with the solvent and are hidden with it.
set PS_ION_RESNAMES {SOD CLA POT CES RUB}
# Structural/bound metals: always drawn.  A zinc in an active site or a magnesium on a
# nucleotide is part of the structure, not part of the bath, and hiding it with the solvent
# silently drops the thing an active-site figure is usually about.
set PS_METAL_RESNAMES {ZN ZN2 MG MGA CAL CA MN MN2 FE FE2 CU CU1 NI CO CD HG K2}
# A species present in at least this many residue copies is bulk solvent, not a ligand.  Real
# separation is wide -- a non-aqueous solvent box runs to hundreds or thousands of copies while
# a cofactor comes in ones -- so the exact value is not delicate.
set PS_BULK_MIN_COPIES 50
# ColorIDs for -highlight, picked to contrast with the NewCartoon Structure palette
# (purple/yellow/cyan/white/blue) so a highlighted residue never blends into the ribbon
set PS_HIGHLIGHT_COLORS {1 3 7 9 5}
# Cartoon colors for -domains; well separated from each other at a glance
set PS_DOMAIN_COLORS {0 1 7 3 11 9}

proc ps_magick {} {
    foreach cand {magick convert} {
        if {![catch {exec which $cand} p]} { return $p }
    }
    return ""
}

proc ps_die {msg} {
    puts stderr "pestifer_snapshot: ERROR: $msg"
    exit 1
}

proc ps_note {msg} {
    puts "pestifer_snapshot: $msg"
}

# number of atoms matching a selection; 0 if the selection text is invalid
proc ps_nsel {molid text} {
    if {[catch {atomselect $molid $text} s]} { return 0 }
    set n [$s num]
    $s delete
    return $n
}

proc ps_parse_args {argv} {
    array set opt {
        psf "" coor "" frame -1 o "snapshot.png" size "1600x1200" view auto
        style auto solvent 0 bg white renderer TachyonInternal zoom 1.0 rot "" fill 0.92
        highlight "" ss 0 face "" ghost 0 side "" domains ""
    }
    for {set i 0} {$i < [llength $argv]} {incr i} {
        set a [lindex $argv $i]
        if {![string match "-*" $a]} continue
        set key [string range $a 1 end]
        # pestifer's VMD launcher injects flags of its own (--tcl-root ...); skip anything
        # we do not define rather than dying on it
        if {![info exists opt($key)]} continue
        incr i
        set opt($key) [lindex $argv $i]
    }
    return [array get opt]
}

proc ps_load {molid_var psf coor frame} {
    upvar 1 $molid_var molid
    set ext [string tolower [file extension $coor]]
    switch -- $ext {
        .pdb    { set coortype pdb }
        .coor   { set coortype namdbin }
        .dcd    { set coortype dcd }
        default { ps_die "unrecognized coordinate extension '$ext' (want .pdb, .coor or .dcd)" }
    }
    set molid [mol new $psf type psf waitfor all]
    if {$molid < 0} { ps_die "could not read $psf" }
    if {[catch {mol addfile $coor type $coortype waitfor all molid $molid} e]} {
        ps_die "could not read $coor: $e"
    }
    set nframes [molinfo $molid get numframes]
    if {$nframes == 0} { ps_die "no frames loaded from $coor" }
    if {$frame < 0 || $frame >= $nframes} { set frame [expr {$nframes - 1}] }
    # drop every other frame so `all` selections and the view fit refer to the rendered one
    for {set f [expr {$nframes - 1}]} {$f >= 0} {incr f -1} {
        if {$f != $frame} { animate delete beg $f end $f $molid }
    }
    ps_note "loaded [molinfo $molid get numatoms] atoms; rendering frame $frame of $nframes"
}

# Adds a representation and returns its index, or -1 if the selection is empty.
proc ps_addrep {molid seltext style color {material Opaque}} {
    if {[ps_nsel $molid $seltext] == 0} { return -1 }
    mol representation $style
    mol color $color
    mol selection $seltext
    mol material $material
    mol addrep $molid
    return [expr {[molinfo $molid get numreps] - 1}]
}

# Map a -solvent mode onto a representation.  "1" is accepted as a synonym for points so the
# old boolean spelling keeps working.
# Indices of the representations holding solvent, so they can be hidden while the view is fitted.
proc ps_solvent_reps {} {
    if {![info exists ::PS_SOLVENT_REPS]} { return {} }
    return $::PS_SOLVENT_REPS
}

proc ps_solvent_style {mode} {
    switch -- $mode {
        1 - points   { return "Points 3.000000" }
        lines        { return "Lines 1.000000" }
        licorice     { return "Licorice 0.100000 12.000000 12.000000" }
        default      { return "Points 3.000000" }
    }
}

# Normalize -solvent into "off" or a mode name, so callers can test it as a string.
proc ps_solvent_mode {v} {
    if {$v eq "" || $v eq "0" || $v eq "off" || $v eq "no" || $v eq "false"} { return "off" }
    return $v
}

# Resnames in *seltext* present in enough residue copies to count as bulk solvent.
proc ps_bulk_species {molid seltext} {
    if {[catch {atomselect $molid $seltext} s]} { return {} }
    if {[$s num] == 0} { $s delete ; return {} }
    array set resname_of {}
    foreach rec [$s get {segname resid resname}] {
        lassign $rec sg rid rn
        set resname_of($sg:$rid) $rn
    }
    $s delete
    array set copies {}
    foreach k [array names resname_of] { incr copies($resname_of($k)) }
    set bulk {}
    foreach rn [array names copies] {
        if {$copies($rn) >= $::PS_BULK_MIN_COPIES} { lappend bulk $rn }
    }
    return [lsort $bulk]
}

proc ps_represent {molid style show_solvent highlight show_ss ghost domains} {
    global PS_GLYCAN_RESNAMES PS_LIPID_EXTRA PS_ION_RESNAMES PS_METAL_RESNAMES PS_HIGHLIGHT_COLORS PS_DOMAIN_COLORS
    set glycan_sel "resname $PS_GLYCAN_RESNAMES"
    set lipid_sel  "lipid or resname $PS_LIPID_EXTRA"
    set ion_sel    "resname $PS_ION_RESNAMES"
    set metal_sel  "resname $PS_METAL_RESNAMES"

    mol delrep 0 $molid

    # Color by Chain only when there is more than one chain to tell apart; a single-chain
    # system rendered "by chain" is one flat color, where secondary structure is informative.
    set protcolor Structure
    if {[ps_nsel $molid protein] > 0} {
        set s [atomselect $molid "protein"]
        set nchain [llength [lsort -unique [$s get chain]]]
        $s delete
        if {$nchain > 1} { set protcolor Chain }
    }
    set bulk_material [expr {$ghost ? "Transparent" : "Opaque"}]
    # -domains replaces the single protein cartoon with one cartoon per named selection, each in
    # its own color.  Anything not claimed by a domain still gets a cartoon, so no residue is
    # dropped just because the domain list is incomplete.
    set domain_sels {}
    foreach d [split $domains ";"] {
        set d [string trim $d]
        if {$d ne ""} { lappend domain_sels $d }
    }
    if {[llength $domain_sels] > 0} {
        set di 0
        foreach d $domain_sels {
            set n [ps_nsel $molid $d]
            if {$n == 0} { ps_note "WARNING: -domains selection matched no atoms: $d" ; continue }
            set c [lindex $PS_DOMAIN_COLORS [expr {$di % [llength $PS_DOMAIN_COLORS]}]]
            ps_addrep $molid $d "NewCartoon 0.300000 20.000000 3.000000 0" "ColorID $c" $bulk_material
            ps_note "domain [expr {$di + 1}] (ColorID $c): $n atoms -- $d"
            incr di
        }
        set rest "protein and not ([join $domain_sels ") and not ("])"
        ps_addrep $molid $rest "NewCartoon 0.300000 20.000000 3.000000 0" $protcolor $bulk_material
    } else {
        ps_addrep $molid "protein" "NewCartoon 0.300000 20.000000 3.000000 0" $protcolor $bulk_material
    }
    ps_addrep $molid "nucleic" "NewCartoon 0.300000 20.000000 3.000000 0" Chain $bulk_material
    ps_addrep $molid $glycan_sel "Licorice 0.250000 20.000000 20.000000" Name
    ps_addrep $molid $metal_sel "VDW 0.600000 20.000000" ResName

    if {$style eq "membrane"} {
        # Thin licorice keeps the tails from swamping the protein; the phosphorus beads mark the
        # head-group planes, which is what makes a bilayer side view readable at all.  Colour by
        # ResName only when the leaflet is actually a mixture -- on a single-species bilayer that
        # collapses to one flat colour, and the beads may as well be a contrasting one.
        set s [atomselect $molid $lipid_sel]
        set nspecies [llength [lsort -unique [$s get resname]]]
        $s delete
        set beadcolor [expr {$nspecies > 1 ? "ResName" : "ColorID 3"}]
        ps_addrep $molid $lipid_sel "Licorice 0.100000 12.000000 12.000000" \
            [expr {$nspecies > 1 ? "ResName" : "ColorID 6"}]
        ps_addrep $molid "name P and ($lipid_sel)" "VDW 0.700000 16.000000" $beadcolor
    }
    # Ions travel with the solvent: counterions scattered through a water box read as stray dots
    # around the solute and, worse, stretch the view fit to the whole box.  They come back with it.
    # Solvent is drawn without hydrogens -- at box scale they add a haze of geometry that obscures
    # the solute and costs render time without conveying anything.
    if {$show_solvent ne "off"} {
        set wstyle [ps_solvent_style $show_solvent]
        set ::PS_SOLVENT_REPS {}
        # neutral gray: solvent is context, and coloring it by element puts a field of red
        # oxygens in competition with the highlights
        foreach r [list [ps_addrep $molid "water and noh" $wstyle "ColorID 6"] \
                        [ps_addrep $molid $ion_sel "VDW 0.500000 16.000000" ResName]] {
            if {$r >= 0} { lappend ::PS_SOLVENT_REPS $r }
        }
    }

    # catch-all: ligands, cofactors, anything unrecognized, so a snapshot never silently
    # omits part of the system
    set covered "protein or nucleic or water or ($ion_sel) or ($metal_sel) or ($lipid_sel) or ($glycan_sel)"
    # Whatever remains is either bulk solvent or a ligand, and they want opposite treatment.  A
    # non-aqueous solvent (DMSO, acetone, acetonitrile) is not matched by VMD's `water` macro, so
    # without this it falls to the catch-all below and buries the solute under thousands of
    # licorice molecules -- and drags the view fit out to the whole box.  Telling them apart by
    # residue copy count needs no list of solvent names, so it holds for any solvent.
    set bulk_sel ""
    set bulk [ps_bulk_species $molid "not ($covered)"]
    if {[llength $bulk] > 0} {
        set bulk_sel "resname $bulk"
        set covered "$covered or ($bulk_sel)"
        if {$show_solvent ne "off"} {
            # Element colors here, unlike water: after dropping hydrogens water is a field of
            # identical oxygens, but a non-aqueous solvent has chemistry worth seeing -- and it
            # is what distinguishes otherwise identical builds (DMSO vs acetone vs acetonitrile).
            set r [ps_addrep $molid "($bulk_sel) and noh" [ps_solvent_style $show_solvent] Name]
            if {$r >= 0} { lappend ::PS_SOLVENT_REPS $r }
        } else {
            ps_note "hiding bulk solvent ($bulk); pass -solvent points to show it"
        }
    }
    if {[ps_addrep $molid "not ($covered)" "Licorice 0.300000 20.000000 20.000000" Name] >= 0} {
        ps_note "note: [ps_nsel $molid "not ($covered)"] atoms fell to the catch-all representation"
    }

    # Cysteine sidechains as ONE licorice rep: drawing only the SG atoms (as spheres) shows
    # where the sulfurs are but not whether they are bonded, which is the thing a reader is
    # actually looking for.  Both partners in a single rep means VMD draws the S-S bond itself,
    # so an intact disulfide and a reduced pair differ by a visible stick rather than by the
    # ~2 A gap between two spheres.  Hydrogens dropped: on a sidechain this small they are
    # clutter, and HG1 on a reduced cysteine is not what carries the distinction either.
    if {$show_ss && [ps_nsel $molid "resname CYS and sidechain and noh"] > 0} {
        ps_addrep $molid "resname CYS and sidechain and noh" \
            "Licorice 0.200000 20.000000 20.000000" "ColorID 4"
    }

    # Highlights last so they draw over the ribbon; each selection gets its own contrasting color.
    set hi 0
    foreach sel [split $highlight ";"] {
        set sel [string trim $sel]
        if {$sel eq ""} continue
        set color [lindex $PS_HIGHLIGHT_COLORS [expr {$hi % [llength $PS_HIGHLIGHT_COLORS]}]]
        set n [ps_nsel $molid $sel]
        if {$n == 0} {
            ps_note "WARNING: -highlight selection matched no atoms: $sel"
            continue
        }
        # Hydrogens off by default: a licorice highlight is read for its heavy-atom shape, and
        # the H cage around it only thickens the rep.  A selection that is ONLY hydrogens was
        # clearly meant, so honour it rather than silently drawing nothing.
        set drawsel "($sel) and noh"
        if {[ps_nsel $molid $drawsel] == 0} {
            set drawsel $sel
            ps_note "highlight [expr {$hi + 1}]: hydrogen-only selection, drawing H"
        }
        ps_addrep $molid $drawsel "Licorice 0.300000 20.000000 20.000000" "ColorID $color"
        ps_note "highlight [expr {$hi + 1}] (ColorID $color): [ps_nsel $molid $drawsel] atoms -- $drawsel"
        incr hi
    }

    # union of everything actually drawn, for the view fit
    set vis {}
    for {set i 0} {$i < [molinfo $molid get numreps]} {incr i} {
        lappend vis "([lindex [molinfo $molid get "{selection $i}"] 0])"
    }
    if {[llength $vis] == 0} { return "all" }
    set vis [join $vis " or "]
    # Frame on the solute even when solvent is drawn: fitting to the water box shrinks the thing
    # the figure is about, whereas letting the solvent run past the edge keeps it as context.
    if {$show_solvent ne "off"} {
        set solvent_sel "water or ($ion_sel)"
        if {[llength $bulk] > 0} { set solvent_sel "$solvent_sel or ($bulk_sel)" }
        set vis "($vis) and not ($solvent_sel)"
    }
    return $vis
}

# Centre the view on the *visible* atoms.  `display resetview` frames every atom in the
# molecule, so a solvated system renders its solute as a speck adrift in a mostly empty image.
proc ps_center {molid seltext} {
    if {[catch {atomselect $molid $seltext} svis]} {
        ps_note "WARNING: could not select the visible atoms; using the default view"
        return
    }
    if {[$svis num] == 0} { $svis delete ; return }
    set cen [measure center $svis]
    $svis delete
    molinfo $molid set center_matrix [list [transoffset [vecscale -1.0 $cen]]]
}

# Scale so the drawn geometry fills the frame, measured from a cheap preview render.
#
# Working out the scale analytically means knowing how VMD maps view units to pixels *and* how
# much a representation swells past the atom centres it is built from -- a cartoon ribbon, a VDW
# radius, a licorice bond all differ.  Measuring the preview's ink sidesteps both: the ink is
# exactly what will be in the final image, and the correction is a pure pixel ratio, so no view
# geometry enters.  One extra small render is far cheaper than the full-size one.
proc ps_autoscale {molid renderer fill w h} {
    set magick [ps_magick]
    if {$magick eq ""} {
        ps_note "WARNING: no ImageMagick on PATH; skipping the fill-the-frame scaling"
        return
    }
    # preview at reduced size, but keep the final aspect ratio so the fit transfers
    set pw 500
    set ph [expr {int(500.0 * $h / $w)}]
    if {$ph < 1} { set ph 1 }
    set tmp "_pssnap_preview.tga"
    # Measure the solute alone: the scale should frame what the figure is about, and including
    # the solvent's ink would shrink it back to a speck in the middle of the box.
    set hidden {}
    foreach r [ps_solvent_reps] {
        if {$r < [molinfo $molid get numreps]} { mol showrep $molid $r 0 ; lappend hidden $r }
    }
    display resize $pw $ph
    display update
    if {[catch {render $renderer $tmp} e]} {
        ps_note "WARNING: preview render failed ($e); skipping the fill-the-frame scaling"
        foreach r $hidden { mol showrep $molid $r 1 }
        display resize $w $h
        return
    }
    if {[catch {exec $magick $tmp -fuzz 2% -trim -format "%w %h %X %Y" info:} bbox]} {
        ps_note "WARNING: could not measure the preview ($bbox); skipping the scaling"
        file delete $tmp
        foreach r $hidden { mol showrep $molid $r 1 }
        display resize $w $h
        return
    }
    file delete $tmp
    foreach r $hidden { mol showrep $molid $r 1 }
    lassign $bbox iw ih ix iy
    display resize $w $h
    if {![string is integer -strict $iw] || ![string is integer -strict $ih] || $iw < 1 || $ih < 1} {
        ps_note "WARNING: preview contained no geometry; skipping the fill-the-frame scaling"
        return
    }
    # Scaling happens about the view centre, not the centre of the ink, so an off-centre object
    # (a two-domain fusion whose centroid sits inside the larger domain) grows *past* the frame
    # edge on its long side and gets clipped -- and clipped pixels are gone, the trim-and-recentre
    # afterwards cannot recover them.  Fit the ink's greatest distance from the centre instead of
    # its width, which bounds the far edge whatever the asymmetry; the recentring then removes the
    # slack that leaves on the near side.
    set ix [expr {[string is integer -strict $ix] ? $ix : 0}]
    set iy [expr {[string is integer -strict $iy] ? $iy : 0}]
    set cx [expr {$pw / 2.0}] ; set cy [expr {$ph / 2.0}]
    set dx [expr {max(abs($ix - $cx), abs($ix + $iw - $cx))}]
    set dy [expr {max(abs($iy - $cy), abs($iy + $ih - $cy))}]
    if {$dx <= 0 || $dy <= 0} { return }
    set fx [expr {$fill * $cx / $dx}]
    set fy [expr {$fill * $cy / $dy}]
    set f  [expr {$fx < $fy ? $fx : $fy}]
    scale by $f
    ps_note [format "fit: ink %dx%d at +%d+%d of %dx%d px (max offset %.0f,%.0f) -> scale by %.3f" \
                 $iw $ih $ix $iy $pw $ph $dx $dy $f]
}


# Any unit vector perpendicular to v, for the degenerate 180-degree case.
proc ps_perp {v} {
    lassign $v x y z
    if {abs($x) < 0.9} { return [vecnorm [veccross $v {1 0 0}]] }
    return [vecnorm [veccross $v {0 1 0}]]
}

# Rotate the view so *seltext* lies along *target* in the view frame ({0 0 1} = toward the
# camera, {1 0 0} = off to the side).
#
# The camera looks down -z, so +z is toward the viewer: we want the vector from the whole
# structure's centre to the feature's centre pointing along +z *after* the current view
# rotation.  Composing the correction onto the existing rotate_matrix (rather than replacing
# it) keeps whatever -view/-rot already established, so this refines an orientation instead
# of discarding it.
proc ps_orient {molid seltext fitsel target what} {
    if {$seltext eq ""} return
    if {[catch {atomselect $molid $seltext} f]} {
        ps_note "WARNING: -$what selection is invalid: $seltext" ; return
    }
    if {[$f num] == 0} {
        ps_note "WARNING: -$what selection matched no atoms: $seltext" ; $f delete ; return
    }
    if {[catch {atomselect $molid $fitsel} a]} { set a [atomselect $molid all] }
    set cf [measure center $f]
    set ca [measure center $a]
    $f delete ; $a delete
    set v [vecsub $cf $ca]
    if {[veclength $v] < 1e-6} {
        ps_note "note: -$what selection is concentric with the structure; leaving the view as is"
        return
    }
    set R0 [lindex [molinfo $molid get rotate_matrix] 0]
    # direction in the *current* view frame
    set vv [vecnorm [coordtrans $R0 $v]]
    # subtract the origin's image: coordtrans applies translation too, and a pure rotation
    # matrix has none, but be explicit so a stray translation cannot skew the direction
    set o [coordtrans $R0 {0 0 0}]
    set vv [vecnorm [vecsub [coordtrans $R0 $v] $o]]
    set z $target
    set dot [vecdot $vv $z]
    if {$dot > 0.9999} { return }
    if {$dot < -0.9999} {
        set axis [ps_perp $z] ; set ang 180.0
    } else {
        set axis [vecnorm [veccross $vv $z]]
        set ang [expr {57.2957795130823 * acos($dot)}]
    }
    molinfo $molid set rotate_matrix [list [transmult [transabout $axis $ang] $R0]]
    ps_note [format "$what: rotated %.1f deg" $ang]
}

proc ps_view {molid visible view rot zoom bg w h renderer fill face side} {
    color Display Background $bg
    axes location Off
    display projection Orthographic
    display depthcue off
    display shadows off
    display ambientocclusion off
    display resize $w $h
    display resetview
    ps_center $molid $visible

    # resetview looks down -z with x right and y up; rotating -90 about x puts z up on the
    # screen, the standard membrane side view (bilayer normal vertical)
    switch -- $view {
        xy {}
        xz { rotate x by -90 }
        yz { rotate x by -90 ; rotate y by 90 }
        default { ps_die "-view must be xy, xz, yz or auto (got '$view')" }
    }
    if {$rot ne ""} {
        foreach {rx ry rz} $rot break
        foreach {axis deg} [list x $rx y $ry z $rz] {
            if {$deg ne "" && $deg != 0} { rotate $axis by $deg }
        }
    }
    ps_orient $molid $face $visible {0 0 1} face
    ps_orient $molid $side $visible {1 0 0} side
    # autoscale after the rotations, so the fit is measured on the orientation being rendered
    display update
    ps_autoscale $molid $renderer $fill $w $h
    if {$zoom != 1.0} { scale by $zoom }
    # REQUIRED in -dispdev text: without an explicit update the scene is never handed to the
    # renderer and every Tachyon render silently writes a blank, background-colored image.
    display update
}

proc ps_render {renderer out w h bg} {
    set outext [string tolower [file extension $out]]
    # the internal Tachyon renderers write TGA; convert afterwards for any other format
    set target $out
    set convert 0
    if {$renderer ne "snapshot" && $outext ne ".tga"} {
        set convert 1
        set target "[file rootname $out].tga"
    }
    ps_note "rendering ${w}x${h} with $renderer -> $target"
    if {[catch {render $renderer $target} e]} { ps_die "render failed: $e" }
    if {![file exists $target]} { ps_die "renderer produced no file: $target" }

    set magick [ps_magick]
    if {$magick eq ""} {
        if {$convert} { ps_note "WARNING: no ImageMagick on PATH; leaving $target as-is" }
        return $target
    }
    # Centre the drawn geometry on the canvas.  The view is centred on the visible atoms'
    # centroid, which is not the centre of the *ink* (a lopsided system, or a cartoon that
    # swells asymmetrically, sits off-centre).  Trimming and re-extending to the same canvas
    # size only crops and pads -- no resampling -- so nothing is softened.
    if {[catch {exec $magick $target -fuzz 2% -trim +repage \
                    -background $bg -gravity center -extent ${w}x${h} $out} e]} {
        ps_note "WARNING: could not centre the image ($e); writing it uncentred"
        if {$convert} {
            if {[catch {exec $magick $target $out} e2]} { ps_die "image conversion failed: $e2" }
        } else {
            return $target
        }
    }
    if {$target ne $out} { file delete $target }
    return $out
}

proc ps_main {argv} {
    global PS_GLYCAN_RESNAMES PS_LIPID_EXTRA
    array set opt [ps_parse_args $argv]

    if {$opt(psf) eq ""}  { ps_die "no -psf given" }
    if {$opt(coor) eq ""} { ps_die "no -coor given" }
    foreach f [list $opt(psf) $opt(coor)] {
        if {![file exists $f]} { ps_die "file not found: $f" }
    }
    if {![regexp {^([0-9]+)x([0-9]+)$} $opt(size) -> w h]} {
        ps_die "-size must look like 1600x1200 (got '$opt(size)')"
    }

    ps_load molid $opt(psf) $opt(coor) $opt(frame)

    set n_lipid  [ps_nsel $molid "lipid or resname $PS_LIPID_EXTRA"]
    set n_glycan [ps_nsel $molid "resname $PS_GLYCAN_RESNAMES"]
    ps_note "content: protein [ps_nsel $molid protein], nucleic [ps_nsel $molid nucleic],\
             lipid $n_lipid, glycan $n_glycan, water [ps_nsel $molid water]"

    if {$opt(style) eq "auto"} {
        if {$n_lipid > 0}        { set opt(style) membrane
        } elseif {$n_glycan > 0} { set opt(style) glycan
        } else                   { set opt(style) protein }
        ps_note "style: auto -> $opt(style)"
    }
    if {$opt(view) eq "auto"} {
        set opt(view) [expr {$n_lipid > 0 ? "xz" : "xy"}]
        ps_note "view: auto -> $opt(view)"
    }

    set opt(solvent) [ps_solvent_mode $opt(solvent)]
    set visible [ps_represent $molid $opt(style) $opt(solvent) $opt(highlight) $opt(ss) $opt(ghost) $opt(domains)]
    ps_view $molid $visible $opt(view) $opt(rot) $opt(zoom) $opt(bg) $w $h $opt(renderer) $opt(fill) $opt(face) $opt(side)
    set written [ps_render $opt(renderer) $opt(o) $w $h $opt(bg)]
    ps_note "wrote $written"
}

ps_main $argv
exit 0
