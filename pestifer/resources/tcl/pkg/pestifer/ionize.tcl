## PestiferIonize -- add ions to a non-water solvated system by replacing solvent
## molecules with ions, modeled on VMD's ``autoionize`` (which only works for water).
##
## VMD ``autoionize`` places ions by carving space out of ``water and noh``; a non-aqueous
## solvent box therefore cannot be ionized with it.  ``ionize_by_replacement`` does the same
## job for an arbitrary solvent: it selects solvent molecules (by their key atom) that sit far
## enough from the solute and from one another, deletes them, and puts an ion at each freed
## position via psfgen.
##
## Author: Cameron F. Abrams, <cfa22@drexel.edu>
package provide PestiferIonize 1.0
package require psfgen

namespace eval ::PestiferIonize:: {
    namespace export *
}

# ionize_by_replacement psf pdb -solvent RESN -keyatom NAME -ions {RESN N ...} -topology FILE
#                        [-o PREFIX] [-from D] [-between D] [-seg NAME] [-maxtries N]
#
# Writes <prefix>.psf / <prefix>.pdb.  Errors (does not exit) if it cannot place the ions.
proc ::PestiferIonize::ionize_by_replacement {psf pdb args} {
    array set opt {-from 5.0 -between 5.0 -seg ION -maxtries 10 \
                   -solvent {} -keyatom {} -o ionized -topology {} -ions {}}
    array set opt $args

    set solvent $opt(-solvent)
    set keyatom $opt(-keyatom)
    set prefix  $opt(-o)
    set from    $opt(-from)
    set between $opt(-between)
    set segname $opt(-seg)
    set maxTries $opt(-maxtries)
    set topfile $opt(-topology)
    set ionspec $opt(-ions)

    set nIons 0
    foreach {rn cnt} $ionspec { incr nIons $cnt }
    if {$nIons == 0} {
        puts "PestiferIonize) No ions requested; nothing to do."
        return
    }
    puts "PestiferIonize) Placing $nIons ion(s) by replacing $solvent molecules (key atom $keyatom)."

    # load a molecule for the geometric selection
    set molid [mol new $psf type psf waitfor all]
    mol addfile $pdb type pdb waitfor all molid $molid

    # pick nIons solvent key atoms that are >from from the solute and >between from each other
    set nTries 0
    set ionList {}
    while {1} {
        set ionList {}
        set sel [atomselect $molid "name $keyatom and resname $solvent and not (within $from of (not resname $solvent))"]
        set candIndex [$sel get index]
        $sel delete
        set candSize [llength $candIndex]
        while {[llength $ionList] < $nIons} {
            if {$candSize == 0} { break }
            set thisIon [lindex $candIndex [expr {int($candSize * rand())}]]
            set clash 0
            if {[llength $ionList]} {
                set tempsel [atomselect $molid "index [concat $ionList] and within $between of (index $thisIon)"]
                set clash [$tempsel num]
                $tempsel delete
            }
            if {$clash == 0} { lappend ionList $thisIon }
        }
        if {[llength $ionList] == $nIons} { break }
        incr nTries
        if {$nTries >= $maxTries} {
            mol delete $molid
            error "PestiferIonize) could not place $nIons ions in $solvent after $maxTries tries; lower -from/-between, reduce the ion count, or use a larger box"
        }
    }

    # record ion positions and which solvent residues to delete
    set sel [atomselect $molid "index $ionList"]
    set ionPos   [$sel get {x y z}]
    set delSegids [$sel get segid]
    set delResids [$sel get resid]
    $sel delete
    mol delete $molid

    # rebuild the system with psfgen: delete the chosen solvent residues, add an ion segment
    resetpsf
    readpsf $psf
    coordpdb $pdb
    foreach segid $delSegids resid $delResids {
        delatom $segid $resid
    }
    if {$topfile ne {}} { topology $topfile }
    set resid 1
    segment $segname {
        first NONE
        last NONE
        foreach {rn cnt} $ionspec {
            for {set i 0} {$i < $cnt} {incr i} {
                residue $resid $rn
                incr resid
            }
        }
    }
    # place one ion at each freed position (monatomic ion: atom name == resname)
    set resid 1
    set p 0
    foreach {rn cnt} $ionspec {
        for {set i 0} {$i < $cnt} {incr i} {
            coord $segname $resid $rn [lindex $ionPos $p]
            incr resid
            incr p
        }
    }
    writepsf $prefix.psf
    writepdb $prefix.pdb
    resetpsf
    puts "PestiferIonize) Added $nIons ion(s) to segment $segname; wrote $prefix.psf / $prefix.pdb."
}
