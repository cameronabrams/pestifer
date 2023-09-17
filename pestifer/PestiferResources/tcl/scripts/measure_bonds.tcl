# measure_bonds.tcl -- 
# reads lines from input file, each is CHAIN, RES1, RES2
# measures distance between C on RES1 and N on RES2
# appends that distance to the lines
# writes all lines back out to file

vmdcon -info "measure_bonds.tcl: args: $argv"
set rfzr 0 ;# receiver flexibility zone radius
for { set i 0 } { $i < [llength $argv] } { incr i } {
    if { [lindex $argv $i] == "-pdb"} {
       incr i
       set coor [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-psf"} {
       incr i
       set psf [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-opdb"} {
       incr i
       set fixed [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-i"} {
       incr i
       set infile [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-o"} {
       incr i
       set outfile [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-rfzr"} {
       incr i
       set rfzr [lindex $argv $i]
    }
}

mol new $psf
mol addfile $coor
set a [atomselect top "all"]
$a set occupancy 0
set ca [atomselect top "name CA or (water) or (glycan and element C O)"]
$ca set occupancy 1

set fp [open $infile "r"]
set lines [split [read $fp] \n]
close $fp
set RES0 {}
set RES1 {}
set RES2 {}
set SEG {}
foreach l $lines {
    if { [llength $l] == 4 } {
       lappend SEG [lindex $l 0]
       lappend RES0 [lindex $l 1]
       lappend RES1 [lindex $l 2]
       lappend RES2 [lindex $l 3]
    }
}
puts "SEG $SEG"
puts "RES0 $RES0"
puts "RES1 $RES1"
puts "RES2 $RES2"

set bl {}
set CC {}
set NN {}
foreach seg $SEG h $RES0 i $RES1 j $RES2 {
    set labelus [atomselect top "segname $seg and resid $h to $i"]
    $labelus set occupancy 0
    set ci [string index $i end]
    if {[string is alpha $ci]} {
        set resid [string range $i 0 end-1]
        set theC [atomselect top "segname $seg and resid $resid and insertion $ci and name C"]
    } else {
        set theC [atomselect top "segname $seg and resid $i and name C"]
    }
    set ni [string index $j end]
    if {[string is alpha $ni]} {
        set resid [string range $j 0 end-1]
        set theN [atomselect top "segname $seg and resid $resid and insertion $ni and name N"]
    } else {
        set theN [atomselect top "segname $seg and resid $j and name N"]
    }
    if {[expr $rfzr > 0]} {
        set theN_idx [$theN get index]
        set flex_zone [atomselect top "segname $seg and name CA and same residue as (exwithin $rfzr of index $theN_idx)"]
        $flex_zone set occupancy 0
    }
    puts "resid $i on segname $seg has [$theC num] Cs with serial [$theC get serial]"
    puts "resid $j on segname $seg has [$theN num] Ns with serial [$theN get serial]"
    set ii [$theC get index]
    set jj [$theN get index]
    lappend bl [measure bond [list $ii $jj]]
    lappend CC [$theC get serial]
    lappend NN [$theN get serial]
}
puts "CC $CC"
puts "NN $NN"

set fp [open $outfile "w"]
puts $fp "# segname c-serial n-serial distance(A)"
foreach s $SEG i $CC j $NN b $bl {
    puts $fp "$s $i $j [format %.4f $b]"
}
close $fp
$a writepdb $fixed
vmdcon -info "measure_bonds.tcl: $outfile created and $fixed created"
exit
