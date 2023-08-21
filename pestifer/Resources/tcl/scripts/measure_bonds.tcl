# measure_bonds.tcl -- 
# reads lines from input file, each is CHAIN, RES1, RES2
# measures distance between C on RES1 and N on RES2
# appends that distance to the lines
# writes all lines back out to file

vmdcon -info "measure_bonds.tcl: args: $argv"
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
       set fixed [lindex $argv $i]
    }
    if { [lindex $argv $i] == "-i"} {
       incr i
       set infile [lindex $argv $i]
    }
}

mol new $psf
mol addfile $coor
set a [atomselect top "all"]
$a set occupancy 0
set ca [atomselect top "name CA"]
$ca set occupancy 1

set fp [open $infile "r"]
set lines [split [read $fp] \n]
close $fp
set RES0 {}
set RES1 {}
set RES2 {}
set CH {}
foreach l $lines {
    if { [llength $l] == 4 } {
       lappend CH [lindex $l 0]
       lappend RES0 [lindex $l 1]
       lappend RES1 [lindex $l 2]
       lappend RES2 [lindex $l 3]
    }
}
puts "CH $CH"
puts "RES0 $RES0"
puts "RES1 $RES1"
puts "RES2 $RES2"

set bl {}
set CC {}
set NN {}
foreach c $CH h $RES0 i $RES1 j $RES2 {
    set labelus [atomselect top "chain $c and resid $h to $i"]
    $labelus set occupancy 0
    set ci [string index $i end]
    if {[string is alpha $ci]} {
        set resid [string range $i 0 end-1]
        set theC [atomselect top "chain $c and resid $resid and insertion $ci and name C"]
    } else {
        set theC [atomselect top "chain $c and resid $i and name C"]
    }
    set ni [string index $j end]
    if {[string is alpha $ni]} {
        set resid [string range $j 0 end-1]
        set theN [atomselect top "chain $c and resid $resid and insertion $ni and name N"]
    } else {
        set theN [atomselect top "chain $c and resid $j and name N"]
    }
    puts "resid $i on chain $c has [$theC num] Cs with serial [$theC get serial]"
    puts "resid $j on chain $c has [$theN num] Ns with serial [$theN get serial]"
    set ii [$theC get index]
    set jj [$theN get index]
    lappend bl [measure bond [list $ii $jj]]
    lappend CC [$theC get serial]
    lappend NN [$theN get serial]
}
puts "CC $CC"
puts "NN $NN"

set fp [open $infile "w"]
foreach c $CH i $CC j $NN b $bl {
    puts $fp "$c $i $j [format %.4f $b]"
}
close $fp
$a writepdb $fixed
vmdcon -info "measure_bonds.tcl: $infile modified and $fixed created"
exit
