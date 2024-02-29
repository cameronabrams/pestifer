# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# General VMD startup script for pestifer
#
# This name of this file is set to the vmd_startup_script attribute
# of the Config instance of a pestifer run, by default (config.py)
#
set RESPATH "."
set quiet 1
for { set a 0 } { $a < [llength $argv] } { incr a } {
  set arg [lindex $argv $a]
  if { $arg == "-respath" } {
    incr a
    set RESPATH [lindex $argv $a]
  }
  if { $arg == "-verbose" } {
    incr a
    set quiet 0
  }
}


if {$RESPATH=="."} {
  set RESPATH [exec pestifer wheretcl --proc-dir]
  if {$quiet == 0} {
    vmdcon -info "No RESPATH passed in; detected $RESPATH"
  }
}

set sources [glob "${RESPATH}/*.tcl"]

if {$quiet == 0} {
  vmdcon -info "Sourcing ${RESPATH}/util.tcl"
}
# must source this first
source ${RESPATH}/util.tcl
# now remove it from the globbed list; vmd does not have "lremove"
set lidx [lsearch $sources "${RESPATH}/util.tcl"]
set tmp_sources [list]
for {set i 0} {$i < [llength $sources]} {incr i} {
  if {$i != $lidx} {
    lappend tmp_sources [lindex $sources $i]
  }
}
set sources $tmp_sources

foreach s $sources {
  source $s
  if {$quiet == 0} {
    vmdcon -info "Sourcing $s"
  }
}
