# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# General VMD startup script for pestifer
#
# This name of this file is set to the vmd_startup_script attribute
# of the Config instance of a pestifer run, by default (config.py)
#
set PESTIFER_TCLROOT "."
set quiet 1
for { set a 0 } { $a < [llength $argv] } { incr a } {
  set arg [lindex $argv $a]
  if { $arg == "--tcl-root" } {
    incr a
    set PESTIFER_TCLROOT [lindex $argv $a]
  }
  if { $arg == "--verbose" } {
    incr a
    set quiet 0
  }
}

if {$PESTIFER_TCLROOT=="."} {
  set PESTIFER_TCLROOT [exec pestifer wheretcl --root]
  if {$quiet == 0} {
    vmdcon -info "No PESTIFER_TCLROOT passed in; detected $PESTIFER_TCLROOT"
  }
}

set PESTIFER_TCLPKG ${PESTIFER_TCLROOT}/pkg
set packages [glob -type d $PESTIFER_TCLPKG]
foreach pkg $packages {
  lappend auto_path $pkg
}

package require PestiferUtil
namespace import PestiferUtil::*

# these should all be put into packages

set PESTIFER_TCLPROC ${PESTIFER_TCLROOT}/proc

set sources [glob "${PESTIFER_TCLPROC}/*.tcl"]

foreach s $sources {
  source $s
  if {$quiet == 0} {
    vmdcon -info "Sourcing $s"
  }
}
