# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# General VMD startup script for pestifer-managed VMD sessions or
# VMD sessions in initiated with the statement
#
# source [pestifer_init]
#
# and run in a conda environment in which pestifer is installed
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
  set PESTIFER_TCLROOT [exec pestifer --no-banner wheretcl --root]
  if {$quiet == 0} {
    vmdcon -info "No PESTIFER_TCLROOT passed in; detected $PESTIFER_TCLROOT"
  }
}

# set up package path
set PESTIFER_TCLPKG ${PESTIFER_TCLROOT}/pkg
set packages [glob -type d $PESTIFER_TCLPKG]
foreach pkg $packages {
  lappend auto_path $pkg
}

# import the PestiferUtil general purpose package
package require PestiferUtil
namespace import PestiferUtil::*

# source the current atomselect macro definitions
if {[file exists ${PESTIFER_TCLROOT}/macros.tcl]} {
  source ${PESTIFER_TCLROOT}/macros.tcl
}
