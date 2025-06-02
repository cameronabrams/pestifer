# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# General VMD startup script for pestifer-managed VMD sessions or
# VMD sessions in initiated with the statement
#
# pestifer_init
#
# and run in a conda environment in which pestifer is installed
# 

# we are calling this script at vmd startup so we can pass arguments to it
global auto_path

if {[info exists argv]} {
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
}
global quiet
if {![info exists quiet]} {
  vmdcon -info "No quiet argument passed in; defaulting to 0"
  set quiet 1
}

global PESTIFER_TCLROOT
if {![info exists PESTIFER_TCLROOT]} {
  set PESTIFER_TCLROOT [exec pestifer --no-banner wheretcl --root]
  if {$quiet == 0} {
    vmdcon -info "No PESTIFER_TCLROOT passed in; detected $PESTIFER_TCLROOT"
  }
}

# set up package path
set PESTIFER_TCLPKG ${PESTIFER_TCLROOT}/pkg
set packages [glob -type d $PESTIFER_TCLPKG]
foreach pkg $packages {
  vmdcon -info "Adding $pkg to auto_path"
  lappend auto_path $pkg
}

# import the PestiferUtil general purpose package
package require PestiferUtil
namespace import PestiferUtil::*

# source the current atomselect macro definitions
if {[file exists ${PESTIFER_TCLROOT}/macros.tcl]} {
  source ${PESTIFER_TCLROOT}/macros.tcl
}
