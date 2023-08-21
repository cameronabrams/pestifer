set RESPATH "."
for { set a 0 } { $a < [llength $argv] } { incr a } {
  set arg [lindex $argv $a]
  if { $arg == "-respath" } {
    incr a
    set RESPATH [lindex $argv $a]
  }
}
source ${RESPATH}/proc/axeq.tcl
source ${RESPATH}/proc/crot.tcl
source ${RESPATH}/proc/declashloop.tcl
source ${RESPATH}/proc/checkpierce.tcl
source ${RESPATH}/proc/dcdlog.tcl