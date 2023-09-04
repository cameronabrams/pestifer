set RESPATH "."
for { set a 0 } { $a < [llength $argv] } { incr a } {
  set arg [lindex $argv $a]
  if { $arg == "-respath" } {
    incr a
    set RESPATH [lindex $argv $a]
  }
}
source ${RESPATH}/axeq.tcl
source ${RESPATH}/caco.tcl
source ${RESPATH}/crot.tcl
source ${RESPATH}/declashloop.tcl
source ${RESPATH}/checkpierce.tcl
source ${RESPATH}/dcdlog.tcl
source ${RESPATH}/saverestore.tcl
