if {![info exists TCL_DIR]} {
  # see if user set an environment variable
  if {[info exists env(TCL_DIR)]} {
      set TCL_DIR $env(TCL_DIR)
      puts "TCL_DIR $TCL_DIR"
  } else {
      set TCL_DIR $env(HOME)/Git/packlearn/
  }
}

source ${TCL_DIR}/list_handlers.tcl

load ${TCL_DIR}/modules/lib/bondstruct.so cfa_bondstruct

set il [list 0 1 2 3]
set bondcount 0
set ia [intListToArray $il]
set bs [new_bondstruct $ia [llength $il] $bondcount]

puts "new bondstruct at $bs"

exit
