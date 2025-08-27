
pestifer_init

proc abc { arg } {
    puts "arg: $arg"
}
catch {set a "hello world"}
for { set i 0 } { $i < [llength $argv] } { incr i } {
   if { [lindex $argv $i] == "-a"} {
      incr i
      set a [lindex $argv $i]
   }
}
set a_arg [deprotect_str_arg $a]
abc $a_arg
exit
