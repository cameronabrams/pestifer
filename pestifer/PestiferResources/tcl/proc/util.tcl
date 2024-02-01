# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Utility procs for VMD/Psfgen TcL scripts
#

# Logs a message to both the VMD console and to an open
# log file
proc double_log { logf message } {
   vmdcon -info "$message"
   puts $logf "$message"
   flush $logf
}

# outputs the message if num is zero
proc checknum { num msg } {
  if { $num == 0 } { 
    vmdcon -info "$msg"
  }
}

# returns -1 if x is negative, 1 otherwise
proc sign { x } {
    if { $x < 0 } {
        return -1
    } else {
        return 1
    }
}

# returns the unit vector of an input vector
proc hatvec { a_vec } {
   set l [veclength $a_vec]
   return [vecscale $a_vec [expr 1.0/$l]]
}

# Wrap x into periodic domain [xlo,xhi]
proc wrap_domain { x xlo xhi } {
   set dsize [expr ($xhi) - ($xlo)]
   if { [expr $x < $xlo] } {
      set y [expr $x + $dsize]
   } elseif { [expr $x > $xhi] } {
      set y [expr $x - $dsize]
   } else {
      set y $x
   }
   return $y
}