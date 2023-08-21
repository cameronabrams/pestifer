proc checknum { num msg } {
  if { $num == 0 } { 
    puts "$msg"
    #exit
  }
}
# A very simple Metropolis algorithm to alter conformation of raw model-built
# loops to minimize clashes
proc declash_loop { molid c loop maxcycles } {
  set nr [llength $loop]
  set loopsel [atomselect $molid "chain $c and resid [join $loop]"]
  set residue_numbers [[atomselect $molid "[$loopsel text] and name CA"] get residue]
  set env [atomselect $molid "same residue as exwithin 4.0 of (chain $c and resid [join $loop])"]
  set residuenum_end [lindex $residue_numbers end]
  for { set i 0 } { $i < $nr } { incr i } {
    # rotate phi angle and psi angle to minimize number of contacts between residue and 
    # its environment
    set rsel [atomselect $molid "chain $c and resid [lindex $loop $i] to [lindex $loop end]"]
    set residuenum1 [lindex $residue_numbers $i]
    set CON_STRUCT [measure contacts 2.0 $rsel $env]
    set CON [llength [lindex $CON_STRUCT 0]]
#    puts "LAYLOOP) ${c}[lindex $loop $i] INIT $CON"
    for { set t 0 } { $t < $maxcycles && $CON > 0 } { incr t } {
      set SAVEPOS [$loopsel get {x y z}]
      set rphi [expr (1-2*rand())*120.0]
      #set rpsi [expr (1-2*rand())*60.0]
      Crot_phi_toCterm $residuenum1 $residuenum_end $c $molid $rphi
      if { $i > [expr $nr - 1] } {
        Crot_psi_toCterm $residuenum1 $residuenum_end $c $molid $rphi
      }
      $env update
      set TRICON_STRUCT [measure contacts 2.0 $rsel $env]
      set TRICON  [llength [lindex $TRICON_STRUCT 0]]
      if { [expr $TRICON < $CON] } {
        # accept this move
        set CON $TRICON
        puts "LAYLOOP) ${c}:[lindex $loop 0]-[lindex $loop $i] $t $CON"
      } else {
        # reject this move
        $loopsel set {x y z} $SAVEPOS
      }
    }
    $rsel delete
  }
  $env delete
  $loopsel delete
}
