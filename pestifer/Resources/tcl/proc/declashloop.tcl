
# A very simple Metropolis algorithm to alter conformation of raw model-built
# loops to minimize clashes
proc declash_loop { molid segname loop maxcycles } {
  set nr [llength $loop]
  set loopsel [atomselect $molid "segname $segname and resid [join $loop]"]
  set residue_numbers [[atomselect $molid "[$loopsel text] and name CA"] get residue]
  set env [atomselect $molid "same residue as exwithin 4.0 of (segname $segname and resid [join $loop])"]
  set residuenum_end [lindex $residue_numbers end]
  puts "DECLASH_LOOP) molid $molid segname $segname loop $loop maxcycles $maxcycles"
  puts "DECLASH_LOOP) loopsel [$loopsel num] atoms; residue_numbers $residue_numbers"
  puts "DECLASH_LOOP) env [$env num] atoms"
  for { set i 0 } { $i < $nr } { incr i } {
    # rotate phi angle and psi angle to minimize number of contacts between residue and 
    # its environment
    set rsel [atomselect $molid "segname $segname and resid [lindex $loop $i] to [lindex $loop end]"]
    set residuenum1 [lindex $residue_numbers $i]
    set CON_STRUCT [measure contacts 2.0 $rsel $env]
    set CON [llength [lindex $CON_STRUCT 0]]
    puts "DECLASH_LOOP) [lindex $loop $i] INIT $CON"
    for { set t 0 } { $t < $maxcycles && $CON > 0 } { incr t } {
      set SAVEPOS [$loopsel get {x y z}]
      set rphi [expr (1-2*rand())*120.0]
      #set rpsi [expr (1-2*rand())*60.0]
      Crot_phi_toCterm $residuenum1 $residuenum_end $segname $molid $rphi
      if { $i > [expr $nr - 1] } {
        Crot_psi_toCterm $residuenum1 $residuenum_end $segname $molid $rphi
      }
      $env update
      set TRICON_STRUCT [measure contacts 2.0 $rsel $env]
      set TRICON  [llength [lindex $TRICON_STRUCT 0]]
      if { [expr $TRICON < $CON] } {
        # accept this move
        set CON $TRICON
        puts "DECLASH_LOOP) ${segname}:[lindex $loop 0]-[lindex $loop $i] $t $CON"
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
