# matrix inversion routine from http://wiki.tcl.tk/14921 (Keith Vetter)
proc Inverse3 {matrix} {
  if {[llength $matrix] != 3 ||
    [llength [lindex $matrix 0]] != 3 || 
    [llength [lindex $matrix 1]] != 3 || 
    [llength [lindex $matrix 2]] != 3} {
     error "wrong sized matrix"
  }
  set inv {{? ? ?} {? ? ?} {? ? ?}}
 
  # Get adjoint matrix : transpose of cofactor matrix
  for {set i 0} {$i < 3} {incr i} {
    for {set j 0} {$j < 3} {incr j} {
      lset inv $i $j [_Cofactor3 $matrix $i $j]
    }
  }
  # Now divide by the determinant
  set det [expr {double([lindex $matrix 0 0]   * [lindex $inv 0 0]
                 + [lindex $matrix 0 1] * [lindex $inv 1 0]
                 + [lindex $matrix 0 2] * [lindex $inv 2 0])}]
  if {$det == 0} {
    error "non-invertable matrix"
  }
    
  for {set i 0} {$i < 3} {incr i} {
    for {set j 0} {$j < 3} {incr j} {
      lset inv $i $j [expr {[lindex $inv $i $j] / $det}]
    }
  }
  return $inv
}

proc _Cofactor3 {matrix i j} {
  array set COLS {0 {1 2} 1 {0 2} 2 {0 1}}
  foreach {row1 row2} $COLS($j) break
  foreach {col1 col2} $COLS($i) break
    
  set a [lindex $matrix $row1 $col1]
  set b [lindex $matrix $row1 $col2]
  set c [lindex $matrix $row2 $col1]
  set d [lindex $matrix $row2 $col2]
 
  set det [expr {$a*$d - $b*$c}]
  if {($i+$j) & 1} { set det [expr {-$det}]}
  return $det
}

# returns the position of an amide N given the upstream position
# of a CA, C, and O.  This is necessary because psfgen does not know
# how to use IC's to generate de novo the position of the N of
# a modeled in residue off of a residue with an arbitrary position.
proc cacoIn_nOut { resid segname molid } {
  set rn {? ? ?}

  set ca [atomselect $molid "segname $segname and resid $resid and name CA"]
  set c  [atomselect $molid "segname $segname and resid $resid and name C"]
  set o  [atomselect $molid "segname $segname and resid $resid and name O"]

  set r1 [lindex [$ca get {x y z}] 0]
  set r2 [lindex [$c  get {x y z}] 0]
  set r3 [lindex [$o  get {x y z}] 0]

  set R21 [vecsub $r2 $r1]
  set r21 [vecscale $R21 [expr 1.0/[veclength $R21]]]
  set R32 [vecsub $r3 $r2]
  set r32 [vecscale $R32 [expr 1.0/[veclength $R32]]]
  set c [veccross $r21 $r32]

  set MAT { {? ? ?} {? ? ?} {? ? ?} }
  lset MAT 0 0 [lindex $c 0]
  lset MAT 0 1 [lindex $c 1]
  lset MAT 0 2 [lindex $c 2]
  lset MAT 1 0 [lindex $r21 0]
  lset MAT 1 1 [lindex $r21 1]
  lset MAT 1 2 [lindex $r21 2]
  lset MAT 2 0 [lindex $r32 0]
  lset MAT 2 1 [lindex $r32 1]
  lset MAT 2 2 [lindex $r32 2]

#  puts "MAT $MAT"
#  exit
  set b { ? ? ? }
  lset b 0 0
  lset b 1 [expr -cos(3.14159*114.44/180.0)] 
  lset b 2 [expr cos(3.14159*123.04/180.0)]
  
  set AMAT [Inverse3 $MAT]

  set rnd { ? ? ? }
  lset rnd 0 [expr [lindex $AMAT 0 0] * [lindex $b 0] + [lindex $AMAT 0 1] * [lindex $b 1] + [lindex $AMAT 0 2] * [lindex $b 2]]
  lset rnd 1 [expr [lindex $AMAT 1 0] * [lindex $b 0] + [lindex $AMAT 1 1] * [lindex $b 1] + [lindex $AMAT 1 2] * [lindex $b 2]]
  lset rnd 2 [expr [lindex $AMAT 2 0] * [lindex $b 0] + [lindex $AMAT 2 1] * [lindex $b 1] + [lindex $AMAT 2 2] * [lindex $b 2]]

  set rn [vecadd $r2 [vecscale $rnd 1.355]]
  vmdcon -info "caco $rn"
  return $rn
}
