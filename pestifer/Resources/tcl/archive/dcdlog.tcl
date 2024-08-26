# logs frames to a logging molecule
proc log_addframe { molid logid } {
   if { $logid != "-1" } {
     animate dup $logid
     [atomselect $logid all] set x [[atomselect $molid all] get x]
     [atomselect $logid all] set y [[atomselect $molid all] get y]
     [atomselect $logid all] set z [[atomselect $molid all] get z]
#     [atomselect $molid all] writepdb "tmp.pdb"
#     animate read pdb tmp.pdb $logid
     vmdcon -info "Molid $molid - logging molecule $logid has [molinfo $logid get numframes] frames."
#     exec rm -f tmp.pdb
   }
}