# pestifer: testtcl.tcl
####################### Created Tue May 27 17:49:17 2025 #######################
mol new 6pti.pdb
set a [atomselect top all]
set data [ backup $a [ list chain x y z resid resname name ] ]
[atomselect top "name CA"] set x 0
[atomselect top "name CA"] set y -1
[atomselect top "name CA"] set z 1
[atomselect top "name CA"] writepdb bad.pdb
restore $a [ list chain x y z resid resname name ]  $data
[atomselect top "name CA"] writepdb good.pdb
exit
########################### END PESTIFER VMD SCRIPT ############################
######################## Thank you for using pestifer! #########################
