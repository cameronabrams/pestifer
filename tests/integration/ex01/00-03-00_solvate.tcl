# pestifer.scripters: 00-03-00_solvate.tcl
####################### Created Mon Aug 25 22:57:05 2025 #######################
package require solvate
package require autoionize
psfcontext mixedcase
mol new 00-01-00_psfgen-build.psf
mol addfile 00-02-00_md-minimize.pdb waitfor all
solvate 00-01-00_psfgen-build.psf 00-02-00_md-minimize.pdb -minmax { { -19.66 -20.808 -22.180999999999997 } { 38.327 30.366 29.744 } } -o 00-03-00_solvate_solv
autoionize -psf 00-03-00_solvate_solv.psf -pdb 00-03-00_solvate_solv.pdb -neutralize -o 00-03-00_solvate
exit
###################### END PESTIFER.SCRIPTERS VMD SCRIPT #######################
################### Thank you for using pestifer.scripters! ####################
