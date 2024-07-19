These pdb files are lipids for which packmol-memgen does not have CHARMM-compatible pdb files
and which are not in packmol-memgen's charmmlipid2amber subpackage database.

> Cameron F. Abrams, <cfa22@drexel.edu>


Below is an example psfgen script that verifies these pdb files
are compatible with the CHARMM force-field and parameter files
```
package require psfgen
topology ../toppar/top_all36_lipid.rtf
topology ../toppar/stream/lipid/toppar_all36_lipid_miscellaneous.str

segment A {
    pdb SOPE.pdb
}
coordpdb SOPE.pdb A 

writepsf test.psf 
writepdb test.pdb
```

and

```
package require psfgen
topology ../toppar/top_all36_lipid.rtf
topology ../toppar/stream/lipid/toppar_all36_lipid_miscellaneous.str

segment A {
    pdb SOPS.pdb
}
coordpdb SOPS.pdb A 

writepsf test.psf 
writepdb test.pdb
```

Using either of these systems in a NAMD script results in a stable simulation.

```
structure test.psf
coordinates test.pdb
parameters ../toppar/par_all36_lipid.prm
parameters ../toppar/stream/lipid/toppar_all36_lipid_miscellaneous.str

set temperature 300
paraTypeCharmm True
exclude scaled1-4
1-4scaling 1.0
cutoff 10.0
switching True
switchdist 9.0
pairlistdist 11.5
outputenergies 100
temperature $temperature
nonbondedFreq 1
fullElectFrequency 2
stepspercycle 4
dielectric 80
timestep 1.0
rigidbonds none
outputName test-out
dcdfreq 100
firsttimestep 0

minimize 100
run 1000
```