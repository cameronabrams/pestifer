###################### PESTIFER SCRIPT step1-layloops.tcl ######################
####################### Created Fri Aug 18 17:32:43 2023 #######################
mol new step1-psfgen.psf
mol addfile step1-psfgen.pdb waitfor all
set mLL [molinfo top get id]
lay_loop $mLL G [list 185A 185B 185C 185D 185E 185F 185G 185H 185I] 200
lay_loop $mLL J [list 185A 185B 185C 185D 185E 185F 185G 185H 185I] 200
lay_loop $mLL O [list 185A 185B 185C 185D 185E 185F 185G 185H 185I] 200
lay_loop $mLL G [list 400 401 402 403 404 405 406 407 408 409 410] 200
lay_loop $mLL J [list 400 401 402 403 404 405 406 407 408 409 410] 200
lay_loop $mLL O [list 400 401 402 403 404 405 406 407 408 409 410] 200
lay_loop $mLL B [list 548 549 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568] 200
lay_loop $mLL F [list 548 549 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568] 200
lay_loop $mLL L [list 548 549 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568] 200
set X [atomselect $mLL all]
$X writepdb step1-layloops.pdb
exit
########################### END PESTIFER VMD SCRIPT ############################
######################## Thank you for using pestifer! #########################
