wheretcl
--------

Pestifer has a pretty handy library of TcL packages.  If you want to peruse the source, pestifer will tell you where to find them:

.. code-block:: console

   $ pestifer wheretcl --pkg-dir

If you want to use any of the procs defined in those packages in your own VMD script, the easiest thing to do is to put this ``proc`` definition in your own VMD startup file:

.. code-block:: tcl

      proc pestifer_init { } {
         set status 0
         if {[catch {exec which pestifer} results options]} {
            set details [dict get $options -errorcode]
            if {[lindex $details 0] eq "CHILDSTATUS"} {
               set status [lindex $details 2]
            } else {
               return -options $options -level 0 $results
            }
         }
         if { $status == 0 } {
            set pestifer_tcl_root [exec pestifer --no-banner wheretcl --root]
            vmdcon -info "Source ${pestifer_tcl_root}/vmdrc.tcl"
            return ${pestifer_tcl_root}/vmdrc.tcl
         } else {
            vmdcon -info "Pestifer is not available in your current environment."
         }
      }

Then, you can use it in a source command in any VMD script or TcL session you like:

.. code-block:: tcl

   source [pestifer_init]

This of course requires that your VMD session was launched from a shell running a python virtual environment in which ``pestifer`` is installed.