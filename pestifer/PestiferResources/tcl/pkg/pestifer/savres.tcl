# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#

package provide SavRes 1.0

namespace eval ::SavRes:: {
    namespace export *
}

proc SavRes::backup { atomsel attributes } {
    set data [list]
    foreach attr $attributes {
        lappend data [$atomsel get $attr]
    }
    return $data
}

proc SavRes::restore { atomsel attributes data } {
    foreach attr $attributes values $data {
        $atomsel set $attr $values
    }
}