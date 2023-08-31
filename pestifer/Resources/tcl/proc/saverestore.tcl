proc backup { atomsel attributes } {
    set data [list]
    foreach attr $attributes {
        lappend data [$atomsel get $attr]
    }
    return $data
}

proc restore { atomsel attributes data } {
    foreach attr $attributes values $data {
        $atomsel set $attr $values
    }
}