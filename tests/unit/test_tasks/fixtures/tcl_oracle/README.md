# Tcl orientation oracle (test-only)

These files are **not** part of the pestifer distribution.  `tests/**` is excluded from both
the wheel and the sdist, so nothing here ships to users.

- `bilayer_orient.tcl` -- the retired VMD script that oriented a transmembrane protein so that
  z is the membrane normal.
- `pkg/orient/` -- the third-party `Orient` package it requires (from the VMD folks).
- `pkg/la1.0/` -- the third-party `La` linear-algebra package that `Orient` requires
  (Hume Integration Software).

Orientation now happens in Python: `pestifer.objs.rottrans.RotTrans` with `movetype='ALIGN'`,
built by `MakeMembraneSystemTask._orientation_align`.  No pestifer workflow loads any of this.

They are retained solely as the **oracle** for the regression tests in
`../../test_make_membrane_system.py`, which run the old Tcl path and the new Python path over
the same structure and compare coordinates.  That comparison is the evidence that the
coordinate-transform port preserved behavior; deleting these files would delete the proof.

Because they no longer live under `${PESTIFER_TCLROOT}/pkg/`, `vmdrc.tcl` does not glob them
onto `auto_path`.  A test that uses the oracle must therefore do two things: copy
`bilayer_orient.tcl` into its working directory and ingest it with `usescript(..., local=True)`,
and emit `lappend auto_path <this dir>/pkg` into the script *before* that ingest so the
script's `package require Orient` resolves.  See `_stage_orient_oracle` in the test module.
