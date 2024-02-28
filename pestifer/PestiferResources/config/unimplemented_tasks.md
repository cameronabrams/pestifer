packmol-memgen
--------------

pestifer should

1. orient the protein in the box (n_ter in or out, COM z at 0, or some other way to locate embedding orientation)
2. call packmol-memgen with notrun and preoriented
3. modify the packmol script to get rid of spurious commands and set positions correctly
4. run packmol

Note: Two instances of `np.float` changed to `np.float64` in packmol_memgen/lib/pdbremix/v3numpy.py

