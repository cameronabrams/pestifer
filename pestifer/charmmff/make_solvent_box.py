# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Build a pre-equilibrated periodic **solvent box** for a single CHARMM residue, for use
as a ``kind: box`` entry in the PDB repository's ``solvent`` collection.  Such a box is
what VMD's ``solvate`` plugin needs to bulk-solvate with anything other than its built-in
TIP3P water (``-spsf``/``-spdb``/``-ws``/``-ks``); see
``docs/design/solvent-collection.md``.

This module holds the deterministic, dependency-light core -- box geometry and cubic
packing -- separately from the psfgen/NAMD build orchestration, so the geometry can be
unit-tested without a force field or MD engine.
"""
import logging
import os
import shutil

import numpy as np

logger = logging.getLogger(__name__)

# Avogadro's number (1/mol); 1 cm^3 = 1e24 A^3
_N_AVOGADRO = 6.02214076e23
_CM3_PER_A3 = 1.0e-24


def box_edge_for_density(nmol: int, molweight: float, density_gcc: float) -> float:
    """
    Edge length (Å) of a cube holding ``nmol`` molecules of molar mass ``molweight``
    (g/mol) at bulk density ``density_gcc`` (g/cm^3).

    Volume = (nmol * molweight / N_A) / density, converted from cm^3 to Å^3.  This sizes
    the *initial* box; the NPT equilibration relaxes it to the true equilibrium edge,
    which is what gets recorded for ``solvate -ws``.
    """
    if nmol <= 0 or molweight <= 0 or density_gcc <= 0:
        raise ValueError('nmol, molweight, and density_gcc must all be positive')
    mass_g = nmol * molweight / _N_AVOGADRO
    volume_cm3 = mass_g / density_gcc
    volume_a3 = volume_cm3 / _CM3_PER_A3
    return volume_a3 ** (1.0 / 3.0)


def _random_rotations(n: int, rng: np.random.Generator) -> np.ndarray:
    """Return ``n`` uniformly-random 3x3 rotation matrices (via random unit quaternions)."""
    # Marsaglia/Shoemake: uniform quaternions -> rotation matrices, no scipy dependency
    u1, u2, u3 = rng.random(n), rng.random(n), rng.random(n)
    q = np.empty((n, 4))
    q[:, 0] = np.sqrt(1 - u1) * np.sin(2 * np.pi * u2)
    q[:, 1] = np.sqrt(1 - u1) * np.cos(2 * np.pi * u2)
    q[:, 2] = np.sqrt(u1) * np.sin(2 * np.pi * u3)
    q[:, 3] = np.sqrt(u1) * np.cos(2 * np.pi * u3)
    w, x, y, z = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    R = np.empty((n, 3, 3))
    R[:, 0, 0] = 1 - 2 * (y * y + z * z); R[:, 0, 1] = 2 * (x * y - z * w); R[:, 0, 2] = 2 * (x * z + y * w)
    R[:, 1, 0] = 2 * (x * y + z * w); R[:, 1, 1] = 1 - 2 * (x * x + z * z); R[:, 1, 2] = 2 * (y * z - x * w)
    R[:, 2, 0] = 2 * (x * z - y * w); R[:, 2, 1] = 2 * (y * z + x * w); R[:, 2, 2] = 1 - 2 * (x * x + y * y)
    return R


def cubic_lattice_sites(nmol: int, edge: float) -> np.ndarray:
    """
    ``nmol`` cell-centered sites of the smallest cubic lattice that holds them, spanning
    an ``edge``-length cube.  Returns an ``(nmol, 3)`` array of centers in ``[0, edge)``.
    """
    ncell = int(np.ceil(round(nmol ** (1.0 / 3.0), 6)))
    while ncell ** 3 < nmol:
        ncell += 1
    spacing = edge / ncell
    grid = (np.arange(ncell) + 0.5) * spacing
    sites = np.array([(x, y, z) for x in grid for y in grid for z in grid])
    return sites[:nmol]


def pack_cubic(mol_coords: np.ndarray, nmol: int, edge: float, seed=None,
               jitter: float = 0.0) -> np.ndarray:
    """
    Place ``nmol`` randomly-oriented copies of a molecule (whose atom coordinates are
    ``mol_coords``, an ``(natom, 3)`` array) on the cubic lattice of
    :func:`cubic_lattice_sites`.

    Each copy is centered on its atoms' centroid, given a uniformly-random orientation,
    and translated to a lattice site (plus optional ``jitter``).  Returns an
    ``(nmol, natom, 3)`` array of placed coordinates.  Overlaps left for the downstream
    NPT relaxation to resolve.
    """
    rng = np.random.default_rng(seed)
    centered = mol_coords - mol_coords.mean(axis=0)
    sites = cubic_lattice_sites(nmol, edge)
    rots = _random_rotations(nmol, rng)
    placed = np.einsum('nij,aj->nai', rots, centered)   # rotate each copy
    if jitter:
        sites = sites + rng.uniform(-jitter, jitter, size=sites.shape)
    placed += sites[:, None, :]
    return placed


def _segid(i: int) -> str:
    """A 4-character base-36 segid for the i-th molecule (covers 36^4 = 1.6M copies)."""
    digits = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    s = ''
    for _ in range(4):
        i, r = divmod(i, 36)
        s = digits[r] + s
    return s


def write_box_pdb(placed: np.ndarray, template_lines: list, output_pdb: str,
                  max_res_per_seg: int = 9999) -> list:
    """
    Write an ``nmol``-molecule box PDB by rewriting one molecule's ATOM lines
    (``template_lines``, ``natom`` fixed-format PDB records) once per placed copy.

    The molecules are grouped into segments of at most ``max_res_per_seg`` residues (so a
    solvent box is a handful of segments of many residues each, like water in a membrane
    build -- not one segment per molecule, which would exhaust chain IDs).  Each molecule
    becomes one residue (``resSeq`` 1..k within its segment) with a shared ``segID``
    (``W0``, ``W1``, ...).  ``placed`` is the ``(nmol, natom, 3)`` array from
    :func:`pack_cubic`.  Returns a list of ``(segid, nresidues)`` for the segments written,
    in order.
    """
    nmol, natom, _ = placed.shape
    if len(template_lines) != natom:
        raise ValueError(f'template has {len(template_lines)} atom lines but placed has {natom} atoms')
    tmpl = [ln.rstrip('\n').ljust(80) for ln in template_lines]
    segments = []   # (segid, nres)
    serial = 0
    _chain_pool = ([chr(c) for c in range(ord('A'), ord('Z') + 1)]
                   + [chr(c) for c in range(ord('a'), ord('z') + 1)]
                   + [chr(c) for c in range(ord('0'), ord('9') + 1)])
    with open(output_pdb, 'w') as f:
        for i in range(nmol):
            seg_index, res_in_seg = divmod(i, max_res_per_seg)
            sid = f'W{seg_index}'
            chain = _chain_pool[seg_index % len(_chain_pool)]
            resid = res_in_seg + 1
            if resid == 1:
                segments.append([sid, 0])
            segments[-1][1] += 1
            for a in range(natom):
                serial += 1
                x, y, z = placed[i, a]
                ln = tmpl[a]
                rec = (ln[:6]
                       + f'{serial % 100000:5d}'
                       + ln[11:21]
                       + chain                          # chainID must be in the manager's pool
                       + f'{resid:4d}'
                       + ln[26:30]
                       + f'{x:8.3f}{y:8.3f}{z:8.3f}'
                       + ln[54:72]
                       + f'{sid:<4s}'
                       + ln[76:80])
                f.write(rec.rstrip() + '\n')
        f.write('END\n')
    return [(sid, n) for sid, n in segments]


def atom_lines_of(pdb_path: str) -> list:
    """Return the ATOM/HETATM record lines of a PDB file (in order)."""
    with open(pdb_path) as f:
        return [ln for ln in f if ln.startswith(('ATOM', 'HETATM'))]


def coords_of(atom_lines: list) -> np.ndarray:
    """Parse the (natom, 3) coordinate array from fixed-format PDB ATOM lines."""
    return np.array([[float(ln[30:38]), float(ln[38:46]), float(ln[46:54])] for ln in atom_lines])


def _resi_to_molblock(topo) -> str:
    """
    Serialize a CHARMM RESI's atom+bond graph to an MDL V2000 molblock with zeroed
    coordinates, atoms in RESI order.  This is the connectivity-only input Open Babel's
    ``--gen3d`` expands into 3D coordinates; keeping RESI atom order means the generated
    coordinates map 1:1 back to the CHARMM atom names with no re-matching.
    """
    atoms = topo.atoms
    bonds = topo.bonds
    name_to_idx = {a.name: i + 1 for i, a in enumerate(atoms)}
    # per-atom valence = sum of its bond orders.  Pinning it in the molblock's ``vvv`` field
    # stops obabel from adding implicit hydrogens to atoms it perceives as under-valent --
    # e.g. a sulfoxide S=O that CHARMM lists as a single bond (DMSO would otherwise gain 2 H).
    valence = {a.name: 0 for a in atoms}
    for b in bonds:
        valence[b.name1] += int(b.degree)
        valence[b.name2] += int(b.degree)
    lines = [topo.resname, '  pestifer', '']
    lines.append(f'{len(atoms):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000')
    for a in atoms:
        v = valence[a.name]
        # V2000 atom line: coords, symbol, then dd ccc sss hhh bbb vvv HHH rrr iii mmm nnn eee
        lines.append(f'{0.0:10.4f}{0.0:10.4f}{0.0:10.4f} {a.element:<3s}'
                     f'{0:2d}{0:3d}{0:3d}{0:3d}{0:3d}{v:3d}{0:3d}{0:3d}{0:3d}{0:3d}{0:3d}{0:3d}')
    for b in bonds:
        lines.append(f'{name_to_idx[b.name1]:3d}{name_to_idx[b.name2]:3d}{int(b.degree):3d}  0  0  0  0')
    lines.append('M  END')
    return '\n'.join(lines) + '\n'


def single_molecule_from_graph(topo, output_pdb: str) -> str:
    """
    Build a single-molecule PDB for ``topo`` (a CHARMM RESI with no usable internal
    coordinates) by generating 3D coordinates from its atom+bond graph with Open Babel's
    ``--gen3d``, then writing them out under the RESI's own atom names.  Returns
    ``output_pdb``.

    The rough ``gen3d`` geometry only has to seed psfgen/``coordpdb``; the box is minimized
    and NPT-equilibrated afterward, so small imperfections wash out.  Requires ``obabel`` on
    PATH (already a pestifer dependency via the ligand-parametrization path).
    """
    from ..core.command import Command

    resname = topo.resname
    molblock = _resi_to_molblock(topo)
    molpath = f'{resname}-seed.mol'
    xyzpath = f'{resname}-seed.xyz'
    with open(molpath, 'w') as f:
        f.write(molblock)
    c = Command(f'obabel {molpath} -O {xyzpath} --gen3d')
    rc = c.run(quiet=True)
    if rc != 0 or not os.path.exists(xyzpath):
        raise RuntimeError(
            f'obabel --gen3d failed to build coordinates for {resname} (rc={rc}); '
            f'{c.stderr or c.stdout}')

    with open(xyzpath) as f:
        xyz_lines = f.read().splitlines()
    # xyz: line 0 = atom count, line 1 = comment, then "<element> x y z" in molblock order
    coords = []
    for ln in xyz_lines[2:]:
        parts = ln.split()
        if len(parts) >= 4:
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
    if len(coords) != len(topo.atoms):
        raise RuntimeError(
            f'obabel returned {len(coords)} atoms for {resname} but the RESI has '
            f'{len(topo.atoms)}; cannot map coordinates to atom names')

    with open(output_pdb, 'w') as f:
        for i, (a, (x, y, z)) in enumerate(zip(topo.atoms, coords), start=1):
            name = a.name
            # PDB fixed columns: atom name 13-16 (<4-char names start in col 14),
            # altLoc 17, resName 18-21, chainID 22 (blank), resSeq 23-26
            namecol = f'{name:<4s}' if len(name) >= 4 else f' {name:<3s}'
            f.write(f'ATOM  {i:5d} {namecol} {resname:<4s} {1:4d}    '
                    f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00      A   {a.element:>2s}\n')
        f.write('END\n')
    return output_pdb


def make_solvent_box(resname, DB, RM=None, nmol: int = 216, density: float = 1.0,
                     temperature: float = 300.0, pressure: float = 1.0,
                     minimize_steps: int = 1000, npt_steps: int = 50000,
                     seed: int = None, key_atom: str = None, refic_idx: int = 0) -> dict:
    """
    Build a pre-equilibrated periodic solvent box for ``resname`` and return a summary.

    Pipeline: build one molecule from the residue's internal coordinates (psfgen) ->
    pack ``nmol`` randomly-oriented copies into a cube sized for ``density`` -> psfgen the
    box (one segment per molecule) -> minimize + NPT-equilibrate under PBC -> measure the
    equilibrated edge.  Leaves ``<resname>-box.psf``/``.pdb`` (the equilibrated box) and an
    ``info.yaml`` (``kind: box``) in the cwd; the caller installs the entry into the
    ``solvent`` collection.
    """
    import yaml

    from ..core.config import Config
    from ..core.controller import Controller
    from ..scripters.psfgen import PsfgenScripter
    from ..scripters.namd import NAMDScripter
    from ..psfutil.psfcontents import PSFContents
    from ..util.util import cell_to_xsc, cell_from_xsc

    topo = DB.get_resi(resname)
    meta = topo.metadata
    charmm_topfile = os.path.basename(meta['charmmfftopfile'])
    extension = os.path.splitext(charmm_topfile)[1][1:]

    # Controller/config gives us the psfgen + namd scripters and drives the equilibration
    tasklist = [
        {'continuation': {'psf': f'{resname}-box.psf', 'pdb': f'{resname}-box.pdb', 'xsc': f'{resname}-box.xsc'}},
        {'md': {'ensemble': 'minimize', 'nsteps': 0, 'minimize': minimize_steps,
                'dcdfreq': 0, 'xstfreq': 0, 'temperature': temperature}},
        {'md': {'ensemble': 'NPT', 'nsteps': npt_steps, 'dcdfreq': max(npt_steps, 1),
                'xstfreq': 1000, 'temperature': temperature, 'pressure': pressure}},
    ]
    config = Config(quiet=True, RM=RM).configure_new()
    C = Controller().configure(config=config,
                               userspecs={'title': f'Solvent box for {resname}', 'tasks': tasklist},
                               terminate=False)

    W: PsfgenScripter = C.tasks[0].scripters['psfgen']
    if extension in ('rtf', 'str') and charmm_topfile not in W.charmmff_config['standard'][extension]:
        W.charmmff_config['standard'][extension].append(charmm_topfile)

    # 1. single molecule: coordinates come from, in order of preference,
    #    (a) an existing PDB-repository entry (e.g. water/ions have no ICs to build from),
    #    (b) the RESI's internal coordinates via psfgen, or
    #    (c) Open Babel's --gen3d applied to the RESI's atom+bond graph (for residues that
    #        ship neither a PDB entry nor a usable IC table, e.g. many CGenFF solvents).
    if DB.pdbrepository is None:
        DB.provision_pdbrepository()

    def _build_single_from_pdb(pdb_path):
        W.newscript('single', additional_topologies=[charmm_topfile])
        W.addline(f'segment A {{ first NONE; last NONE; residue 1 {resname} }}')
        W.addline(f'coordpdb {pdb_path} A')
        W.writescript(f'{resname}-single')
        W.runscript()

    if resname in DB.pdbrepository:
        _build_single_from_pdb(DB.pdbrepository.checkout(resname).get_pdb(conformerID=0))
    else:
        W.newscript('single', additional_topologies=[charmm_topfile])
        if topo.to_psfgen(W, refic_idx=refic_idx) != -1:
            W.writescript(f'{resname}-single')
            W.runscript()
        else:
            # no usable internal coordinates: build a seed molecule from the RESI's
            # atom+bond graph with obabel --gen3d (fresh script; the failed one is discarded)
            logger.info(f'{resname} has no usable internal coordinates; '
                        f'generating a seed molecule with obabel --gen3d')
            _build_single_from_pdb(single_molecule_from_graph(topo, f'{resname}-seed.pdb'))

    # 2. molar mass + geometry, then pack
    psf = PSFContents(f'{resname}-single.psf')
    molweight = float(sum(a.atomicwt for a in psf.atoms))
    atomlines = atom_lines_of(f'{resname}-single.pdb')
    single_coords = coords_of(atomlines)
    edge = box_edge_for_density(nmol, molweight, density)
    logger.info(f'{resname} box: {nmol} molecules, MW {molweight:.2f}, target {density} g/cc -> initial edge {edge:.2f} A')
    placed = pack_cubic(single_coords, nmol, edge, seed=seed)
    segments = write_box_pdb(placed, atomlines, f'{resname}-boxin.pdb')

    # 3. box psfgen: a few segments of many residues each (like water in a membrane build)
    W.newscript('box', additional_topologies=[charmm_topfile])
    for sid, nres in segments:
        W.addline(f'segment {sid} {{')
        W.addline('  first NONE')
        W.addline('  last NONE')
        for r in range(1, nres + 1):
            W.addline(f'  residue {r} {resname}')
        W.addline('}')
        W.addline(f'coordpdb {resname}-boxin.pdb {sid}')
    W.writescript(f'{resname}-box', guesscoord=True, regenerate=True)
    W.runscript()

    # 4. initial cubic cell centered on the packed box ([0, edge) per axis)
    box = np.array([[edge, 0.0, 0.0], [0.0, edge, 0.0], [0.0, 0.0, edge]])
    origin = np.array([edge / 2.0, edge / 2.0, edge / 2.0])
    cell_to_xsc(box, origin, f'{resname}-box.xsc')

    # 5. minimize + NPT-equilibrate under PBC
    na: NAMDScripter = C.tasks[0].scripters['namd']
    if extension in ('rtf', 'str') and charmm_topfile not in na.charmmff_config['standard'][extension]:
        na.charmmff_config['standard'][extension].append(charmm_topfile)
    result = C.do_tasks()
    for k, v in result.items():
        if v['result'] != 0:
            raise RuntimeError(f'solvent-box equilibration task {v["taskname"]} failed (result {v["result"]})')

    # 6. promote the equilibrated NPT output as the box coordinates.  NAMD's wrapAll leaves
    # every molecule whole and periodic-consistent with the equilibrated cell, which is
    # exactly what VMD's ``solvate -ws <edge>`` needs (atoms may poke just past an edge --
    # those are the periodic images and tile cleanly).  The initial packed pdb is discarded.
    final_task = C.tasks[-1]
    final_pdb = f'{final_task.basename}.pdb'
    final_xsc = f'{final_task.basename}.xsc'
    shutil.copyfile(final_pdb, f'{resname}-box.pdb')

    # read the equilibrated edge from the final NPT xsc
    final_box, _ = cell_from_xsc(final_xsc)
    final_edge = float(final_box[0][0])
    final_density = nmol * molweight / (_N_AVOGADRO * (final_edge ** 3) * _CM3_PER_A3)
    if key_atom is None:
        key_atom = psf.atoms[0].atomname
    logger.info(f'{resname} box equilibrated: edge {final_edge:.3f} A, density {final_density:.3f} g/cc')

    info = {
        'kind': 'box', 'resname': resname, 'nmol': nmol,
        'box_edge': float(f'{final_edge:.4f}'), 'density': float(f'{final_density:.4f}'),
        'key_atom': key_atom, 'defined-in': charmm_topfile,
        'parameters': sorted(set(na.charmmff_config['standard']['prm'] + na.charmmff_config['standard']['str'])),
        'psf': f'{resname}-box.psf', 'pdb': f'{resname}-box.pdb',
    }
    with open('info.yaml', 'w') as f:
        yaml.dump(info, f)
    return info


def build_solvent_entry(args):
    """
    Driver for ``make-pdb-collection solvent``: build a ``kind: box`` solvent entry and
    assemble it into an installable entry directory ``<output_dir>/<RESN>/`` containing
    ``info.yaml`` + the equilibrated box psf/pdb.  Install it with
    ``modify-package pdb-repo add-entry <output_dir>/<RESN> --collection solvent``.

    The heavy NAMD build runs in a scratch working directory so its many intermediate
    files stay out of the entry directory.
    """
    from ..core.resourcemanager import ResourceManager

    resname = args.resname
    charmmff_config = {'release': args.charmmff_release} if getattr(args, 'charmmff_release', '') else {}
    RM = ResourceManager(charmmff_config=charmmff_config)
    CC = RM.charmmff_content
    CC.provision()
    if resname not in CC:
        raise ValueError(f'RESI {resname} not found in the CHARMM force field; cannot build a solvent box')

    outdir = args.output_dir or 'solvent'
    entrydir = os.path.abspath(os.path.join(outdir, resname))
    workdir = os.path.abspath(os.path.join(outdir, f'{resname}-work'))
    for d in (entrydir, workdir):
        os.makedirs(d, exist_ok=True)

    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        info = make_solvent_box(
            resname, CC, RM=RM, nmol=args.nmol, density=args.density,
            temperature=args.temperature, pressure=args.pressure,
            minimize_steps=args.minimize_steps, npt_steps=args.npt_steps,
            seed=args.seed, key_atom=(args.key_atom or None), refic_idx=args.refic_idx)
        for fname in ('info.yaml', info['psf'], info['pdb']):
            shutil.copyfile(fname, os.path.join(entrydir, fname))
    finally:
        os.chdir(cwd)
    if getattr(args, 'cleanup', True):
        shutil.rmtree(workdir, ignore_errors=True)

    logger.info(f"solvent box entry for {resname} written to {os.path.relpath(entrydir, cwd)}/ "
                f"(edge {info['box_edge']} A, density {info['density']} g/cc)")
    print(f"Solvent box entry for {resname}: {os.path.relpath(entrydir, cwd)}/")
    print(f"  edge {info['box_edge']} A, density {info['density']} g/cc, key atom {info['key_atom']}")
    print("Install it into the solvent collection with:")
    print(f"    pestifer modify-package pdb-repo add-entry {os.path.relpath(entrydir, cwd)} --collection solvent")
    return info
