# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
This module implements the ``make-pdb-collection`` subcommand.  It generates a PDB collection from the CHARMM force field residue database.  The purpose is to generate 3D structures for all standalone (small-molecule) residues defined in the CHARMM force field, which can be used for molecular dynamics simulations and system preparation.
"""
from copy import deepcopy
import glob
import logging
import os
import shutil
import yaml
import numpy as np
from   itertools            import product

from   pidibble.pdbparse    import PDBParser

from pestifer.objs import resid
from pestifer.scripters.namd import NAMDScripter
from pestifer.scripters.psfgen import PsfgenScripter

from  .charmmffcontent      import CHARMMFFContent
from ..core.config          import Config
from ..core.controller      import Controller
from ..core.resourcemanager import ResourceManager
from ..util.stringthings    import my_logger
from ..psfutil.psfcontents  import PSFContents

logger = logging.getLogger(__name__)

# Fluid area per acyl chain (A^2).  A single-molecule athermal-MC confinement cylinder is sized as
# n_chains * this value (times the geometric inflation factor).  Calibrated from DMPC: a 2-chain lipid
# gets cylinder_apl = 60 (2 x 30), which lands the melted ensemble near the fluid DMPC APL ~67.  This
# is what lets a 4-chain cardiolipin be confined at ~120 rather than crushed into a 2-chain cylinder.
_PER_CHAIN_APL = 30.0

# Target chain order parameter (the -mean 1/2(3cos^2 theta - 1) over tail C-H vs the membrane normal
# computed by athermal_mc.ensemble_chain_order) for each bilayer phase.  Ld = None: the pure athermal
# MC (trans bias 0) already IS the fluid ensemble -- its order is the disordered floor, so there is
# nothing to tune (you cannot go more disordered).  Lo tunes the run_mc trans bias UP to an ordered
# target, extending the chains toward the liquid-ordered (cholesterol-rich raft) state.  The Lo value
# is PROVISIONAL -- calibrate against reference bilayers / the ex17 acceptance test; see
# docs/design/lipid-conformer-generation.md (Phase section).
_PHASE_ORDER_TARGET = {'Ld': None, 'Lo': 0.30}


def phase_order_target(phase):
    """Chain-order target for a bilayer ``phase`` (``'Ld'`` | ``'Lo'``), or ``None`` for no tuning.

    ``None`` -- returned for no phase and for ``'Ld'`` (the athermal fluid ensemble needs no ordering
    bias) -- means "do not tune": the generator samples once at trans bias 0.
    """
    if phase is None:
        return None
    try:
        return _PHASE_ORDER_TARGET[phase]
    except KeyError:
        raise ValueError(f'unknown bilayer phase {phase!r}; expected one of '
                         f'{sorted(_PHASE_ORDER_TARGET)}')


def _bisect_to_order(order_of, target, bounds, tol=0.02, max_iters=6):
    """Bisect a scalar knob to hit a chain-order ``target`` for a monotonic ``order_of(knob)``.

    The monotonic direction (order rising or falling with the knob) is inferred from the endpoints,
    so this serves either parametrization -- e.g. a trans-bias strength (order rises with it) or a
    cylinder looseness (order falls).  Returns ``(best_knob, best_order, bracketed)``; ``bracketed``
    is False when ``target`` lies outside the order range spanned by ``bounds`` (then the closest
    bound is returned and the caller should warn).
    """
    a, b = bounds
    oa, ob = order_of(a), order_of(b)
    best = min(((abs(oa - target), a, oa), (abs(ob - target), b, ob)), key=lambda t: t[0])
    if target < min(oa, ob) - 1e-9 or target > max(oa, ob) + 1e-9:
        return best[1], best[2], False
    increasing = ob >= oa                      # order rises from a to b
    for _ in range(max_iters):
        if best[0] <= tol:
            break
        m = 0.5 * (a + b)
        om = order_of(m)
        if abs(om - target) < best[0]:
            best = (abs(om - target), m, om)
        if (om < target) == increasing:        # target lies on the b side of m
            a = m
        else:                                  # target lies on the a side of m
            b = m
    return best[1], best[2], True


def auto_cylinder_apl(chain_info, *, substream: str = '', per_chain: float = _PER_CHAIN_APL) -> float:
    """On-the-fly confinement-cylinder APL from the lipid's acyl-chain count.

    ``chain_info`` is the per-chain list produced by the head/tail network analysis
    (:meth:`~pestifer.charmmff.charmmfftop.CharmmResiTopology.generic_lipid_annotate`), so
    ``len(chain_info)`` is the number of acyl chains: 2 for a di-acyl PL, 4 for a cardiolipin, 1 for a
    lyso lipid.  The cylinder is sized ``n_chains * per_chain`` (A^2) so its cross-section tracks the
    lipid's real footprint rather than a global default.  Detergents annotate with an empty
    ``chain_info`` but do carry one chain, so 0 chains falls back to 1 (never 0)."""
    n = len(chain_info) if chain_info else 0
    if n == 0:
        n = 1   # detergents (and any lipid the labeler left un-split) have at least one chain
    return float(n) * per_chain


def _sample_and_write_mc_conformers(resid: str, psf_file: str, pdb_file: str,
                                    heads: list, tails: list, par_basenames: list,
                                    DB: CHARMMFFContent, *, nsamples: int, digits: int,
                                    cylinder_apl: float, cylinder_inflation: float,
                                    n_equil: int, n_decorr: int, seed: int,
                                    max_angle: float, radius_scale: float,
                                    target_order: float = None, order_tol: float = 0.02,
                                    max_tune_iters: int = 6,
                                    bias_bounds=(0.0, 40.0)) -> dict:
    """Athermal-MC sampler: melt the tails of a built lipid into a fluid-like conformer set.

    Given the built, minimized single-molecule ``psf_file`` + ``pdb_file`` (the same artifacts
    the vacuum-MD sampler would relax), identify the acyl-tail torsions, run cylinder-confined
    athermal MC (:mod:`pestifer.charmmff.athermal_mc`), and write ``{resid}-NN.pdb`` conformers
    that the shared ``info.yaml`` block below then measures.  This is the drop-in replacement for
    the vacuum ``NVT`` sampling that produced over-thin rod conformers.

    Parameters
    ----------
    psf_file, pdb_file : str
        The built topology and (minimized/oriented) coordinates.
    heads, tails : list[str]
        Head reference atom name(s) and terminal tail-carbon name(s), as computed for the lipid.
    par_basenames : list[str]
        CHARMM parameter/stream file basenames whose NONBONDED records supply per-atom ``Rmin/2``.
    DB : CHARMMFFContent
        Used to resolve ``par_basenames`` to absolute paths.
    nsamples, digits, cylinder_apl, n_equil, n_decorr, seed
        MC controls; ``digits`` zero-pads the conformer filenames to match ``do_psfgen``.

    Returns
    -------
    int
        Number of conformer PDBs written.
    """
    from .charmmffprm import CharmmParamFile
    from .athermal_mc import (build_lipid_mc, run_mc, cylinder_radius_for_apl,
                              ensemble_chain_order)

    psf = PSFContents(psf_file, parse_topology=['bonds'])
    name_of = {a.serial: a.atomname for a in psf.atoms.data}
    type_of = {a.serial: a.atomtype for a in psf.atoms.data}
    mass_of = {a.serial: a.atomicwt for a in psf.atoms.data}

    # read the template PDB in file order; MC coordinate rows are indexed the same way
    template, coords, serials = [], [], []
    with open(pdb_file) as fh:
        for ln in fh:
            if ln.startswith(('ATOM', 'HETATM')):
                serials.append(int(ln[6:11]))
                coords.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
                template.append(ln)
    coords = np.array(coords)
    row_of = {s: i for i, s in enumerate(serials)}

    masses = [mass_of[s] for s in serials]
    atomtypes = [type_of[s] for s in serials]
    bonds = [(row_of[b.serial1], row_of[b.serial2]) for b in psf.bonds]
    head_idx = [i for i, s in enumerate(serials) if name_of[s] in heads]
    tail_idx = [i for i, s in enumerate(serials) if name_of[s] in tails]

    params = CharmmParamFile()
    for base in set(par_basenames):
        try:
            local = DB.copy_charmmfile_local(base)   # materialize into cwd, return basename
            params.merge(CharmmParamFile.from_file(local))
        except Exception as exc:
            logger.debug(f'MC: skipping parameter file {base}: {exc}')
    missing = sorted(t for t in set(atomtypes) if t not in params.nonbonded)
    if missing:
        raise RuntimeError(f'MC conformer generation for {resid}: missing vdW (Rmin/2) '
                           f'parameters for atom types {missing}')
    rmin_half = [params.nonbonded[t].Rmin_half for t in atomtypes]

    # The cylinder (fixed at cylinder_inflation) is the *packing* knob -- it sets the footprint/APL,
    # not the chain order (an athermal MC stays fluid however tight the cylinder, because nothing
    # drives the chains trans).  Ordering is a separate knob: the run_mc trans bias.  So build the
    # molecule ONCE and vary only the bias per probe.
    R = cylinder_radius_for_apl(cylinder_apl, cylinder_inflation)
    mol = build_lipid_mc(coords, None, masses, bonds, head_idx, tail_idx, rmin_half,
                         cylinder_radius=R, radius_scale=radius_scale)
    if not mol.rotatable:
        logger.warning(f'MC {resid}: no rotatable tail torsions found; conformers will be copies')

    probed = {}   # bias -> (samples, chain_order)

    def order_of(bias):
        samples = run_mc(mol, nsamples=nsamples, n_equil=n_equil, n_decorr=n_decorr,
                         max_angle=max_angle, seed=seed, torsion_bias=bias)
        order = ensemble_chain_order(samples, None, masses, bonds)
        probed[bias] = (samples, order)
        logger.info(f'MC {resid}: probe trans bias {bias:.2f} -> chain order {order:.3f}')
        return order

    if target_order is None:
        used_bias = 0.0
        order_of(used_bias)
        samples, order = probed[used_bias]
        logger.info(f'MC {resid}: {len(mol.rotatable)} tail torsions, '
                    f'{int(mol.confined.sum())} confined atoms, cylinder radius {R:.2f} A '
                    f'(APL {cylinder_apl:.0f} x inflation {cylinder_inflation:.2f}); '
                    f'chain order {order:.3f} (no bias)')
    else:
        used_bias, tuned_order, bracketed = _bisect_to_order(
            order_of, target_order, bounds=bias_bounds, tol=order_tol, max_iters=max_tune_iters)
        samples, order = probed[used_bias]
        if not bracketed:
            logger.warning(f'MC {resid}: chain-order target {target_order:.3f} not reachable within '
                           f'trans-bias {bias_bounds}; using the closest bound '
                           f'(bias {used_bias:.2f}, order {order:.3f})')
        logger.info(f'MC {resid}: tuned trans bias to {used_bias:.2f} at fixed cylinder '
                    f'(radius {R:.2f} A, APL {cylinder_apl:.0f} x inflation {cylinder_inflation:.2f}); '
                    f'chain order {order:.3f} vs target {target_order:.3f} '
                    f'({len(mol.rotatable)} tail torsions)')

    for f, s in enumerate(samples):
        with open(f'{resid}-{f:0{digits}d}.pdb', 'w') as out:
            for row, ln in enumerate(template):
                x, y, z = s[row]
                out.write(f'{ln[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{ln[54:]}')
            out.write('END\n')
    return {'nconf': len(samples), 'torsion_bias': float(used_bias),
            'chain_order': (float(order) if order == order else None)}


def do_psfgen(resid: str, DB: CHARMMFFContent, RM: ResourceManager = None,
              lenfac: float = 1.2, minimize_steps: int = 500,
              sample_steps: int = 5000, nsamples: int = 10,
              sample_temperature: float = 300,
              refic_idx: int = 0, force_constant: float = 1.0,
              borrow_ic_from: str = None, sampler: str = 'md',
              cylinder_apl: float = None, cylinder_inflation: float = 1.9,
              mc_n_equil: int = 20000, mc_n_decorr: int = 3000, mc_seed: int = 0,
              mc_max_angle: float = np.pi / 6, mc_radius_scale: float = 0.8,
              phase: str = None):
    """ 
    Generate a PDB file for a residue defined by the CHARMM force field using psfgen, and sample it using NAMD.  Also generate the ``info.yaml`` file for this residue.

    This function uses the :class:`pdibble.pdbparse.PDBParser` to read in generated PDB files and :class:`~pestifer.psfutil.psfcontents.PSFContents` to read in generated PSF files to make measurements needed for the ``info.yaml`` files.

    Parameters
    ----------
    resid : str
        The residue ID for which to generate the PDB file.
    DB : :class:`~pestifer.charmmff.charmmffcontent.CHARMMFFContent`
        The database containing CHARMM residue information.
    lenfac : float, optional
        The lengthening factor to apply to the residue (default is 1.2).
    minimize_steps : int, optional
        The number of minimization steps to perform (default is 500).
    sample_steps : int, optional
        The total number of sampling steps to perform (default is 5000).
    nsamples : int, optional
        The number of samples to generate (default is 10).
    sample_temperature : float, optional
        The temperature to use for sampling (default is 300 K).
    refic_idx : int, optional
        The index of the reference IC to use (default is 0).
    force_constant : float, optional
        The force constant to use for sampling (default is 1.0).
    borrow_ic_from : str, optional
        The residue ID from which to borrow internal coordinates (default is None).

    """
    if nsamples > sample_steps:
        raise ValueError(f'nsamples ({nsamples}) must be less than or equal to sample_steps ({sample_steps})')
    digits = len(str(nsamples))
    topo = DB.get_resi(resid)
    synonym = topo.synonym
    meta = topo.metadata
    charmm_topfile = meta.get('charmmtopfile', None)
    logger.debug(f'do_psfgen for {resid} with metadata {meta}')
    if borrow_ic_from:
        logger.debug(f'borrowing ICs from {borrow_ic_from}')
        take_my_ics = DB.get_resi(borrow_ic_from)
        topo.copy_ICs_from(take_my_ics)
    charmm_topfile, stream, substream = meta['charmmfftopfile'], meta['streamID'], meta['substreamID']
    charmm_topfile = os.path.basename(charmm_topfile)
    logger.debug(f'charmm_topfile: {charmm_topfile}, stream: {stream}, substream: {substream}')

    # Sterols (cholesterol substream) are dominated by a rigid fused-ring system; the grid packer
    # only ever needs a single representative conformer for them.  Force single-conformer generation
    # -- skip both the MD and the athermal-MC tail sampling and emit one minimized structure --
    # regardless of the requested sampler.
    is_sterol = (substream == 'cholesterol')
    if is_sterol and sampler != 'single':
        logger.info(f'{resid} is a sterol ({substream}); generating a single rigid conformer '
                    f"(skipping '{sampler}' sampling)")
        sampler = 'single'
    if sampler == 'single':
        nsamples, digits = 1, 1
    # substream_overrides = DB.overrides.get('substreams', {})
    # ss_override = substream_overrides.get(resid, None)
    # if ss_override is not None:
    #     substream = ss_override
    #     logger.debug(f'Overriding substream for {resid} from {meta["substreamID"]} to {substream}')
    #     meta['substreamID'] = substream
    if charmm_topfile is None:
        return -2
    if topo.error_code != 0:
        topo.show_error(logger.warning)
        # my_logger(f'Parsing error for {resid}',logger.warning,just='^',frame='!',fill='!') 
        with open(f'{resid}-topo.rtf','w') as f:
            topo.to_file(f)
        # return -1

    # for b in topo.bonds:
    #     logger.debug(f'{b.name1}-{b.name2}')
    my_logger(f'{os.path.basename(charmm_topfile)}', logger.info, just='^', frame='*', fill='*')

    tasklist = [
        {'continuation': {'psf': f'{resid}-init.psf', 'pdb': f'{resid}-init.pdb'}},
         {'md': {'ensemble': 'minimize', 'nsteps': 0, 'minimize': minimize_steps, 'dcdfreq': 0, 'xstfreq': 0, 'temperature': 100}},
    ]
    heads = []
    chain_info = []
    if stream == 'lipid':
        topo.lipid_annotate()
        ano = topo.annotation
        heads, tails, shortest_paths, chain_info = ano.get('heads', []), ano.get('tails', []), ano.get('shortest_paths', {}), ano.get('chain_info', [])
        if substream == 'cholesterol' or substream == 'model':
            pass
        elif substream == 'detergent':
            logger.debug('detergent')
            logger.debug(f'heads: {heads}. tails: {tails}, shortest_paths: {shortest_paths}, chain_info: {chain_info}')
            pass
        else:
            logger.debug(f'heads: {heads}. tails: {tails}, shortest_paths: {shortest_paths}, chain_info: {chain_info}')
            groups = {}
            distances = {}
            dvalue = {}
            harmonics = {}
            colvars = []
            distance = []
            for chain in chain_info:
                logger.debug(f'chain: {chain}')
                head = chain.get('head', None)
                heads.append(head)
                tail = chain.get('tail', None)
                tails.append(tail)
                groups[f'head_{head}'] = {'atomnames': [head]}
                groups[f'tail_{tail}'] = {'atomnames': [tail]}
                distances[f'head_{head}_tail_{tail}'] = {'groups': [f'head_{head}', f'tail_{tail}']}
                dvalue[f'head_{head}_tail_{tail}'] = chain.get('tail_distance',-1) * lenfac
                colvars.append(f'head_{head}_tail_{tail}')
                distance.append(dvalue[f'head_{head}_tail_{tail}'])
                harmonics[f'head_{head}_tail_{tail}_stretch'] = {
                    'colvars': [f'head_{head}_tail_{tail}'],
                    'forceConstant': force_constant,
                    'distance':[dvalue[f'head_{head}_tail_{tail}']]}
            for i_chain in chain_info:
                for j_chain in chain_info:
                    if i_chain == j_chain or j_chain['chain_id'] < i_chain['chain_id']:
                        continue
                    i_tail = i_chain.get('tail', None)
                    j_tail = j_chain.get('tail', None)
                    distances[f'tail_{i_tail}_tail_{j_tail}'] = {'groups': [f'tail_{i_tail}', f'tail_{j_tail}']}
                    harmonics[f'tail_{i_tail}_tail_{j_tail}_attract'] = {
                        'colvars': [f'tail_{i_tail}_tail_{j_tail}'],
                        'forceConstant': force_constant,
                        'distance':[4.0]}
            # run a non-equilibrium MD simulation to bring the tails together to make a 'standard' conformation
            base_md = {'ensemble': 'NVT', 'nsteps': 15000, 'dcdfreq': 100, 'xstfreq': 100, 'temperature': 100}
            colvar_specs = {'groups': groups, 'distances': distances, 'harmonics': harmonics}
            base_md['colvar_specs'] = deepcopy(colvar_specs)
            assert 'minimize' not in base_md
            logger.debug(f'base_md {base_md}')
            # The colvar-driven stretch pushes tails apart into a 'standard' (rod-like) shape --
            # exactly what the athermal-MC sampler must avoid, so skip it for sampler='mc'.
            if sampler != 'mc':
                tasklist.append({'md': base_md})
    else:
        my_logger(f'{resid} is from stream {stream}', logger.debug, just='^', frame='!', fill='!')
        with open(f'{resid}-topo.rtf', 'w') as f:
            topo.to_file(f)
    logger.debug(f'heads: {heads}')
    if len(heads) > 0:
        # reorient molecule so head is at highest z position if molecule is rotated about its COM
        tasklist.append({'manipulate': {'mods': {'orient': [f'z,{heads[0]}']}}})
    # The vacuum sampling MD (below) generates the conformer ensemble for sampler='md'.  For
    # sampler='mc' the ensemble comes from athermal MC on the minimized+oriented structure
    # instead (run after do_tasks); for sampler='single' (sterols) the minimized+oriented
    # structure IS the sole conformer.  Neither schedules a sampling MD.
    if sampler not in ('mc', 'single'):
        # do a conformer-generation MD simulation with the external forces dialed down a bit
        base_md = {'ensemble': 'NVT', 'nsteps': sample_steps, 'dcdfreq': sample_steps // nsamples, 'xstfreq': 100, 'temperature': sample_temperature}
        if substream not in ['cholesterol', 'detergent']:
            for cv, spec in colvar_specs['harmonics'].items():
                if 'forceConstant' in spec:
                    spec['forceConstant'] *= 0.1
            base_md['colvar_specs'] = deepcopy(colvar_specs)
        tasklist.append({'md': base_md})
        logger.debug(f'base_md (2): {base_md}')
    # All NAMD runs for PDB collection are single-molecule vacuum runs; restrict to 1 core.
    for task_dict in tasklist:
        if 'md' in task_dict:
            task_dict['md']['single-core'] = True
    userspecs = {'title': f'Build PDBCollection entry for {resid}', 'tasks': tasklist}
    config = Config(quiet=True, RM=RM).configure_new()
    C = Controller().configure(config=config, userspecs=userspecs, terminate=False)
    # First we de-novo generate a pdb/psf file using seeded internal coordinates
    W: PsfgenScripter = C.tasks[0].scripters['psfgen']
    logger.debug(f'charmm_topfile from resiDB entry {topo.resname}: {charmm_topfile}')
    extension = os.path.splitext(charmm_topfile)[1][1:]
    if extension not in ['rtf', 'str']:
        raise ValueError(f'Invalid CHARMM topology file extension: {extension}. Expected .rtf or .str.')
    if not charmm_topfile in W.charmmff_config['standard'][extension]:
        W.charmmff_config['standard'][extension].append(charmm_topfile)
    if charmm_topfile.endswith('detergent.str'):
        needed='toppar_all36_lipid_sphingo.str'
        W.charmmff_config['standard']['str'].append(needed)
    if charmm_topfile.endswith('initosol.str'):
        needed='toppar_all36_carb_glycolipid.str'
        W.charmmff_config['standard']['str'].append(needed)

    W.newscript('init')
    if topo.to_psfgen(W, refic_idx=refic_idx) == -1:
        my_logger(f'No valid IC\'s for {resid}', logger.warning, just='^', frame='!', fill='!')
        with open(f'{resid}-topo.rtf', 'w') as f:
            topo.to_file(f)
        return -1
    W.writescript(f'{resid}-init')
    W.runscript()
    with open('init.log', 'r') as f:
        loglines = f.read().split('\n')
    for l in loglines:
        if 'Warning: failed to guess coordinates' in l:
            my_logger(f'psfgen failed for {resid}', logger.warning, just='^', frame='!', fill='!')
            with open(f'{resid}-topo.rtf', 'w') as f:
                topo.to_file(f)
            return -3

    # now run the tasks to minimize, stretch, orient along z, equilibrate, and sample
    tasks = C.tasks
    par = []
    for task in tasks:
        if task.taskname == 'md':
            needed = []
            na: NAMDScripter = C.tasks[0].scripters['namd']
            if extension == 'str':
                # this is a stream file (combined topology/parameter file); charmm_topfile is ABSOLUTE, so we need to add its RELATIVE path to the list of parameter files
                if not charmm_topfile in na.charmmff_config['standard'][extension]:
                    na.charmmff_config['standard'][extension].append(charmm_topfile)
            if charmm_topfile.endswith('sphingo.str'):
                needed = ['toppar_all36_carb_imlab.str',
                        'toppar_all36_lipid_lps.str']
            if charmm_topfile.endswith('detergent.str'):
                needed = ['toppar_all36_lipid_sphingo.str',
                        'toppar_all36_lipid_cholesterol.str',
                        'toppar_all36_carb_glycolipid.str']
            if charmm_topfile.endswith('inositol.str'):
                needed = ['toppar_all36_carb_glycolipid.str']
            if charmm_topfile.endswith('cardiolipin.str'):
                needed = ['toppar_all36_lipid_bacterial.str']
            if charmm_topfile.endswith('lps.str'):
                needed = ['toppar_all36_carb_imlab.str',
                        'toppar_all36_lipid_bacterial.str']

            for n in needed:
                if n not in na.charmmff_config['standard']['str']:
                    na.charmmff_config['standard']['str'].append(n)

            par = na.charmmff_config['standard']['prm'] + na.charmmff_config['standard']['str']

    result = C.do_tasks()
    for k, v in result.items():
        logger.debug(f'{k}: {v}')
        if v['result'] != 0:
            with open(f'{resid}-topo.rtf', 'w') as f:
                topo.to_file(f)
            return -1
        
    # now sample: write the {resid}-NN.pdb conformer set, either from the vacuum-MD trajectory
    # (sampler='md'), from cylinder-confined athermal MC on the minimized structure
    # (sampler='mc'), or -- for sterols (sampler='single') -- the single minimized+oriented
    # structure itself; the shared info.yaml block below measures whichever set was produced.
    if sampler == 'single':
        state = C.tasks[-1].get_current_artifact('state')
        src = state.pdb.name if state and state.pdb else None
        if not src or not os.path.exists(src):
            logger.warning(f'single-conformer {resid}: no minimized structure produced')
            return -1
        shutil.copyfile(src, f'{resid}-{0:0{digits}d}.pdb')
    elif sampler == 'mc':
        state = C.tasks[-1].get_current_artifact('state')
        mc_psf = state.psf.name if state and state.psf else f'{resid}-init.psf'
        mc_pdb = state.pdb.name if state and state.pdb else None
        if not mc_pdb or not os.path.exists(mc_pdb):
            logger.warning(f'MC {resid}: no minimized structure to sample from')
            return -1
        # size the confinement cylinder to the lipid's real footprint: auto from the chain count
        # (n_chains x per-chain area) unless an explicit cylinder_apl was requested
        if cylinder_apl is None:
            cylinder_apl = auto_cylinder_apl(chain_info, substream=substream)
            logger.info(f'MC {resid}: {len(chain_info)} acyl chain(s) -> auto cylinder_apl '
                        f'{cylinder_apl:.0f} A^2 (inflation {cylinder_inflation:.2f})')
        try:
            mc_result = _sample_and_write_mc_conformers(
                resid, mc_psf, mc_pdb, heads, tails, par, DB,
                nsamples=nsamples, digits=digits, cylinder_apl=cylinder_apl,
                cylinder_inflation=cylinder_inflation, n_equil=mc_n_equil,
                n_decorr=mc_n_decorr, seed=mc_seed,
                max_angle=mc_max_angle, radius_scale=mc_radius_scale,
                target_order=phase_order_target(phase))
            # the tuner resolves a trans-bias strength (0 for Ld) to hit the phase's order target
            torsion_bias = mc_result['torsion_bias']
            chain_order = mc_result['chain_order']
        except Exception as exc:
            logger.warning(f'MC conformer generation for {resid} failed: {exc}')
            return -1
    else:
        dcd = None
        # prefer 00-04-00_md-NVT.dcd, but if it doesn't exist, use 00-03-00_md-NVT.dcd, but
        # only for sterols
        possible_dcds = glob.glob('*.dcd')
        possible_dcds.sort(key=os.path.getmtime, reverse=True)
        for d in possible_dcds:
            if os.path.exists(d):
                dcd = d
                break
        if dcd is None:
            return -1
        if dcd == possible_dcds[-1] and substream != 'cholesterol':
            # we don't do any steered md for sterols, so the final
            # dcd file is 00-03-00_md-NVT.dcd; otherwise, if there is no
            # 00-04-00_md-NVT.dcd, we bail
            return -1
        W = C.tasks[0].scripters['vmd']
        W.newscript('sample')
        W.addline(f'mol new {resid}-init.psf')
        W.addline(f'mol addfile {dcd} waitfor all')
        W.addline(f'set a [atomselect top all]')
        W.addline(f'set ref [atomselect top all]')
        W.addline(f'$ref frame 0')
        W.addline(f'$ref move [vecscale -1 [measure center $ref]]')
        W.addline(r'for { set f 0 } { $f < [molinfo top get numframes] } { incr f } {')
        W.addline( '    $a frame $f')
        W.addline( '    $a move [measure fit $a $ref]')
        W.addline(f'    $a writepdb {resid}-[format %0{digits}d $f].pdb')
        W.addline(r'}')
        W.writescript()
        W.runscript()

    # and now we write the ancillary info file in YAML format
    psf = PSFContents(f'{resid}-init.psf')
    q = psf.get_charge()
    if np.abs(q) < 1.e-3:
        qstr = 0.0
    else:
        qstr = float(f'{q:.3f}')
    psfatoms = psf.atoms
    info = {
        'defined-in': charmm_topfile,
        'parameters': list(set(par)),
        'reference-atoms': {
            'heads': [dict(serial=p.serial, name=p.atomname) for p in psfatoms.data if p.atomname in heads],
            'tails': [dict(serial=p.serial, name=p.atomname) for p in psfatoms.data if p.atomname in tails],
        },
        'charge': qstr,
        'sampler': sampler,   # effective sampler ('md'|'mc'|'single'); sterols are forced to 'single'
        'conformers': []
        }

    q = '?' * digits
    pdbs = [os.path.basename(x) for x in glob.glob(f'{resid}-{q}.pdb')]
    logger.debug(f'found {len(pdbs)} pdbs: {pdbs}')
    for pdb in pdbs:
        entry = {}
        p = PDBParser(filepath=pdb).parse()
        pdbatoms = p.parsed['ATOM']
        entry['pdb'] = pdb
        head_z = np.array([x.z for x in pdbatoms if x.name in heads])
        tail_z = np.array([x.z for x in pdbatoms if x.name in tails])
        logger.debug(f'head_z: {head_z}, tail_z: {tail_z}')
        length = np.array([np.abs(x-y) for x,y in product(head_z, tail_z)]).min()
        entry['head-tail-length'] = float(f'{length:.3f}')
        max_internal_length = 0.0
        for i,j in product(pdbatoms,pdbatoms):
            internal_length = np.sqrt((i.z-j.z)**2+(i.y-j.y)**2+(i.x-j.x)**2)
            if internal_length > max_internal_length:
                max_internal_length = internal_length
        entry['max-internal-length'] = float(f'{max_internal_length:.3f}')
        info['conformers'].append(entry)

    # provenance: how this conformer set was generated, with *resolved* values (e.g. the auto-sized
    # cylinder_apl), so a later run can reproduce or audit it.  Written here so both the CLI
    # (make-pdb-collection) and the on-the-fly autocache paths record it identically.
    if sampler == 'mc':
        info['generation'] = {'sampler': 'mc', 'phase': phase, 'n_chains': len(chain_info),
                              'cylinder_apl': float(cylinder_apl),
                              'cylinder_inflation': float(cylinder_inflation),
                              'torsion_bias': float(torsion_bias), 'chain_order': chain_order,
                              'mc_n_equil': int(mc_n_equil), 'mc_n_decorr': int(mc_n_decorr),
                              'mc_seed': int(mc_seed), 'mc_max_angle': float(mc_max_angle),
                              'mc_radius_scale': float(mc_radius_scale)}
    elif sampler == 'single':
        info['generation'] = {'sampler': 'single'}
    else:
        info['generation'] = {'sampler': 'md', 'sample_steps': int(sample_steps),
                              'sample_temperature': float(sample_temperature)}

    info['synonym'] = synonym
    with open('info.yaml','w') as f:
        f.write(yaml.dump(info))
    return 0

def do_cleanup(resname, dirname):
    """ 
    Remove all files in the directory except for the ``init.tcl``, ``info.yaml``, and ``*.psf`` files.
    """
    cwd = os.getcwd()
    os.chdir(dirname)
    files = os.listdir('.')
    files.remove('init.tcl')
    files.remove('info.yaml')
    files.remove(f'{resname}-init.psf')
    for f in glob.glob(f'{resname}-*.pdb'):
        files.remove(f)
    for f in files:
        os.remove(f)
    os.chdir(cwd)

def do_resi(resi: str, DB: CHARMMFFContent, RM: ResourceManager = None,
            outdir: str = 'data', faildir: str = 'fails', force: bool = False,
            lenfac: float = 1.2, cleanup: bool = True, minimize_steps: int = 500,
            sample_steps: int = 5000, nsamples: int = 10,
            sample_temperature: float = 300, refic_idx: int = 0,
            force_constant: float = 1.0, borrow_ic_from: str = None,
            sampler: str = 'md', cylinder_apl: float = None, cylinder_inflation: float = 1.9,
            mc_n_equil: int = 20000, mc_n_decorr: int = 3000, mc_seed: int = 0,
            mc_max_angle: float = np.pi / 6, mc_radius_scale: float = 0.8,
            phase: str = None):
    """
    Manager function for :func:`do_psfgen`.  Makes sure it operates in the correct subdirectories and handles success/failure cases.

    Parameters
    ----------
    resi : str
        The residue ID for which to generate the PDB file.
    DB : :class:`~pestifer.charmmff.charmmffcontent.CHARMMFFContent`
        The database containing CHARMM residue information.
    lenfac : float, optional
        The lengthening factor to apply to the residue (default is 1.2).
    minimize_steps : int, optional
        The number of minimization steps to perform (default is 500).
    sample_steps : int, optional
        The total number of sampling steps to perform (default is 5000).
    nsamples : int, optional
        The number of samples to generate (default is 10).
    sample_temperature : float, optional
        The temperature to use for sampling (default is 300 K).
    refic_idx : int, optional
        The index of the reference IC to use (default is 0).
    force_constant : float, optional
        The force constant to use for sampling (default is 1.0).
    borrow_ic_from : str, optional
        The residue ID from which to borrow internal coordinates (default is None).

    """
    cwd = os.getcwd()
    successdir = os.path.join(outdir, resi)
    failuredir = os.path.join(faildir, resi)
    if (not os.path.exists(successdir)) or force:
        if os.path.exists(successdir):
            shutil.rmtree(successdir)
        if os.path.exists('tmp'): shutil.rmtree('tmp')
        os.mkdir('tmp')
        os.chdir('tmp')
        result = do_psfgen(resi, DB, RM=RM, lenfac=lenfac, minimize_steps=minimize_steps,
                            sample_steps=sample_steps, nsamples=nsamples,
                            sample_temperature=sample_temperature,
                            refic_idx=refic_idx, force_constant=force_constant,
                            borrow_ic_from=borrow_ic_from, sampler=sampler,
                            cylinder_apl=cylinder_apl, cylinder_inflation=cylinder_inflation,
                            mc_n_equil=mc_n_equil, mc_n_decorr=mc_n_decorr, mc_seed=mc_seed,
                            mc_max_angle=mc_max_angle, mc_radius_scale=mc_radius_scale,
                            phase=phase)
        os.chdir(cwd)
        if result == 0:
            if cleanup: 
                do_cleanup(resi, 'tmp')
            shutil.move('tmp', os.path.join(outdir, resi))
        elif result == -2:
            my_logger(f'RESI {resi} is not found', logger.warning, just='^', frame='*', fill='*')
        else:
            if os.path.exists(failuredir):
                shutil.rmtree(failuredir)
            shutil.move('tmp', failuredir)
    else:
        logger.info(f'RESI {resi} built previously; use \'--force\' to recalculate')

def make_pdb_collection(args):
    """
    Make a PDBCollection from the CHARMMFFContent.
    This function will create a PDBCollection from the CHARMMFFContent, either
    for a specific RESI or for all RESIs in a specified stream.
    If a specific RESI is provided, it will create a collection member for that RESI
    and save it in the specified output directory.
    If a stream ID is provided, it will create a collection for all RESIs in that stream.
    If no stream ID is provided, it will create a collection for all RESIs in the CHARMMFFContent.
    The output directory will be created if it does not exist, and if the ``--fail-dir`` argument is provided,
    it will create a directory for failed RESIs in that directory.
    If the ``--force`` argument is set, it will force the recalculation of the RESI, even if it has been built previously.
    If the ``--cleanup`` argument is set, it will remove all files in the RESI directory except for the init.tcl, info.yaml, and psf files.
    The ``--lenfac``, ``--minimize-steps``, ``--sample-steps``, ``--nsamples``, ``--sample-temperature``, ``--refic-idx``, and ``--force-constant`` arguments
    will be passed to the :func:`do_resi` function to control the sampling and equilibration of the RESI.
    """
    streamID = args.streamID # if provided, we will make a collection from RESIs in this stream
    substreamID = args.substreamID
    resname = args.resname # if provided, we will only make a collection member for this RESI
    topfile = args.topfile
    charmmff_config = {'release': args.charmmff_release} if getattr(args, 'charmmff_release', '') else {}
    RM = ResourceManager(charmmff_config=charmmff_config)
    CC = RM.charmmff_content
    CC.provision()
    
    if topfile is not None:
        if not os.path.exists(topfile):
            raise FileNotFoundError(f'Topology file {topfile} does not exist')
        CC.add_topology(topfile)

    if not resname == '' and not resname in CC:
        logger.warning(f'RESI {resname} not found in CHARMMFFContent; will not build PDB collection for it')
        exit(-1)

    outdir = args.output_dir
    if outdir is None or outdir == '':
        if streamID is not None:
            outdir = f'{streamID}'
        else:
            outdir = 'data'
    faildir = args.fail_dir
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not os.path.exists(faildir):
        os.mkdir(faildir)
    if os.path.exists('tmp'):
        shutil.rmtree('tmp')
    
    # sampler selection (default 'md' legacy); for 'mc', cylinder_apl=None means auto-size from the
    # per-lipid acyl-chain count (an explicit --cylinder-apl overrides).  Sterols are always forced to
    # 'single' inside do_psfgen regardless of this choice.
    sampler_kwargs = dict(
        sampler=getattr(args, 'sampler', 'md'),
        cylinder_apl=getattr(args, 'cylinder_apl', None),
        cylinder_inflation=getattr(args, 'cylinder_inflation', 1.9),
    )

    if resname is not None and resname != '':
        my_logger(f'RESI {resname}', logger.info, just='^', frame='*', fill='*')
        do_resi(resname, CC, RM=RM, outdir=outdir, faildir=faildir, force=args.force,
                 cleanup=args.cleanup, lenfac=args.lenfac,
                 minimize_steps=args.minimize_steps, sample_steps=args.sample_steps,
                 nsamples=args.nsamples, sample_temperature=args.sample_temperature,
                 refic_idx=args.refic_idx, force_constant=args.force_constant,
                 borrow_ic_from=args.take_ic_from, **sampler_kwargs)
    else:
        active_resnames = CC.get_resnames_of_streamID(streamID, substreamID=substreamID)
        logger.debug(f'stream {streamID} substream {substreamID} active_resnames: {active_resnames}')
        nresi = len(active_resnames)
        for i, r in enumerate(active_resnames):
            my_logger(f'RESI {r} ({i+1}/{nresi})', logger.info, just='^', frame='*', fill='*')
            do_resi(r, CC, RM=RM, outdir=outdir, faildir=faildir, force=args.force,
                    cleanup=args.cleanup, lenfac=args.lenfac,
                    minimize_steps=args.minimize_steps, sample_steps=args.sample_steps,
                    nsamples=args.nsamples, sample_temperature=args.sample_temperature,
                    refic_idx=args.refic_idx, force_constant=args.force_constant,
                    **sampler_kwargs)

    # if the faildir is empty, remove it
    if len(os.listdir(faildir)) == 0:
        os.rmdir(faildir)
    else:
        logger.warning(f'Failures in {faildir}; see the files there for details')

