# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
This module defines functions that facilitate building PDB files of RESI's using ICs, and equilibrating them and generating samples that can be used as input to packmol.
"""
import glob
import logging
import os
import shutil
import yaml
import numpy as np

from itertools import product

from pidibble.pdbparse import PDBParser
from ..core.resourcemanager import ResourceManager
from .charmmffcontent import CHARMMFFResiDatabase
from ..core.config import Config
from ..core.controller import Controller
from ..psfutil.psfcontents import PSFContents
from ..core.scripters import PsfgenScripter, VMDScripter
from ..core.stringthings import my_logger

logger=logging.getLogger(__name__)

def do_psfgen(resid,DB,lenfac=1.2,minimize_steps=500,sample_steps=5000,nsamples=10,sample_temperature=300,refic_idx=0,force_constant=1.0,borrow_ic_from=None):
    """ 
    Generate a PDB file for a RESI using psfgen, and sample it.
    """
    if nsamples>sample_steps:
        raise ValueError(f'nsamples ({nsamples}) must be less than or equal to sample_steps ({sample_steps})')
    digits=len(str(nsamples))
    topo=DB.get_resi(resid)
    synonym=topo.synonym
    meta=topo.metadata
    charmm_topfile=meta.get('charmmtopfile',None)
    logger.debug(f'do_psfgen for {resid} with metadata {meta}')
    if borrow_ic_from:
        logger.debug(f'borrowing ICs from {borrow_ic_from}')
        take_my_ics=DB.get_resi(borrow_ic_from)
        topo.copy_ICs_from(take_my_ics)
    charmm_topfile,stream,substream=meta['charmmtopfile'],meta['streamID'],meta['substreamID']
    charmm_topfile=os.path.basename(charmm_topfile)
    logger.debug(f'charmm_topfile: {charmm_topfile}, stream: {stream}, substream: {substream}')
    substream_overrides=DB.overrides.get('substreams',{})
    ss_override=substream_overrides.get(resid,None)
    if ss_override is not None:
        substream=ss_override
        logger.debug(f'Overriding substream for {resid} from {meta["substreamID"]} to {substream}')
        meta['substreamID']=substream
    if charmm_topfile is None:
        return -2
    if topo.error_code!=0:
        topo.show_error(logger.warning)
        # my_logger(f'Parsing error for {resid}',logger.warning,just='^',frame='!',fill='!') 
        with open(f'{resid}-topo.rtf','w') as f:
            topo.to_file(f)
        # return -1

    # for b in topo.bonds:
    #     logger.debug(f'{b.name1}-{b.name2}')
    my_logger(f'{os.path.basename(charmm_topfile)}',logger.info,just='^',frame='*',fill='*')
    
    tasklist=[
            {'restart':{'psf': f'{resid}-init.psf','pdb': f'{resid}-init.pdb'}},
            {'md': {'ensemble':'minimize','nsteps':0,'minimize':minimize_steps,'dcdfreq':0,'xstfreq':0,'temperature':100}},
            ]
    heads=[]
    if stream=='lipid':
        topo.lipid_annotate()
        ano=topo.annotation
        heads,tails,shortest_paths=ano.get('heads',[]),ano.get('tails',[]),ano.get('shortest_paths',{})
        if substream=='cholesterol' or substream=='model':
            pass
        elif substream=='detergent':
            logger.debug('detergent')
            logger.debug(f'heads: {heads}. tails: {tails}, shortest_paths: {shortest_paths}')
            pass
        else:
            if shortest_paths and any([len(v)>0 for v in shortest_paths.values()]) and heads and tails: # and len(tails)==2:
                # run a non-equilibrium MD simulation to bring the tails together to make a 'standard' conformation
                base_md={'ensemble':'NVT','nsteps':15000,'dcdfreq':100,'xstfreq':100,'temperature':100}
                groups={'repeller':{'atomnames': [heads[0]]}}
                distances={}
                harmonics={}
                for i in range(len(tails)):
                    groups[f'tail{i+1}']={'atomnames': [tails[i]]}
                    # if len(tails)<4:
                    distances[f'repeller_tail{i+1}']={'groups': ['repeller',f'tail{i+1}']}
                dists={}
                for i in range(len(tails)):
                    dists[i]=shortest_paths[heads[0]][tails[i]]*lenfac
                for i in range(len(tails)):
                    for j in range(i+1,len(tails)):
                        distances[f'tail{i+1}_tail{j+1}']={'groups': [f'tail{i+1}',f'tail{j+1}']}
                        harmonics[f'tail{i+1}_tail{j+1}_attract']={
                                'colvars': [f'tail{i+1}_tail{j+1}'],
                                'forceConstant': force_constant,
                                'distance':[4.0]}
                name='repeller_tail'+''.join(f'{i+1}' for i in range(len(tails)))
                colvars=[]               
                distance=[]
                for i in range(len(tails)):
                    colvars.append(f'repeller_tail{i+1}')
                    distance.append(dists[i])
                harmonics[name]={'colvars':colvars,'distance':distance,'forceConstant':force_constant}
                colvar_specs={'groups':groups,'distances':distances,'harmonics':harmonics}
                base_md['colvar_specs']=colvar_specs.copy()
                assert 'minimize' not in base_md
                logger.debug(f'base_md {base_md}')
                tasklist.append({'md':base_md})
    else:
        my_logger(f'{resid} is from stream {stream}',logger.debug,just='^',frame='!',fill='!')
        with open(f'{resid}-topo.rtf','w') as f:
            topo.to_file(f)
    logger.debug(f'heads: {heads}')
    if len(heads)>0:
        # reorient molecule so head is at highest z position if molecule is rotated about its COM
        tasklist.append({'manipulate':{'mods':{'orient':[f'z,{heads[0]}']}}})
    # do a conformer-generation MD simulation with the external forces dialed down a bit
    base_md={'ensemble':'NVT','nsteps':sample_steps,'dcdfreq':sample_steps//nsamples,'xstfreq':100,'temperature':sample_temperature}
    if substream not in ['cholesterol','detergent']:
        for cv,spec in colvar_specs['harmonics'].items():
            if 'forceConstant' in spec:
                spec['forceConstant']*=0.1
        base_md['colvar_specs']=colvar_specs.copy()
    tasklist.append({'md':base_md})

    config=Config(quiet=True)
    C=Controller(config=config,userspecs={'title':f'Build PDBCollection entry for {resid}','tasks':tasklist})
    # First we de-novo generate a pdb/psf file using seeded internal coordinates
    W=PsfgenScripter(C.config)
    logger.debug(f'charmm_topfile from resiDB entry {topo.resname}: {charmm_topfile}')
    if not charmm_topfile in W.charmmff_config['standard']['topologies']:
        W.charmmff_config['standard']['topologies'].append(charmm_topfile)
    if charmm_topfile.endswith('detergent.str'):
        needed='toppar_all36_lipid_sphingo.str'
        W.charmmff_config['standard']['topologies'].append(needed)
    if charmm_topfile.endswith('initosol.str'):
        needed='toppar_all36_carb_glycolipid.str'
        W.charmmff_config['standard']['topologies'].append(needed)

    W.newscript('init')
    if topo.to_psfgen(W,refic_idx=refic_idx)==-1:
        my_logger(f'No valid IC\'s for {resid}',logger.warning,just='^',frame='!',fill='!') 
        with open(f'{resid}-topo.rtf','w') as f:
            topo.to_file(f)
        return -1
    W.writescript(f'{resid}-init')
    W.runscript()
    with open('init.log','r') as f:
        loglines=f.read().split('\n')
    for l in loglines:
        if 'Warning: failed to guess coordinates' in l:
            my_logger(f'psfgen failed for {resid}',logger.warning,just='^',frame='!',fill='!') 
            with open(f'{resid}-topo.rtf','w') as f:
                topo.to_file(f)
            return -3

    # now run the tasks to minimize, stretch, orient along z, equilibrate, and sample
    tasks=C.tasks
    par=[]
    for task in tasks:
        if task.taskname=='md':
            needed=[]
            if 'toppar' in charmm_topfile: 
                # this is a combined topology/parameter file; charmm_topfile is ABSOLUTE, so we need to add its RELATIVE path to the list of parameter files
                if not charmm_topfile in task.scripters['namd'].charmmff_config['standard']['parameters']:
                    task.scripters['namd'].charmmff_config['standard']['parameters'].append(charmm_topfile)
            if charmm_topfile.endswith('sphingo.str'):
                needed=['toppar_all36_carb_imlab.str',
                        'toppar_all36_lipid_lps.str']
            if charmm_topfile.endswith('detergent.str'):
                needed=['toppar_all36_lipid_sphingo.str',
                        'toppar_all36_lipid_cholesterol.str',
                        'toppar_all36_carb_glycolipid.str']
            if charmm_topfile.endswith('inositol.str'):
                needed=['toppar_all36_carb_glycolipid.str']
            if charmm_topfile.endswith('cardiolipin.str'):
                needed=['toppar_all36_lipid_bacterial.str']
            if charmm_topfile.endswith('lps.str'):
                needed=['toppar_all36_carb_imlab.str',
                        'toppar_all36_lipid_bacterial.str']

            for n in needed:
                if n not in task.scripters['namd'].charmmff_config['standard']['parameters']:
                    task.scripters['namd'].charmmff_config['standard']['parameters'].append(n)

            par=task.scripters['namd'].charmmff_config['standard']['parameters']

    result=C.do_tasks()
    for k,v in result.items():
        logger.debug(f'{k}: {v}')
        if v['result']!=0:
            with open(f'{resid}-topo.rtf','w') as f:
                topo.to_file(f)
            return -1
        
    # now sample
    dcd=None
    # prefer 00-04-00_md-NVT.dcd, but if it doesn't exist, use 00-03-00_md-NVT.dcd, but
    # only for sterols
    possible_dcds=['00-04-00_md-NVT.dcd','00-03-00_md-NVT.dcd']
    for d in possible_dcds:
        if os.path.exists(d):
            dcd=d
            break
    if dcd is None:
        return -1
    if dcd==possible_dcds[-1] and substream!='cholesterol':
        # we don't do any steered md for sterols, so the final
        # dcd file is 00-03-00_md-NVT.dcd; otherwise, if there is no 
        # 00-04-00_md-NVT.dcd, we bail
        return -1
    W=VMDScripter(C.config)
    W.newscript('sample')
    W.addline(f'mol new {resid}-init.psf')
    W.addline(f'mol addfile {dcd} waitfor all')
    W.addline(f'set a [atomselect top all]')
    # W.addline(f'set b [atomselect top noh]')
    W.addline(f'set ref [atomselect top all]')
    W.addline(f'$ref frame 0')
    # W.addline(f'set bref [atomselect top noh]')
    W.addline(f'$ref move [vecscale -1 [measure center $ref]]')
    # W.addline(f'$bref move [vecscale -1 [measure center $bref]]')
    W.addline(r'for { set f 0 } { $f < [molinfo top get numframes] } { incr f } {')
    W.addline( '    $a frame $f')
    W.addline( '    $a move [measure fit $a $ref]')
    W.addline(f'    $a writepdb {resid}-[format %0{digits}d $f].pdb')
    # W.addline( '    $b frame $f')
    # W.addline( '    $b move [measure fit $b $bref]')
    # W.addline(f'    $b writepdb {resid}-noh-[format %0{digits}d $f].pdb')
    W.addline(r'}')
    W.writescript()
    W.runscript()

    # and now we write the ancillary info file in YAML format
    psf=PSFContents(f'{resid}-init.psf')
    q=psf.get_charge()
    if np.abs(q)<1.e-3:
        qstr=0.0
    else:
        qstr=float(f'{q:.3f}')
    psfatoms=psf.atoms
    info={
        'defined-in':charmm_topfile,
        'parameters':list(set(par)),
        'reference-atoms':{
            'heads':[dict(serial=p.serial,name=p.name) for p in psfatoms if p.name in heads],
            'tails':[dict(serial=p.serial,name=p.name) for p in psfatoms if p.name in tails],
            },
        'charge':qstr,
        'conformers':[]
        }

    q='?'*digits
    pdbs=[os.path.basename(x) for x in glob.glob(f'{resid}-{q}.pdb')]
    logger.debug(f'found {len(pdbs)} pdbs')
    for pdb in pdbs:
        entry={}
        p=PDBParser(PDBcode=os.path.splitext(pdb)[0]).parse()
        pdbatoms=p.parsed['ATOM']
        entry['pdb']=pdb
        head_z=np.array([x.z for x in pdbatoms if x.name in heads])
        tail_z=np.array([x.z for x in pdbatoms if x.name in tails])
        length=np.array([np.abs(x-y) for x,y in product(head_z,tail_z)]).min()
        entry['head-tail-length']=float(f'{length:.3f}')
        max_internal_length=0.0
        for i,j in product(pdbatoms,pdbatoms):
            internal_length=np.sqrt((i.z-j.z)**2+(i.y-j.y)**2+(i.x-j.x)**2)
            if internal_length>max_internal_length:
                max_internal_length=internal_length
        entry['max-internal-length']=float(f'{max_internal_length:.3f}')
        info['conformers'].append(entry)

    info['synonym']=synonym
    with open('info.yaml','w') as f:
        f.write(yaml.dump(info))
    return 0

def do_cleanup(resname,dirname):
    """ 
    Remove all files in the directory except for the init.tcl, info.yaml, and psf files.
    """
    cwd=os.getcwd()
    os.chdir(dirname)
    files=glob.glob('*')
    files.remove('init.tcl')
    files.remove('info.yaml')
    files.remove(f'{resname}-init.psf')
    for f in glob.glob(f'{resname}-*.pdb'):
        files.remove(f)
    for f in files:
        os.remove(f)
    os.chdir(cwd)

def do_resi(resi,DB,outdir='data',faildir='fails',force=False,lenfac=1.2,cleanup=True,minimize_steps=500,sample_steps=5000,nsamples=10,sample_temperature=300,refic_idx=0,force_constant=1.0,borrow_ic_from=None):
    """
    Build a RESI using psfgen and sample it.
    This function checks if the RESI has been built previously, and if not, it will
    create a new directory for the RESI in the specified output directory.
    If the RESI has been built previously, it will skip the build step unless the `force` argument is set to True.
    If the RESI is not found in the CHARMMFFResiDatabase, it will log a warning and return -2.
    If the RESI is found, it will create a new directory for the RESI in the specified output directory,
    and if the `cleanup` argument is set to True, it will remove all files in the directory except for the init.tcl, info.yaml, and psf files.
    """
    cwd=os.getcwd()
    successdir=os.path.join(outdir,resi)
    failuredir=os.path.join(faildir,resi)
    if (not os.path.exists(successdir)) or force:
        if os.path.exists(successdir):
            shutil.rmtree(successdir)
        if os.path.exists('tmp'): shutil.rmtree('tmp')
        os.mkdir('tmp')
        os.chdir('tmp')
        result=do_psfgen(resi,DB,lenfac=lenfac,minimize_steps=minimize_steps,sample_steps=sample_steps,nsamples=nsamples,sample_temperature=sample_temperature,refic_idx=refic_idx,force_constant=force_constant,borrow_ic_from=borrow_ic_from)
        os.chdir(cwd)
        if result==0:
            if cleanup: do_cleanup(resi,'tmp')
            shutil.move('tmp',os.path.join(outdir,resi))
        elif result==-2:
            my_logger(f'RESI {resi} is not found',logger.warning,just='^',frame='*',fill='*') 
        else:
            if os.path.exists(failuredir):
                shutil.rmtree(failuredir)
            shutil.move('tmp',failuredir)
    else:
        logger.info(f'RESI {resi} built previously; use \'--force\' to recalculate')

def make_pdb_collection(args):
    """
    Make a PDBCollection from the CHARMMFFResiDatabase.
    This function will create a PDBCollection from the CHARMMFFResiDatabase, either
    for a specific RESI or for all RESIs in a specified stream.
    If a specific RESI is provided, it will create a collection member for that RESI
    and save it in the specified output directory.
    If a stream ID is provided, it will create a collection for all RESIs in that stream.
    If no stream ID is provided, it will create a collection for all RESIs in the CHARMMFFResiDatabase.
    The output directory will be created if it does not exist, and if the `--fail-dir` argument is provided,
    it will create a directory for failed RESIs in that directory.
    If the `--force` argument is set, it will force the recalculation of the RESI, even if it has been built previously.
    If the `--cleanup` argument is set, it will remove all files in the RESI directory except for the init.tcl, info.yaml, and psf files.
    The `--lenfac`, `--minimize-steps`, `--sample-steps`, `--nsamples`, `--sample-temperature`, `--refic-idx`, and `--force-constant` arguments
    will be passed to the `do_psfgen` function to control the sampling and equilibration of the RESI.
    """
    streamID=args.streamID # if provided, we will make a collection from RESIs in this stream
    substreamID=args.substreamID
    resname=args.resname # if provided, we will only make a collection member for this RESI
    topfile=args.topfile # if provided, we will use this topology file instead of the one in the CHARMMFFResiDatabase; stream name is extracted
    loglevel_numeric=getattr(logging,args.log_level.upper())
    if args.log_file:
        if os.path.exists(args.log_file):
            shutil.copyfile(args.log_file,args.log_file+'.bak')
        logging.basicConfig(filename=args.log_file,filemode='w',format='%(asctime)s %(name)s %(message)s',level=loglevel_numeric)
    console=logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter=logging.Formatter('%(levelname)s> %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    RM=ResourceManager()
    CC=RM.charmmff_content
    DB=CHARMMFFResiDatabase(CC)
    if streamID is not None:
        DB.add_stream(streamID)
    
    if topfile is not None:
        if not os.path.exists(topfile):
            raise FileNotFoundError(f'Topology file {topfile} does not exist')
        DB.add_topology(topfile)

    if not resname == '' and not resname in DB:
        logger.warning(f'RESI {resname} not found in CHARMMFFResiDatabase; will not build PDB collection for it')
        exit(-1)

    outdir=args.output_dir
    if outdir is None or outdir=='':
        if streamID is not None:
            outdir=f'{streamID}'
        else:
            outdir='data'
    faildir=args.fail_dir
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not os.path.exists(faildir):
        os.mkdir(faildir)
    if os.path.exists('tmp'):
        shutil.rmtree('tmp')
    
    if resname is not None and resname != '':
        my_logger(f'RESI {resname}',logger.info,just='^',frame='*',fill='*')
        do_resi(resname,DB,outdir=outdir,faildir=faildir,force=args.force,cleanup=args.cleanup,lenfac=args.lenfac,minimize_steps=args.minimize_steps,sample_steps=args.sample_steps,nsamples=args.nsamples,sample_temperature=args.sample_temperature,refic_idx=args.refic_idx,force_constant=args.force_constant,borrow_ic_from=args.take_ic_from)
    else:
        active_resnames=DB.get_resnames_of_streamID(streamID,substreamID=substreamID)
        logger.debug(f'active_resnames: {active_resnames}')
        nresi=len(active_resnames)
        for i,r in enumerate(active_resnames):
            my_logger(f'RESI {r} ({i+1}/{nresi})',logger.info,just='^',frame='*',fill='*')
            do_resi(r,DB,outdir=outdir,faildir=faildir,force=args.force,cleanup=args.cleanup,lenfac=args.lenfac,minimize_steps=args.minimize_steps,sample_steps=args.sample_steps,nsamples=args.nsamples,sample_temperature=args.sample_temperature,refic_idx=args.refic_idx,force_constant=args.force_constant)

    # if the faildir is empty, remove it
    if len(os.listdir(faildir))==0:
        os.rmdir(faildir)
    else:
        logger.warning(f'Failures in {faildir}; see the files there for details')

