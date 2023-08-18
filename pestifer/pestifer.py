
import sys
import argparse as ap
import os
import shutil
import random
import logging
logger=logging.getLogger(__name__)
from datetime import date, datetime

from .controller import Controller



def _main():
    parser=ap.ArgumentParser()
    parser.add_argument('config',help='input configuration file (yaml)')
    parser.add_argument('-log',help='log file name',default='pestifer.log')
    parser.add_argument('--loglevel',default='debug',help='logging level (info)')
    args=parser.parse_args()

    loglevel=args.loglevel
    loglevel_numeric=getattr(logging, loglevel.upper())
    if os.path.exists(args.log):
        shutil.copyfile(args.log,args.log+'.bak')
    logging.basicConfig(filename=args.log,filemode='w',format='%(asctime)s %(message)s',level=loglevel_numeric)
    logger.info(f'pestifer runtime begins')

    C=Controller(args.config)
    C.do_steps()
    
    logger.info('pestifer runtime ends.')

def cli():

    _main()
    
def old_main():
    seed=random.randint(0,100000)
    parser=ap.ArgumentParser()
    print('pestifer.py {} / python {}'.format(date.today(),sys.version.replace('\n',' ').split(' ')[0]))
    i=1
    Molecules=[]
    Mut=[]
    Clv=[]
    Uss=[]
    Mis=[]
    UIC=[]
    # defaults
    psfgen='mkpsf.tcl'
    CTopo=['top_all36_prot.rtf','top_all35_ethers.rtf','top_all36_cgenff.rtf','top_all36_lipid.rtf',
           'top_all36_na.rtf','stream/carb/toppar_all36_carb_glycopeptide.str']
    # default local topologies: these are specially modified charmm str files that get rid of things that PSFGEN can't handle
    LocTopo=['top_all36_carb.rtf','toppar_water_ions.str']
    StdParamFiles=['par_all36_prot.prm','par_all36_carb.prm','par_all36_lipid.prm','par_all36_na.prm','par_all36_cgenff.prm','stream/carb/toppar_all36_carb_glycopeptide.str']
    LocalParamFiles=['toppar_water_ions.str']
    PDBAliases=['residue HIS HSD','atom ILE CD1 CD','residue NAG BGNA','atom BGNA C7 C',
                        'atom BGNA O7 O','atom BGNA C8 CT','atom BGNA N2 N','residue SIA ANE5',
                        'atom ANE5 C10 C','atom ANE5 C11 CT','atom ANE5 N5 N','atom ANE5 O1A O11',
                        'atom ANE5 O1B O12','atom ANE5 O10 O','atom VCG C01 C1','atom VCG C01 C1','atom VCG C02 C2',
                        'atom VCG C03 C3','atom VCG C04 C4','atom VCG C05 C5','atom VCG C06 C6','atom VCG C07 C7',
                        'atom VCG C08 C8','atom VCG C09 C9','residue EIC LIN']
    PostMod={}
    PostMod['center_protein']=True
    prefix='x01_'
    fixConflicts=True
    PostMod['do_loop_mc']=False
    PostMod['do_gly_mc']=False
    PostMod['Crot']=[]

    parser.add_argument('-inpdb',default=[],nargs='+',metavar='<?.pdb>',type=str,help='Name(s) of pdb file to parse; First is treated as the base molecule')
    parser.add_argument('-incif',default=[],nargs='+',metavar='<?.cif>',type=str,help='Name(s) of mmCIF file to parse; First is treated as the base molecule')

    parser.add_argument('-ba','--biological-assembly',metavar='#',default=0,type=int,help='Biological assembly to construct; one may be selected from those defined in PDB/mmCIF metadata.  If not specified, the explicit model is built.')
    parser.add_argument('-charmmtopo',metavar='<name> ...',nargs='+',default=[],help='Additional (standard) CHARMM topology files in your CHARMM directory')
    parser.add_argument('-loctopo',metavar='<name> ...',nargs='+',default=[],help='Additional (local) CHARMM topology files in the $PSFGEN_BASEDIR/charmm directory')
    parser.add_argument('-charmmparam',metavar='<name> ...',nargs='+',default=[],help='Additional (standard) CHARMM parameter files in your CHARMM directory')
    parser.add_argument('-locparam',metavar='<name> ...',nargs='+',default=[],help='Additional (local) CHARMM parameter files in the $PSFGEN_BASEDIR/charmm directory')
    parser.add_argument('-prefix',metavar='<str>',default='x01_',help='Output PDB/PSF prefix; each file name will have the format <prefix><pdbcode>.pdb/psf, where <pdbcode> is the 4-letter PDB code of the base molecule.')
    parser.add_argument('-psfgen',metavar='<name>',default='mkpsf.tcl',help='name of TcL script generated as input to VMD/psfgen')
    parser.add_argument('-ignore',metavar='X ...',nargs='+',default=[],type=str,help='Specify chain(s) to ignore from asymmetric unit.')
    parser.add_argument('-includeonly',metavar='X ...',nargs='+',default=[],type=str,help='Specify only chain(s) to include from asymmetric unit')
    parser.add_argument('-modsfile',metavar='<name> ...',nargs='+',default=[],type=ModsFile,help='One (or more) modifications file(s) to rule them all.')
    parser.add_argument('-mut',metavar='C_OrrrN [C_OrrrN] ...',nargs='+',default=[],type=Mutation,help='One or more point-mutation specifications.  Format: C is chainID, O is one-letter residue code to mutate FROM, rrr is sequence number (can be any number of digits), and N is one-letter residue code to mutate TO.  Multiple mutation instances can be specified with one -mut.  Mutations are automatically replicated if there are BIOMT transformations.')
    parser.add_argument('-mutfile',metavar='<name>',default='',help='Input file listing mutation specifications')
    parser.add_argument('-delete',metavar='C_Orrr [C_Orrr] ...',nargs='+',default=[],type=Deletion,help='One or more single-residue deletion specifications.  Format: C is chainID, O is one-letter residue code, rrr is sequence number (can be any number of digits).  Multiple deletion instances can be specified with one -del.  Deletions are automatically replicated if there are BIOMT transformations.')
    parser.add_argument('-deletefile',metavar='<name>',default='',help='Input file listing deletion specifications')
    parser.add_argument('-clv',metavar='PrrrC [PrrrC] ...',nargs='+',default=[],type=Cleavage,help='One or more cleavage-site specifications.  Format: P is parent chain ID, rrr is residue number immediately N-terminal to the cleavage site, and C is the daughter chain ID that will begin immediately C-terminal to cleavage site. Multiple cleavage instances can be specified after one -clv.')
    parser.add_argument('-clvfile',metavar='<name>',default='',help='input file listing all cleavages (as an alternative to issuing multiple -clv arguments)')
    parser.add_argument('-gra',metavar='<str>,A:XXX-YYY,ZZZ,C:BBB ...',nargs='+',default=[],type=Graft,help='One or more graft specifications; graft resids XXX-YYY of chain A in pdb <str> to chain C of base molecule by overlapping resid ZZZ of chain A of graft and resid BBB of chain C of base.  Grafts are automatically replicated if there are BIOMT transformations.')
    parser.add_argument('-grafile',metavar='<name>',default='',help='Input file listing all grafts (as an alternative to issuing multiple -gra arguments)')
    parser.add_argument('-att',metavar='<str>,A:XXX-YYY,ZZZ,B:QQQ,C:BBB ...',nargs='+',default=[],type=Attach,help='One or more attachment specifications. Format: attach resids XXX-YYY of chain A using resid ZZZ (between XXX and YYY) in pdb <str> to chain C of base molecule at resid BBB by aligning resid QQQ of chain B from source to resid BBB of chain C of base.  CURRENTLY UNIMPLEMENTED!!!')
    parser.add_argument('-attfile',metavar='<name>',default='',help='Input file listing all attachments (as an alternative to issuing multiple -att arguments)')
    parser.add_argument('-crot',metavar='<str>,A,XXX[,YYY],### ...',default=[],nargs='+',type=Crot,help='One or more torsion rotation specifications. Format: <str> is one of phi, psi, omega, chi1, or chi2.  A is the chainID, XXX is the resid of owner of torson, and YYY (if given) marks the end of the sequence C-terminal to XXX that is reoriented by a backbone rotation. ### is the degrees of rotation.  C-rotations are automatically replicated if there are BIOMT transformations.')
    parser.add_argument('-crotfile',metavar='<name>',default='',help='Input file listing all torsion rotations requested (as an alternative to issuing multiple -crot arguments)')
    parser.add_argument('-ssbond',metavar='X_###-Y_### ...',default=[],nargs='+',type=SSBond,help='One or more disulfide bond specification(s) not already in input PDB/CIF: Format: X,Y are chainIDs and ### are resids; if residues are not CYS in wt or by mutations, there is no effect.  Because SSBonds can join chains together, they are NOT automatially replicated if there are BIOMT transformations.')
    parser.add_argument('-xssbond','--exclude-ssbond',metavar='X_###-Y_### ...',default=[],nargs='+',type=SSBond,help='One or more disulfide bond specification(s) in input PDB/CIF that you want to exclude: Format: X,Y are chainIDs and ### are resids; if residues are not CYS in wt or by mutations, there is no effect.  Because SSBonds can join chains together, they are NOT automatially replicated if there are BIOMT transformations.')
    parser.add_argument('-ssfile',metavar='<name>',default='',help='input file listing all disulfide bonds to add that are not already in the PDB file (as an alternative to issuing multiple -ssbond arguments)')
    parser.add_argument('-xssfile','--exclude-ssfile',metavar='<name>',default='',help='input file listing all disulfide bonds already in the PDB file that you want to exclude (as an alternative to issuing multiple -xssbond arguments)')
    parser.add_argument('-link',metavar='string',default=[],action='append',type=Link,help='PDB-format LINK record; must have exact spacing; multiple "-link" options can be supplied.')
    parser.add_argument('-linkfile',metavar='<name>',default='',help='Input file with PDB-format LINK records the user would like to enforce that are not in the PDB/CIF file')
    parser.add_argument('-missing',metavar='string',default=[],action='append',type=Missing,help='PDB-format REMARK 465 record; must have exact spacing; multiple "-missing" options can be supplied.')
    parser.add_argument('-missingfile',metavar='<name>',default='',help='Input file with PDB-format REMARK 465 records the user would like to enforce that are not in the PDB/CIF file (may be useful for insertions)')
    parser.add_argument('-pdbalias',metavar='<str>',default=[],nargs='+',help='One or more psfgen-formatted pdbalias with commas for spaces')
    parser.add_argument('-pdbaliasfile',metavar='<str>',default='',help='Input file containing psfgen-formatted pdbaliases')
    parser.add_argument('-logdcd',metavar='<name>.dcd',default='',help='Name of dcd logging file')
    parser.add_argument('-logevery',metavar='<int>',default=1,help='Number of MC accepts between successive frame logging')
    parser.add_argument('-logsaveevery',metavar='<int>',default=1,help='Number of MC accepts between log writes to disk')
 #   parser.add_argument('-rlxloops',action='store_true',help='asks psfgen to use the loopMC module to relax modeled-in loops of residues missing from PDB')
#    parser.add_argument('-rlxmc',action='store_true',help='Asks psfgen to use do_multiflex_mc module to relax modeled-in loops of residues missing from PDB and glycans')
#    parser.add_argument('-loopmcparams',metavar='<param1=val1,param2=val2,...>',default='',help='Loop Monte Carlo parameters')
#    parser.add_argument('-rlxmcparams',metavar='<param1=val1,param2=val2,...>',default='',help='Loop Monte Carlo parameters')
#    parser.add_argument('-rlxgly',action='store_true',help='asks psfgen to use the loopMC module to relax modeled-in glycans missing from PDB')
#    parser.add_argument('-glymcparams',metavar='<param1=val1,param2=val2,...>',default='',help='Glycan Monte Carlo parameters')
    parser.add_argument('-smdclose',action='store_true',help='Asks psfgen to prep for steered MD simulation to close loops')
    parser.add_argument('-smdcloseparams',metavar='<param1=val1,param2=val2,...>',default='',help='Parameters for steered MD for closing loops')
    parser.add_argument('-includeTerminalLoops',action='store_true',default=False)
    parser.add_argument('-namdparams',metavar='<param1=val1,param2=val2,...>',default='',help='Parameters for NAMD runs')
    parser.add_argument('-fixconflicts',action='store_true',help='If set, residues in CONFLICT are mutated to their proper identities')
    parser.add_argument('-fixengineeredmutations',action='store_true',help='If set, residues in ENGINEERED MUTATIONS are mutated to their wild-type identities')
    parser.add_argument('-noc',action='store_true',help='do not center the protein at the origin of the coordinate system')
    parser.add_argument('-ror',default='None,None',metavar='<atomselect string>,<atomselect string>',help='two comma-separated, single-quoted atomselect strings to define two groups of atoms whose centers of mass are aligned against the global z-axis')
    parser.add_argument('-v','--verbosity',action='count',help='output verbosity')
    parser.add_argument('-postscript',metavar='<name>',default='postscript.sh',help='autogenerated shell script to be run after vmd')
    parser.add_argument('-pe',metavar='<int>',default=8,type=int,help='number of processors to indicated in NAMD inputs')
    parser.add_argument('-mutationsVsLinks',metavar='M/L',default='M',help='M if mutations trump links (default behavior), L otherwise')
    parser.add_argument('-mutationsVsSSBonds',metavar='M/S',default='M',help='M if mutations trump disulfides (default behavior), or S otherwise')
    parser.add_argument('-deletionsVsLinks',metavar='D/L',default='D',help='D if deletions trump links (default behavior), L otherwise')
    parser.add_argument('-deletionsVsSSBonds',metavar='D/S',default='D',help='D if deletions trump disulfides (default behavior), or S otherwise')
    args=parser.parse_args()
    UIC=args.ignore
    IOC=args.includeonly
    if len(sys.argv)>1:
        print('Command-line arguments: '+' '.join(sys.argv[1:]))

    if args.verbosity!=None:
        print('### Vebosity level: {}'.format(args.verbosity))
    else:
       args.verbosity=0

    if args.verbosity>0:
        for k,v in vars(args).items():
            if type(v) is list and len(v)>0:
                print('-{:s} '.format(k)+' '.join([str(_) for _ in v]),end=' ')
            elif type(v) is str and len(v)>0:
                print('-{:s} {}'.format(k,v),end=' ')
            elif not type(v) is str and not type(v) is list:
                print('-{:s} {}'.format(k,str(v)),end=' ')
        print()

    Mut=MrgCmdLineAndFileContents(args.mut,args.mutfile,Mutation)
    Clv=MrgCmdLineAndFileContents(args.clv,args.clvfile,Cleavage)
    Gra=MrgCmdLineAndFileContents(args.gra,args.grafile,Graft)
    Att=MrgCmdLineAndFileContents(args.att,args.attfile,Attach)
    Uss=MrgCmdLineAndFileContents(args.ssbond,args.ssfile,SSBond)
    Uxss=MrgCmdLineAndFileContents(args.exclude_ssbond,args.exclude_ssfile,SSBond)
    Usl=MrgCmdLineAndFileContents(args.link,args.linkfile,Link)
    Del=MrgCmdLineAndFileContents(args.delete,args.deletefile,Deletion)
    Mis=MrgCmdLineAndFileContents(args.missing,args.missingfile,Missing)
    if len(args.modsfile)>0:
        for mf in args.modsfile:
            if args.verbosity>0:
                mf.report()
            Mut.extend(mf.show_type(Mutation))
            Clv.extend(mf.show_type(Cleavage))
            Gra.extend(mf.show_type(Graft))
            Att.extend(mf.show_type(Attach))
            Uss.extend(mf.show_type(SSBond))
            Usl.extend(mf.show_type(Link))
            Del.extend(mf.show_type(Deletion))
            Mis.extend(mf.show_type(Missing))
            Uxss.extend(mf.show_excl(SSBond))
    userMods={'userMutations':Mut,'userCleavages':Clv,'userGrafts':Gra,'userAttachments':Att,
               'userSSBonds':Uss,'userXSSBonds':Uxss,'userLinks':Usl,'userDeletions':Del,'userMissing':Mis,
               'fixConflicts':args.fixconflicts,'fixEngineeredMutations':args.fixengineeredmutations,'ignoreChains':UIC,
               'includeOnlyChains':IOC,'includeTerminalLoops':args.includeTerminalLoops,
               'mutationsVsLinks':args.mutationsVsLinks,'mutationsVsSSBonds':args.mutationsVsSSBonds,
               'deletionsVsLinks':args.deletionsVsLinks,'deletionsVsSSBonds':args.deletionsVsSSBonds}

    UPDBAliases=MrgCmdLineAndFileContents([' '.join(_.split(',')) for _ in args.pdbalias],args.pdbaliasfile,str)
    PDBAliases.extend(UPDBAliases)
    CTopo.extend(args.charmmtopo)
    LocTopo.extend(args.loctopo)
    StdParamFiles.extend(args.charmmparam)
    LocalParamFiles.extend(args.locparam)

    prefix=args.prefix
#    PostMod['do_loop_mc']=args.rlxloops
#    PostMod['loop_mc_params']=DictFromString(args.loopmcparams)
 #   PostMod['do_gly_mc']=args.rlxgly
 #   PostMod['gly_mc_params']=DictFromString(args.glymcparams)
    #PostMod['do_multiflex_mc']=args.rlxmc
    #ostMod['multiflex_mc_params']=DictFromString(args.rlxmcparams)
    PostMod['do_preclose_min_smd']=args.smdclose
    PostMod['preclose_params']=DictFromString(args.smdcloseparams)
    PostMod['NAMD_params']=DictFromString(args.namdparams)
    PostMod['Crot']=MrgCmdLineAndFileContents(args.crot,args.crotfile,Crot)
    if len(args.modsfile)>0:
        for mf in args.modsfile:
            PostMod['Crot'].extend(mf.show_type(Crot))
    PostMod['log_dcd_file']=args.logdcd
    PostMod['log_every']=args.logevery
    PostMod['log_save_every']=args.logsaveevery

    psfgen=args.psfgen
    PostMod['center_protein']=~(args.noc)
    if args.ror!='None,None':
        PostMod['reorient_protein']=True
        PostMod['center_protein']=True
        PostMod['reorselstr']=args.ror.split(',')
    #for k,v in PostMod.items():
    #    print('{}:'.format(k),v)

    postscriptname=args.postscript
    npe=args.pe
    #print('-pe {:d}; NAMD will use {:d} processors.'.format(args.pe,npe))
 
    PDBfiles=args.inpdb+args.incif
    Molecules=[]
    if '.cif' in PDBfiles[0]:
        Molecules.append(Molecule(cif=PDBfiles[0],userMods=userMods,requestedBiologicalAssembly=args.biological_assembly))
    elif '.pdb' in PDBfiles[0]:
        Molecules.append(Molecule(pdb=PDBfiles[0],userMods=userMods,requestedBiologicalAssembly=args.biological_assembly))
    # auxiliary PDB/mmCIF molecules -- why are these needed?
    for p in PDBfiles[1:]:
        if '.cif' in p:
            Molecules.append(Molecule(cif=p))
        elif '.pdb' in p:
            Molecules.append(Molecule(pdb=p))
    Base=Molecules[0]
    Base.summarize()
    if len(Clv)>0:
        Base.CleaveChains(Clv)

    ''' Generate the psfgen TcL script '''
    psfgen_fp=open(psfgen,'w')
    psfgen_fp.write('### This is an automatically generated psfgen input file\n')
    psfgen_fp.write('### created using pestifer.py on {} at {}\n'.format(date.today(),datetime.now().strftime('%H:%M:%S')))
    psfgen_fp.write('### pestifer.py is part of the psfgen repository\n')
    psfgen_fp.write('### github.com:cameronabrams/psfgen/scripts\n')
    psfgen_fp.write('### questions to cfa22@drexel.edu\n')
    psfgen_fp.write('### command: python3 '+' '.join(sys.argv)+'\n')
    WriteHeaders(psfgen_fp,CTopo,LocTopo,PDBAliases)
    Loops=Base.WritePsfgenInput(psfgen_fp,prefix=prefix)
    ''' PostMods are commands that operate on the PSF/PDB pair generated above, and are included in the
        TcL script for psfgen.  These are typically commands that alter coordinates, like centering, rotating, and
        adjusting dihedrals to ease future minimization.  The PSF itself is not modified.
        Regardless of whether any modifications are done or not, this will always write 
        a *_mod.pdb coordinate file '''
    post_pdb=WritePostMods(psfgen_fp,Base.psf_outfile,Base.pdb_outfile,PostMod,Loops,Base.getGlycanSegnames())
    Base.Tcl_PrependHeaderToPDB(post_pdb,psfgen_fp)
    psfgen_fp.write('exit\n')
    psfgen_fp.close()

    ''' Generate the shell script that manages all invocations of vmd and namd2 to 
        complete the system build '''
    nummin=1000
    numsteps=2000
    temperature=310
    if 'NAMD_params' in PostMod:
        p=PostMod['NAMD_params']
        nummin=p.get('nummin',nummin)
        numsteps=p.get('numsteps',numsteps)
        temperature=p.get('temperature',temperature)
    currpdb=post_pdb
    currpsf=Base.psf_outfile
    print('Run the script {} to complete the build.'.format(postscriptname))
    print('After running {}, "read CURRPSF CURRPDB CURRCFG < .tmpvar" will set those variables.'.format(postscriptname))
    print('cfapdbpyparse ends.')

    fp=open(postscriptname,'w')
    fp.write(r'#!/bin/bash'+'\n')
    fp.write('# {}: completes the build of {}\n'.format(postscriptname,currpsf))
    fp.write('''
source $PSFGEN_BASEDIR/scripts/utils.sh
nesting_level=1
TASK=1
i=1
ARGC=$#
while [ $i -le $ARGC ] ; do
  if [ "${!i}" = "-task" ]; then
    i=$((i+1))
    TASK=${!i}
  elif [ "${!i}" = "-nesting-level" ]; then
    i=$((i+1))
    nesting_level=${!i}
  else
    echo "${!i} not known. Exiting."
    exit 1
  fi
  i=$((i+1))
done
ind=`indent $nesting_level "#"`
''')
    fp.write('echo "$ind Postscript {} begins."\n'.format(postscriptname))
    fp.write('echo "$ind Completing the task-'+r'${TASK}'+' build of {}"\n'.format(currpsf))
    vmd_instructions(fp,psfgen,logname=r'psfgen${TASK}.log',msg='generates psf={} pdb={}'.format(currpsf,currpdb))
    fp.write("cat {} | sed \'1,/#### BEGIN PATCHES/d;/#### END PATCHES/,$d\' > patches.inp\n".format(psfgen))
    outname=r'postnamd${TASK}-1'
    currcfg=r'run${TASK}-1.namd'
    currlog=r'run${TASK}-1.log'
    namd_instructions(fp,currcfg,currpsf,currpdb,outname,currlog,npe=npe,
                      numminsteps=nummin,numsteps=numsteps,seed=random.randint(0,10000),
                      template='vac.namd',temperature=temperature,msg='first relaxation',
                      stdparamfiles=StdParamFiles,localparamfiles=LocalParamFiles)
    namdbin='{}.coor'.format(outname)
    currpdb='{}.pdb'.format(outname)
    vmd_instructions(fp,r'$PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl',logname=r'namdbin2pdb${TASK}-1.log',
                        args='{} {} {}'.format(currpsf,namdbin,'tmp.pdb'),msg='converting namdbin to pdb')
    fp.write('cat charmm_header.pdb tmp.pdb > {}\n'.format(currpdb))
    logname=r'ringp${TASK}-1.log'
    vmd_instructions(fp,r'$PSFGEN_BASEDIR/scripts/ringp.tcl',logname=logname,args='{} {}'.format(currpsf,currpdb),msg='checking for pierced rings')
    fp.write('npiercings=`grep -c pierces {}`\n'.format(logname))
    fp.write(r'if [[ $npiercings -gt 0 ]]; then'+'\n')
    fp.write(r'  echo "Error: There are $npiercings piercings in '+'{}"\n'.format(currpdb))
    fp.write('  grep pierces {}\n'.format(logname))
    fp.write('  echo "$ind Change your relaxation parameters and try again."\n')
    fp.write('  exit 1\n')
    fp.write('fi\n')
    if 'do_preclose_min_smd' in PostMod and PostMod['do_preclose_min_smd']:
        temperature_close=400
        target_numsteps=20000
        if 'preclose_params' in PostMod:
            p=PostMod['preclose_params']
            temperature_close=p.get('temperature_close',temperature_close)
            target_numsteps=p.get('target_numsteps',target_numsteps)
        fp.write('cat > close_these.inp << EOF\n')
        for l in sorted(Loops, key=lambda x: len(x.residues)):
            if (l.term=='None' and len(l.residues)>2):
                fp.write(l.input_str())
        fp.write('EOF\n')
        fp.write('cp close_these.inp close_these.inp.bak\n')
        # measures to find the initial distances; generated fixed.pdb to fix the N atoms 
        logname=r'close${TASK}.log'
        args='{} {} close_these.inp fixed.pdb'.format(currpsf,currpdb)
        vmd_instructions(fp,r'$PSFGEN_BASEDIR/scripts/measure_bonds.tcl',logname=logname,args=args,msg='creating closing input')
        fp.write('if [ -f cv.inp ]; then rm cv.inp; fi\n')
        fp.write('touch cv.inp\n')
        fp.write('while IFS=" " read -r C L R B; do\n')
        fp.write(r'  cat $PSFGEN_BASEDIR/templates/cv-template.in | sed s/%C%/$C/g |')
        fp.write(r'  sed s/%NAME%/${C}${L}/g | sed s/%I%/$L/g | sed s/%J%/$R/g | sed s/%R0%/$B/g |')
        fp.write('  sed s/%TARGETNUMSTEPS%/{}/ >> cv.inp ;\n'.format(target_numsteps))
        fp.write('done < close_these.inp\n')
        outname=r'postnamd${TASK}-2'
        currcfg=r'run${TASK}-2.namd'
        currlog=r'run${TASK}-2.log'
        extras=['fixedatoms on','fixedatomsfile fixed.pdb','fixedatomscol B','colvars on','colvarsconfig cv.inp']
        namd_instructions(fp,currcfg,currpsf,currpdb,outname,currlog,npe=npe,
                      numminsteps=0,numsteps=int(1.5*target_numsteps),seed=random.randint(0,10000),
                      template='vac.namd',temperature=temperature_close,extras=extras,msg='closing',
                      stdparamfiles=StdParamFiles,localparamfiles=LocalParamFiles)
        namdbin='{}.coor'.format(outname)
        currpdb='{}.pdb'.format(outname)
        vmd_instructions(fp,r'$PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl',logname=r'namdbin2pdb${TASK}-1.log',
                        args='{} {} {}'.format(currpsf,namdbin,'tmp.pdb'),msg='converting namdbin to pdb')
        fp.write('cat charmm_header.pdb tmp.pdb > {}\n'.format(currpdb))
#        fp.write('cat > closure_patches.inp << EOF\n')
        fp.write('cat > the_closing_patches.inp << EOF\n')
        for l in sorted(Loops, key=lambda x: len(x.residues)):
            if (l.term=='None' and len(l.residues)>2):
                #fp.write('# will try to close bond between {} and {} on chain {}...\n'.format(l.residues[-1].resseqnum,l.nextfragntermresid,l.replica_chainID))
                fp.write(l.heal_str())
        fp.write('EOF\n')
        tfp=open('topologies.inp','w')
        CommonPSFGENheader(tfp,CTopo,LocTopo)
        tfp.close()
        fp.write('cat $PSFGEN_BASEDIR/scripts/loop_closure.tcl | sed "/#### LIGATION LIST STARTS/r the_closing_patches.inp"')
        fp.write(' | sed "/#### TOPOLOGY FILE LIST STARTS/r topologies.inp" > do_the_closures.tcl\n')
        newpsf='ligated.psf'
        newpdb='ligated.pdb'
        vmd_instructions(fp,'do_the_closures.tcl',args='{} {} {} {}'.format(currpsf,currpdb,newpsf,newpdb),logname=r'ligations${TASK}.log')
        currpsf=newpsf
        currpdb=newpdb
        currcfg=r'run${TASK}-3.namd'
        currlog=r'run${TASK}-3.log'
        outname=r'postnamd${TASK}-3'
        namd_instructions(fp,currcfg,currpsf,currpdb,outname,currlog,npe=npe,
                      numminsteps=nummin,numsteps=numsteps,seed=random.randint(0,10000),
                      template='vac.namd',temperature=temperature,msg='minimization of ligated peptide bonds',
                      stdparamfiles=StdParamFiles,localparamfiles=LocalParamFiles)
        namdbin='{}.coor'.format(outname)
        currpdb='{}.pdb'.format(outname)
        vmd_instructions(fp,r'$PSFGEN_BASEDIR/scripts/namdbin2pdb.tcl',logname=r'namdbin2pdb${TASK}-1.log',
                        args='{} {} {}'.format(currpsf,namdbin,'tmp.pdb'),msg='converting namdbin to pdb')
        fp.write('cat charmm_header.pdb tmp.pdb > {}\n'.format(currpdb))
 
        logname=r'ringp${TASK}-2.log'
        vmd_instructions(fp,r'$PSFGEN_BASEDIR/scripts/ringp.tcl',logname=logname,args='{} {}'.format(currpsf,currpdb),msg='checking for pierced rings')
        fp.write('npiercings=`grep -c pierces {}`\n'.format(logname))
        fp.write(r'if [[ $npiercings -gt 0 ]]; then'+'\n')
        fp.write(r'  echo "Error: There are $npiercings piercings in '+'{}"\n'.format(currpdb))
        fp.write('  grep pierces {}\n'.format(logname))
        fp.write('  echo "Change your relaxation parameters and try again."\n')
        fp.write('  exit 1\n')
        fp.write('fi\n')

    fp.write('echo {} {} {} > .tmpvar\n'.format(currpsf,currpdb,currcfg))
    fp.write('echo "$ind Postscript {} finishes."\n'.format(postscriptname))
    fp.write('exit 0')
    fp.close()
    os.system('chmod 744 {}'.format(postscriptname))

def solvate():
    pass

def equilibrate():
    pass

def pack_for_production():
    pass

def do_topogromacs():
    pass
