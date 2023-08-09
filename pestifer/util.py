def vmd_instructions(fp,script,logname='tmp.log',args='',msg=''):
    fp.write('echo "VMD) script={} log={} msg: {}"\n'.format(script,logname,msg))
    if args!='':
        fp.write(r'$VMD -dispdev text -e '+script+r' -args '+args+r' > '+logname+' 2>&1\n')
    else:
        fp.write(r'$VMD -dispdev text -e '+script+r' > '+logname+' 2>&1\n')
    fp.write('if [ $? -ne 0 ]; then\n')
    fp.write('   echo "VMD failed.  Check the log file {}. Exiting."\n'.format(logname))
    fp.write('   exit 1\n')
    fp.write('fi\n')

def namd_instructions(fp,cfgname,psf,coor,outname,logname,
                      npe=8,numminsteps=0,numsteps=0,seed=0,template='vac.namd',
                      temperature=310,extras=[],msg='',stdparamfiles=[],localparamfiles=[],
                      stdcharmmdir=r'\$env(HOME)/charmm/toppar',
                      localcharmmdir=r'\$env(PSFGEN_BASEDIR)/charmm'):
    fp.write('cat $PSFGEN_BASEDIR/templates/{}'.format(template))
    fp.write('  | sed s/%OUT%/{}/g'.format(outname))
    fp.write('  | sed s/%NUMMIN%/{}/'.format(numminsteps))
    fp.write('  | sed s/%NUMSTEPS%/{}/'.format(numsteps))
    fp.write('  | sed s/%SEED%/{}/g'.format(seed))
    fp.write('  | sed s/%TEMPERATURE%/{}/g'.format(temperature))
    fp.write('  | sed "/#### SYSTEM CONFIGURATION FILES END/i structure {}"'.format(psf))
    fp.write('  | sed "/#### SYSTEM CONFIGURATION FILES END/i coordinates {}"'.format(coor))
    sentinelline='#### PARAMETER FILES END'
    for st in stdparamfiles:
        fp.write(' | sed "/{}/i parameters {}/{}" '.format(sentinelline,stdcharmmdir,st))
    for st in localparamfiles:
        fp.write(' | sed "/{}/i parameters {}/{}" '.format(sentinelline,localcharmmdir,st))
    sentinelline='#### EXTRAS END'
    for ex in extras:
        fp.write('  | sed "/'+sentinelline+'/i '+ex+'" ')
    fp.write(' > {}\n'.format(cfgname))
    namdp='+p{:d}'.format(npe)
    fp.write('echo "NAMD2) config={} log={} outputname={} msg={}"\n'.format(cfgname,logname,outname,msg))
    fp.write(r'$CHARMRUN '+namdp+r' $NAMD2 '+cfgname+r' > '+logname+'\n')
    fp.write('if [ $? -ne 0 ]; then\n')
    fp.write('   echo "NAMD failed.  Check log file {}. Exiting."\n'.format(logname))
    fp.write('   exit 1\n')
    fp.write('fi\n')
    
def special_update(dict1,dict2):
    for k,v in dict2.items():
        ov=dict1.get(k,None)
        if not ov:
            dict1[k]=v
        else:
            if type(v)==list and type(ov)==list:
                for nv in v:
                    if not nv in ov:
                        ov.append(nv)
            elif type(v)==dict and type(ov)==dict:
                ov.update(v)
            else:
                dict1[k]=v # overwrite
    return dict1

def isidentity(t):
    if t[0][0]==1.0 and t[1][1]==1.0 and t[2][2]==1.0:
        return True
    else:
        return False

def reduce_intlist(L):
    """reduce_intlist generates a "reduced-byte" representation of a list of integers by collapsing runs of adjacent integers into 'i to j' format. Example:

    [1,2,3,4,5,7,8,9,10,12] -> '1 to 5 7 to 10 12'

    :param L: list of integers
    :type L: list
    :return: string of reduced-byte representation
    :rtype: string
    """
    if not L:
        return ''
    ret=f'{L[0]}'
    inrun=False
    for l,r in zip(L[1:-1],L[2:]):
        adj=(r-l)==1
        if adj and not inrun:
            inrun=True
            ret+=f' to '
        elif not adj and inrun:
            ret+=f'{l} {r}'
            inrun=False
        elif not inrun:
            ret+=f' {l}'
    if inrun:
        ret+=f'{r}'
    return ret