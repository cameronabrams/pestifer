"""

.. module:: sel
   :synopsis: functions for generating TcL/VMD atomselect commands
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
# TODO: adapt these into the scriptwriters module

def backup(selname):
    ret=[]
    attr=['chain','x','y','z','resid','resname','name']
    for a in attr:
       ret.append('set {}_orig_{} [${} get {}]'.format(selname,a,selname,a))
    return '\n'.join(ret)

def restore(selname):
    ret=[]
    attr=['chain','x','y','z','resid','resname','name']
    for a in attr:
       ret.append('${} set {} ${}_orig_{}'.format(selname,a,selname,a))
    return '\n'.join(ret)

def residshift(selname,shift):
    ret=[]
    ret.append('set new_resid [list]')
    ret.append(r'foreach oldresid [${} get resid] {{'.format(selname))
    ret.append('     lappend new_resid [expr $oldresid + {:d}]'.format(shift))
    ret.append(r'}')
    ret.append('${} set resid $new_resid'.format(selname))
    return '\n'.join(ret)

''' charmm_namify converts commonly found atom and residue names in PDB files to their 
    appropriate charmm names -- this is only used for NON-PROTEIN SEGMENTS. '''
def charmm_namify(selname,iswater=False):
    ret=[]
    ret.append('set new_resname [list]')
    ret.append(r'foreach r [${} get resname] {{'.format(selname))
    ret.append(r'   if { [ info exists RESDICT($r) ] } {')
    ret.append('      lappend new_resname $RESDICT($r)')
    ret.append(r'   } else {'+'')
    ret.append('      lappend new_resname $r')
    ret.append('   }')
    ret.append('}')
    ret.append('set new_name [list]')
    ret.append(r'foreach r [${} get name] {{'.format(selname))
    ret.append(r'   if { [ info exists ANAMEDICT($r) ] } {')
    ret.append('      lappend new_name $ANAMEDICT($r)')
    ret.append(r'   } else {')
    ret.append('      lappend new_name $r')
    ret.append('   }')
    ret.append('}')
    ret.append('${} set resname $new_resname').format(selname)
    ret.append('${} set name $new_name').format(selname)
    if iswater=='WATER':
        ret.append('${} set name OH2').format(selname)
    return '\n'.join(ret)
