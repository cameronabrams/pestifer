from time import sleep
from pestifer.util.logparsers import PackmolLog, NAMDLog, get_toeol, get_tokens, get_values
import logging
logger=logging.getLogger(__name__)

def test_get_toeol():
    flag='A_LINE> '
    bytes='A_LINE> 1.23e9'
    res=get_toeol(flag,bytes)
    assert res==None
    bytes='A_LINE> 1.23e9\n'
    res=get_toeol(flag,bytes)
    assert res=='1.23e9'
    bytes='A_LINE> 1.23e9 3.45 6.78\n'
    res=get_toeol(flag,bytes)
    assert res=='1.23e9 3.45 6.78'

def test_get_tokens():
    flag='A_LINE> '
    line='A_LINE> 1.23e9 3.45 6.78\n'
    res=get_tokens(flag,line)
    assert res==['1.23e9', '3.45', '6.78']

def test_get_values():
    flag='A_LINE> '
    line='A_LINE> 1.23e9 3.45 6.78\n'
    res=get_values(flag,line)
    assert res==[1.23e9, 3.45, 6.78]

def test_packmol_log_static():
    l=PackmolLog()
    with open('packmol/a.pm-log','r') as f:
        msg=f.read()
    l.update(msg)
    assert len(l.processed_banners)==9
    assert l.metadata['version']=='20.15.1'
    assert l.metadata['initial_objective_function']==308065.45399671065
    assert len(l.metadata['molecule_types_packed'])==4
    assert l.metadata['molecule_types_packed'][0]=='1'
    assert l.metadata['molecule_types_packed'][1]=='2'
    assert l.metadata['molecule_types_packed'][2]=='3'
    assert l.metadata['molecule_types_packed'][3]=='4'
    assert l.metadata['pbc_boxsize'][0]==70.71
    assert l.metadata['pbc_boxsize'][1]==70.71
    assert l.metadata['pbc_boxsize'][2]==95.99
    assert l.metadata['coordinate_file_types']=='pdb'
    assert l.metadata['pbc_reference_box'][0]==0
    assert l.metadata['pbc_reference_box'][1]==0
    assert l.metadata['pbc_reference_box'][2]==-48.00
    assert l.metadata['pbc_reference_box'][3]==70.71
    assert l.metadata['pbc_reference_box'][4]==70.71
    assert l.metadata['pbc_reference_box'][5]==48.00
    assert l.metadata['seed']==27021972
    assert l.metadata['url']=='http://m3g.iqm.unicamp.br/packmol'
    assert l.metadata['output_file']=='patch.pdb'
    assert l.metadata['coordinate_files'][0]=='POPC-00.pdb'
    assert l.metadata['coordinate_files'][1]=='POPC-00.pdb'
    assert l.metadata['coordinate_files'][2]=='TIP3.pdb'
    assert l.metadata['coordinate_files'][3]=='TIP3.pdb'
    assert len(l.metadata['structures'])==4
    assert l.metadata['structures'][0]['maximum_internal_distance']==31.769849621929282
    assert l.metadata['structures'][1]['maximum_internal_distance']==31.769849621929282
    assert l.metadata['structures'][2]['maximum_internal_distance']==1.5175773456400830
    assert l.metadata['structures'][3]['maximum_internal_distance']==1.5175773456400830
    assert l.metadata['structures'][0]['before_adjusting']['function_value_before_moving_molecules'][0]==6132.6785289414929
    assert l.metadata['structures'][0]['before_adjusting']['function_value_after_moving_molecules'][0]==4999.3267934110645
    assert len(l.metadata['structures'][0]['before_adjusting']['function_value_before_moving_molecules'])==1
    assert len(l.metadata['structures'][0]['before_adjusting']['function_value_after_moving_molecules'])==1
    assert l.metadata['structures'][0]['before_adjusting']['restraint_only_function_value']==1.1265504922994881E-003
    assert l.metadata['structures'][0]['before_adjusting']['maximum_violation_of_the_restraints']==5.1528308368157096E-004
    assert l.metadata['structures'][0]['after_adjusting']['restraint_only_function_value']==1.2707914952250293E-002
    assert l.metadata['structures'][0]['after_adjusting']['maximum_violation_of_the_restraints']==6.9373223144957559E-003
    assert l.metadata['structures'][1]['before_adjusting']['function_value_before_moving_molecules'][0]==556490.28355895542
    assert l.metadata['structures'][1]['before_adjusting']['function_value_after_moving_molecules'][0]==518979.69846450235
    assert l.metadata['structures'][1]['before_adjusting']['function_value_before_moving_molecules'][1]==407234.09803705831
    assert l.metadata['structures'][1]['before_adjusting']['function_value_after_moving_molecules'][1]==375844.82748019119 
    assert l.metadata['structures'][1]['before_adjusting']['function_value_before_moving_molecules'][2]==11312.964715378288
    assert l.metadata['structures'][1]['before_adjusting']['function_value_after_moving_molecules'][2]==6056.6893181096848 
    assert len(l.metadata['structures'][1]['before_adjusting']['function_value_before_moving_molecules'])==3

    assert len(l.gencan[1])==101

    assert 'gencan_success' not in l.metadata['structures'][0]
    assert 'gencan_success' not in l.metadata['structures'][1]
    assert 'gencan_success' in l.metadata['structures'][2]
    assert 'gencan_success' in l.metadata['structures'][3]

    assert len(l.molmoves)==57

    l.finalize()

def test_packmol_log_dynamic():
    l=PackmolLog()
    with open('packmol/a.pm-log','r') as f:
        msg=f.read()
    msg_len=len(msg)
    chunk_size=msg_len//100
    p=0
    prog=[0.0]
    while p<msg_len:
        l.update(msg[p:p+chunk_size])
        p+=chunk_size
        if l.measure_progress()>prog[-1]:
            prog.append(l.measure_progress())
        sleep(0.1)
    l.finalize()
    assert len(prog)==97

def test_namd_log_static():
    l=NAMDLog()
    with open('namd/test_namd.testlog','r') as f:
        msg=f.read()
    l.update(msg)
    assert len(l.processed_line_idx)==607
    assert l.metadata['timestep']==2.0
    assert l.metadata['atom_density']==0.102335
    assert l.metadata['total_mass']==258901
    assert l.metadata['number_of_steps']==9000
    assert l.metadata['number_of_atoms']==44200
    assert l.metadata['total_charge']==1.78814e-05
    assert l.metadata['running_for']==8000
    assert l.measure_progress()==1.0
    l.finalize()

def test_namd_log_dynamic():
    l=NAMDLog()
    with open('namd/test_namd.testlog','r') as f:
        msg=f.read()
    msg_len=len(msg)
    chunk_size=msg_len//100
    p=0
    while p<msg_len:
        l.update(msg[p:p+chunk_size])
        p+=chunk_size
        sleep(0.1)
        logger.debug(f'progress {l.measure_progress()}')
    assert l.success()
    l.finalize()