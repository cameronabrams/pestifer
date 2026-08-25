from time import sleep
from pestifer.logparsers import NAMDLogParser
from pestifer.logparsers.logparser import get_toeol, get_tokens, get_values
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

def test_namd_log_static():
    l=NAMDLogParser()
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

def test_namd_log_static_incomplete():
    l=NAMDLogParser()
    with open('namd/test_namd-incomplete.testlog','r') as f:
        msg=f.read()
    l.update(msg)
    assert len(l.processed_line_idx)==3198
    assert l.metadata['timestep']==2.0
    assert l.time_series_data['restart'][-1]==120000
    assert l.time_series_data['performance'][-1]['ns_per_day']==23.2433
    assert l.time_series_data['timing'][-1]['cpu_time']==786.403
    l.finalize()

def test_namd_log_dynamic():
    l=NAMDLogParser()
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
def test_namd_energy_line_garbled_interleave():
    """A multi-rank (srun) run can interleave NAMD stdout, producing a malformed
    ENERGY line (e.g. a token like '0LINE'). The parser must skip it, not crash."""
    import os
    l = NAMDLogParser()
    nfields = len(l.default_etitle_nvt)  # NVT energy line width
    good = 'ENERGY: ' + ' '.join(str(i) for i in range(nfields)) + os.linesep
    # establishes etitle and one good sample
    assert l.process_line(good) == 0
    assert len(l.time_series_data['energy']) == 1
    # one field is corrupted by interleaved output; same token count
    bad_tokens = [str(i) for i in range(nfields)]
    bad_tokens[13] = '0LINE'
    bad = 'ENERGY: ' + ' '.join(bad_tokens) + os.linesep
    # must not raise, and must not append a bogus sample
    assert l.process_line(bad) == 0
    assert len(l.time_series_data['energy']) == 1

def test_namd_density_column_is_g_per_cc():
    """The DENSITY column must agree with NAMD's own startup MASS DENSITY.

    Pestifer computes this column itself -- NAMD's ETITLE has no DENSITY field -- and for a long
    time stored it as the raw amu/A^3 quotient, a factor of AMU_PER_A3_TO_G_PER_CC away from both
    NAMD's figure and `util.density_convergence.volume_to_density()`, which the convergence gates
    use.  The two paths disagreed silently: a plotted density read 0.62 where the gate beside it
    read 1.03.  NAMD's own line is the oracle here, so this pins the units without hard-coding a
    number that could drift with the fixture.
    """
    from pestifer.util.densityprofile import AMU_PER_A3_TO_G_PER_CC
    l = NAMDLogParser()
    with open('namd/test_namd-incomplete.testlog', 'r') as f:
        l.update(f.read())
    l.finalize()
    density = l.dataframes['energy']['DENSITY']
    namd_value = l.metadata['mass_density']          # NAMD's startup "MASS DENSITY = ... g/cm^3"
    assert abs(density.iloc[0] - namd_value) / namd_value < 1e-3
    # and explicitly not the raw quotient it used to be
    raw = namd_value / AMU_PER_A3_TO_G_PER_CC
    assert abs(density.iloc[0] - raw) / namd_value > 0.3
