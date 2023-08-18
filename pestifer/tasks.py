from .basemod import BaseMod
import logging
logger=logging.getLogger(__name__)
import shutil
import os
class Task(BaseMod):
    req_attr=BaseMod.req_attr+['owner','prior','next','specs','taskname']
    yaml_header='generic_task'
    default_specs={}
    exts=['.psf','.pdb']
    def __init__(self,input_dict,owner=None):
        assert len(input_dict)==1
        taskname=list(input_dict.keys())[0]
        specval=list(input_dict.values())[0]
        specs={} if not specval else specval
        for k,v in self.__class__.default_specs.items():
            if not k in specs:
                specs[k]=v
        input_dict.update({
            'owner':owner,
            'prior':None,
            'next':None,
            'specs':specs,
            'taskname':taskname
        })
        super().__init__(input_dict)

    def priors(self):
        prior_basename=self.prior.basename
        return [f'{prior_basename}{ext}' for ext in self.__class__.exts]
    
    def pass_thru(self):
        for fn in self.priors():
            if os.path.exists(fn):
                bn,ext=os.path.splitext(fn)
                shutil.copyfile(fn,f'{self.basename}{ext}')
            else:
                logger.debug(f'Task {self.taskname}: prior task\'s {fn} not found')

class PsfgenTask(Task):
    yaml_header='psfgen'
    opt_attr=Task.opt_attr+[yaml_header]
    default_specs={'cleanup':False}
    def do(self):
        mod_dict=self.specs.get('mods',{})
        base_molecule=self.owner.base_molecule
        pg=self.owner.psfgen
        pg.newscript(self.basename)
        pg.topo_aliases()
        pg.set_molecule(base_molecule)
        pg.describe_molecule(base_molecule,mod_dict)
        pg.transform_postmods()
        pg.endscript()
        pg.writescript()
        pg.runscript()
        pg.cleanup(cleanup=self.specs['cleanup'])

class LayloopsTask(Task):
    yaml_header='layloops'
    opt_attr=Task.opt_attr+[yaml_header]
    default_specs={'cycles':100}
    def do(self):
        logger.debug(f'layloops: specs {self.specs}')
        # self.pass_thru()
        psf,pdb=self.priors()
        vt=self.owner.vmdtcl
        mol=self.owner.base_molecule
        logger.debug(f'layloops: inputs: {psf} {pdb}')
        vt.newscript(self.basename)
        vt.load_psf_pdb(psf,pdb,new_molid_varname='mLL')
        cycles=self.specs.get('cycles',100)
        sac_n=self.specs.get('min_loop_length',4)
        ba=mol.active_biological_assembly
        au=mol.asymmetric_unit
        for S in au.Segments:
            chainID=S.chainID
            if S.segtype=='PROTEIN':
                for b in S.subsegments:
                    if b.state=='MISSING':
                        if (b.bounds[1]-b.bounds[0])>(sac_n-1):
                            reslist=[f'{r.resseqnum}{r.insertion}' for r in S.residues[b.bounds[0]:b.bounds[1]+1]]
                            tcllist='[list '+' '.join(reslist)+']'
                            for biomt in ba.biomt:
                                act_chainID=biomt.chainIDmap[chainID]
                                vt.addline(f'lay_loop $mLL {act_chainID} {tcllist} {cycles}')
        vt.write_pdb(self.basename,'mLL')
        vt.endscript()
        vt.writescript()
        vt.runscript()
        logger.debug(f'layloops is done')

class SolvateTask(Task):
    yaml_header='solvate'
    opt_attr=Task.opt_attr+[yaml_header]
    def do(self):
        logger.debug(f'solvate: specs {self.specs}')

class SMDcloseTask(Task):
    yaml_header='smdclose'
    opt_attr=Task.opt_attr+[yaml_header]
    def do(self):
        logger.debug(f'smdclose: specs {self.specs}')

class RelaxTask(Task):
    yaml_header='relax'
    opt_attr=Task.opt_attr+[yaml_header]
    def do(self):
        logger.debug(f'relax: specs {self.specs}')

    # def do_relax_series(self,task):
    #     logger.debug(f'do_relax_series: specs {task.specs}')
    #     pass

class TerminateStepTask(Task):
    yaml_header='terminate_step'
    opt_attr=Task.opt_attr+[yaml_header]
    def do(self):
        self.pass_thru()

