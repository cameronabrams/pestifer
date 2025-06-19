# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import numpy as np
import os

from pidibble.pdbparse import PDBParser

from ..basetask import BaseTask

from ..util.util import cell_from_xsc, cell_to_xsc

class SolvateTask(BaseTask):
    yaml_header='solvate'

    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        self.next_basename()
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        xsc=self.statevars.get('xsc',None)
        use_minmax=True
        if xsc is not None:
            box,origin=cell_from_xsc(xsc)
            if box is not None and origin is not None:
                use_minmax=False
                basisvec=box.diagonal()
                LL=origin-0.5*basisvec
                UR=origin+0.5*basisvec
        if use_minmax:
            p=PDBParser(PDBcode=os.path.splitext(pdb)[0]).parse()
            x=np.array([a.x for a in p.parsed['ATOM']])
            y=np.array([a.y for a in p.parsed['ATOM']])
            z=np.array([a.z for a in p.parsed['ATOM']])
            minmax=np.array([[x.min(),y.min(),z.min()],[x.max(),y.max(),z.max()]])
            spans=minmax[1]-minmax[0]
            maxspan=spans.max()
            cubic=self.specs.get('cubic',False)
            sympad=np.array([0.,0.,0.])
            if cubic:
                sympad=0.5*(maxspan-spans)
            pad=self.specs.get('pad',10.0)
            LL=minmax[0]-pad-sympad
            UR=minmax[1]+pad+sympad
            basisvec=UR-LL
            origin=0.5*(UR+LL)
            cell_to_xsc(np.diag(basisvec),origin,f'{self.basename}.xsc')

        ll_tcl=r'{ '+' '.join([str(_) for _ in LL.tolist()])+r' }'
        ur_tcl=r'{ '+' '.join([str(_) for _ in UR.tolist()])+r' }'
        box_tcl=r'{ '+ll_tcl+' '+ur_tcl+r' }'
        vt=self.writers['vmd']
        vt.newscript(self.basename)
        vt.addline( 'package require solvate')
        vt.addline( 'package require autoionize')
        vt.addline( 'psfcontext mixedcase')
        vt.addline(f'mol new {psf}')
        vt.addline(f'mol addfile {pdb} waitfor all')
        vt.addline(f'solvate {psf} {pdb} -minmax {box_tcl} -o {self.basename}_solv')
        ai_args=[]
        sc=self.specs.get('salt_con',None)
        if sc is None:
            ai_args.append('-neutralize')
        else:
            ai_args.append(f'-sc {sc}')
        cation=self.specs.get('cation',None)
        if cation is not None:
            ai_args.append(f'-cation {cation}')
        anion=self.specs.get('anion',None)
        if anion is not None:
            ai_args.append(f'-anion {anion}')
        vt.addline(f'autoionize -psf {self.basename}_solv.psf -pdb {self.basename}_solv.pdb {" ".join(ai_args)} -o {self.basename}')
        vt.writescript()
        self.result=vt.runscript(progress_title='solvate')
        if self.result!=0:
            return super().do()
        self.save_state(exts=['psf','pdb','xsc'])
        self.log_message('complete')
        return super().do()
