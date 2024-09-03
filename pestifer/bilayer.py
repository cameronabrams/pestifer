# Author: Cameron F. Abrams, <cfa2@drexel.edu>
import logging
import os
import shutil
import yaml

import numpy as np
import pandas as pd

from scipy.constants import physical_constants, Avogadro
from scipy.spatial import ConvexHull

from .charmmtop import CharmmResiDatabase
from .coord import coorddf_from_pdb
from .psf import PSFContents
from .scriptwriters import PackmolInputWriter
from .tasks import BaseTask
from .util import cell_to_xsc

logger=logging.getLogger(__name__)

g_per_amu=physical_constants['atomic mass constant'][0]*1000
A_per_cm=1.e8
A3_per_cm3=A_per_cm**3
cm3_per_A3=1.0/A3_per_cm3
n_per_mol=Avogadro
def molec_n(MW_g,density_gcc,volume_A3):
    return np.floor(density_gcc/MW_g*cm3_per_A3*n_per_mol*volume_A3)

class BilayerEmbedTask(BaseTask):
    """ A class for handling embedding proteins into bilayers
    
    Attributes
    ----------
    yaml_header(str) 

    Methods
    -------
    do(): 
        Based on specs, writes a packmol input file to generate a membrane-embedded protein and
        then runs packmol

    """
    yaml_header='bilayer'
    def __init__(self,input_dict,taskname,config,writers,prior):
        super().__init__(input_dict,taskname,config,writers,prior)
        self.progress=config.progress
        self.lipid_pdb_path=os.path.join(config.charmmff_pdb_path,'lipid')
        self.water_ion_pdb_path=os.path.join(config.charmmff_pdb_path,'water_ions')
        assert os.path.exists(self.lipid_pdb_path),f'No lipid PDB database found -- bad installation!'

    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        RDB=CharmmResiDatabase()
        RDB.add_stream('lipid')
        RDB.add_topology('toppar_all36_moreions.str',streamnameoverride='water_ions')

        psf=self.statevars.get('psf',None)
        pro_psc=PSFContents(psf)
        pro_charge=pro_psc.get_charge()
        pdb=self.statevars.get('pdb',None)
        logger.debug(f'BilayerEmbedTask will use psf {psf} and pdb {pdb} as inputs')

        self.next_basename('bilayer')

        lipid_specstring=self.specs.get('lipids','POPC')
        ratio_specstring=self.specs.get('mole_fractions','1.0')
        solvent_specstring=self.specs.get('solvents','TIP3')
        solv_molfrac_specstring=self.specs.get('solvent_mole_fractions','1.0')
        cation_name=self.specs.get('cation','POT')
        cation_q=RDB['water_ions'][cation_name].charge
        anion_name=self.specs.get('anion','CLA')
        anion_q=RDB['water_ions'][anion_name].charge
        assert cation_name in RDB['water_ions'],f'Cation {cation_name} not found in available CHARMM topologies'
        assert anion_name in RDB['water_ions'],f'Anion {anion_name} not found in available CHARMM topologies'
        shutil.copy(os.path.join(self.water_ion_pdb_path,f'{cation_name}.pdb'),'./')
        shutil.copy(os.path.join(self.water_ion_pdb_path,f'{anion_name}.pdb'),'./')
        salt_con_M=self.specs.get('salt_con',0.0)

        SAPL=self.specs.get('SAPL',60.0)
        seed=self.specs.get('seed',27021972)
        tolerance=self.specs.get('tolerance',2.0)
        solution_gcc=self.specs.get('solution_gcc',1.0)
        prot_vol_factor=self.specs.get('prot_volume_factor',1.15)
        leaflet_thickness=self.specs.get('leaflet_thickness',22.0)
        nloop=self.specs.get('nloop',200)
        nloop_all=self.specs.get('nloop_all',200)
        boxdim=np.array(list(map(float,self.specs.get('dims',[]))))

        slices={'LOWER-CHAMBER':{},'LOWER-LEAFLET':{},'UPPER-LEAFLET':{},'UPPER-CHAMBER':{}}

        if 'embed' in self.specs:
            # use the protein and embedding specs to size the box
            embed_specs=self.specs['embed']
            prot_rad_scal=embed_specs.get('protein_radius_scaling',1.0)
            xydist=embed_specs.get('xydist',0.0)
            zdist=embed_specs.get('zdist',0.0)
            self.next_basename('embed')
            self.membrane_embed(pdb,[embed_specs['z_head_group'],embed_specs['z_tail_group'],embed_specs['z_ref_group']['text']],embed_specs['z_ref_group']['z_value'],outpdb=f'{self.basename}.pdb')
            self.save_state(exts=['pdb'])
            pdb=self.statevars.get('pdb',None)
            with open(f'{self.basename}-results.yaml','r') as f:
                embed_results=yaml.safe_load(f)
            protein_rad=float(embed_results["maxrad"])
            protein_ll=np.array(list(map(float,embed_results["lower_left"])))
            protein_ur=np.array(list(map(float,embed_results["upper_right"])))
            protein_com=np.array(list(map(float,embed_results["com"])))
            midplane=0.0
            # logger.debug(f'protein ll {protein_ll} ur {protein_ur} com {protein_com}')
            box_ll=protein_ll-np.array([xydist,xydist,zdist])
            box_ur=protein_ur+np.array([xydist,xydist,zdist])
            
            shift=-box_ll
            box_ur+=shift
            box_ll+=shift
            midplane+=shift[2]

            origin=0.5*(box_ur+box_ll)
            boxdim=box_ur-box_ll
            mem_area=boxdim[1]*boxdim[2]
            protein_xyarea=np.pi*(protein_rad*prot_rad_scal)**2
            logger.debug(f'protein radius {protein_rad:.3f} xyarea {protein_xyarea:.3f}')
            mem_area-=protein_xyarea
            slices['LOWER-CHAMBER']['z-lo']= box_ll[2]
            slices['LOWER-CHAMBER']['z-hi']= midplane-leaflet_thickness
            slices['LOWER-LEAFLET']['z-lo']= midplane-leaflet_thickness
            slices['LOWER-LEAFLET']['z-hi']= midplane
            slices['UPPER-LEAFLET']['z-lo']= midplane
            slices['UPPER-LEAFLET']['z-hi']= midplane+leaflet_thickness
            slices['UPPER-CHAMBER']['z-lo']= midplane+leaflet_thickness
            slices['UPPER-CHAMBER']['z-hi']= box_ur[2]
        else:
            box_ll=-0.5*boxdim
            box_ur=0.5*boxdim
            origin=0.5*(box_ur+box_ll)
            mem_area=boxdim[1]*boxdim[2]
            slices['LOWER-CHAMBER']['z-lo']= box_ll[2]
            slices['LOWER-CHAMBER']['z-hi']= origin[2]-leaflet_thickness
            slices['LOWER-LEAFLET']['z-lo']= origin[2]-leaflet_thickness
            slices['LOWER-LEAFLET']['z-hi']= origin[2]
            slices['UPPER-LEAFLET']['z-lo']= origin[2]
            slices['UPPER-LEAFLET']['z-hi']= origin[2]+leaflet_thickness
            slices['UPPER-CHAMBER']['z-lo']= origin[2]+leaflet_thickness
            slices['UPPER-CHAMBER']['z-hi']= box_ur[2]

        span=box_ur-box_ll
        boxV=span.prod()
        logger.debug(f'box corners ll {box_ll} ur {box_ur}')
        logger.debug(f'membrane area {mem_area}')
        logger.debug(f'box volume {boxV:.3f}')

        cell_to_xsc(np.diag(boxdim),origin,f'{self.basename}.xsc')
        self.save_state(exts=['xsc'])

        solvent_names=solvent_specstring.split(':')
        solvent_molfracs=np.array(list(map(float,solv_molfrac_specstring.split(':'))))
        assert len(solvent_names)==len(solvent_molfracs),f'Solvent specs are not congruent'
        solvent={}
        for sn,sm in zip(solvent_names,solvent_molfracs):
            assert sn in RDB['water_ions'],f'solvent {sn} is not found in the available CHARMM topologies'
            solvent[sn]={}
            solvent[sn]['mol-frac']=sm
            solvent_pdb_path=os.path.join(self.water_ion_pdb_path,f'{sn}.pdb')
            assert os.path.exists(solvent_pdb_path),f'{sn}.pdb is not found'
            shutil.copy(solvent_pdb_path,'./')
            solvent[sn]['mass']=RDB['water_ions'][sn].mass()

        leaf_lipspec=lipid_specstring.split('//')
        leaf_ratiospec=ratio_specstring.split('//')
        assert len(leaf_lipspec)==len(leaf_ratiospec),f'lipid names and mole fractions are not congruent'
        if len(leaf_lipspec)==1:  # symmetrical bilayer
            leaf_lipspec.append(leaf_lipspec[0])
            leaf_ratiospec.append(leaf_ratiospec[0])
        
        global_lipid_names=[]
        # build leaflets
        LL=slices['LOWER-LEAFLET']
        UL=slices['UPPER-LEAFLET']
        LL['lipids']=[]
        UL['lipids']=[]
        for li,(leaflet_name,lipid_name_string,lipid_molfrac_string) in enumerate(zip(['LOWER-LEAFLET','UPPER-LEAFLET'],leaf_lipspec,leaf_ratiospec)):
            lipid_names=lipid_name_string.split(':')
            lipid_molfracs=list(map(float,lipid_molfrac_string.split(':')))
            assert len(lipid_names)==len(lipid_molfracs),f'lipid names and mole fractions are not congruent'
            for ll in lipid_names:
                if not ll in global_lipid_names: global_lipid_names.append(ll)
            slices[leaflet_name]['lipids']=[dict(name=name,frac=frac) for name,frac in zip(lipid_names,lipid_molfracs)]

        lipid_data={}
        addl_params=[]
        for l in lipid_names:
            lpath=os.path.join(self.lipid_pdb_path,l)
            assert os.path.exists(lpath),f'No PDB available for lipid {l}'
            for i in range(10):
                lpdb=os.path.join(lpath,f'{l}-0{i}.pdb')
                # lpdb=os.path.join(lpath,f'{l}-noh-0{i}.pdb')
                shutil.copy(lpdb,'./')
            psc=PSFContents(os.path.join(lpath,f'{l}-init.psf'))
            orients=pd.read_csv(os.path.join(lpath,'orientations.dat'),header=0,index_col=None)
            with open(os.path.join(lpath,'parameters.txt'),'r') as f:
                paramlist=f.read().split('\n')
            ldata={'charge':psc.get_charge(),'parameters':paramlist,'odf':orients}
            for p in paramlist:
                if p.endswith('str') and not p in addl_params:
                    addl_params.append(p)
            lipid_data[l]=ldata

        tot_lip_mols=mem_area/SAPL
        bilayer_charge=0.0
        for leafname in ['LOWER-LEAFLET','UPPER-LEAFLET']:
            for component in slices[leafname]['lipids']:
                component['n']=int(np.floor(component['frac']*tot_lip_mols))
                bilayer_charge+=lipid_data[component['name']]['charge']*component['n']
        global_charge=pro_charge+bilayer_charge
        logger.debug(f'Slices: {slices}')
        sg='+' if global_charge>0 else ''
        logger.debug(f'Total charge: {sg}{global_charge:.3f}')



        self.next_basename('packmol')
        # first packmol run embeds protein and builds only the bilayer        
        pm=PackmolInputWriter(self.config)
        pm.newscript(f'{self.basename}-1')
        packmol_output_pdb=f'{self.basename}-1.pdb'
        pm.comment('packmol input automatically generated by pestifer')
        pm.addline(f'output {packmol_output_pdb}')
        pm.addline(f'filetype pdb')
        if seed is not None:
            pm.addline(f'seed {seed}')
        pm.addline(f'pbc {" ".join([f"{_:.3f}" for _ in box_ll])} {" ".join([f"{_:.3f}" for _ in box_ur])}')
        pm.addline(f'tolerance {tolerance}')
        pm.addline(f'nloop {nloop_all}')
        if 'embed' in self.specs:
            pm.addline(f'structure {pdb}')
            pm.addline( 'number 1',indents=1)
            pm.addline(f'fixed {" ".join([f"{x:.3f}" for x in shift])} 0. 0. 0.',indents=1)
            pm.addline( 'end structure')

        for leaflet in [LL,UL]:
            if leaflet is LL: selp='-05.pdb'
            else: selp='-08.pdb'
            for specs in leaflet['lipids']:
                name=specs['name']
                n=specs['n']
                pm.addline(f'structure {name}{selp}')
                pm.addline(f'number {n}',indents=1)
                odf=lipid_data[name]['odf']
                # can adjust to pick a random sample maybe
                osrs=odf[odf['pdb']==f'{name}{selp}'].iloc[0,:]
                length=osrs.top_z-osrs.bottom_z
                zlims=[leaflet['z-lo'],leaflet['z-hi']]
                avail_depth=zlims[1]-zlims[0]
                if length>avail_depth:
                    specs['below-z']=zlims[0]
                    specs['above-z']=zlims[1]
                    margin=length-avail_depth
                    specs['within-z-lo']=zlims[0]-0.5*margin
                    specs['within-z-hi']=zlims[1]+0.5*margin
                    specs['below-z-atom']=osrs.top_serial if leaflet is LL else osrs.bottom_serial
                    specs['above-z-atom']=osrs.bottom_serial if leaflet is LL else osrs.top_serial
                    pm.addline(f'inside box {box_ll[0]:.3f} {box_ll[1]:.3f} {specs["within-z-lo"]:.3f} {box_ur[0]:.3f} {box_ur[1]:.3f} {specs["within-z-hi"]:.3f}',indents=1)
                    pm.addline(f'atoms {specs["below-z-atom"]}',indents=1)
                    pm.addline(f'below plane 0. 0. 1. {specs["below-z"]:.3f}',indents=2)
                    pm.addline( 'end atoms',indents=1)
                    pm.addline(f'atoms {specs["above-z-atom"]}',indents=1)
                    pm.addline(f'above plane 0. 0. 1. {specs["above-z"]:.3f}',indents=2)
                    pm.addline( 'end atoms',indents=1)
                else:
                    constrain_rotation=180.0 if leaflet is LL else 0.0
                    specs['within-z-lo']=zlims[0]
                    specs['within-z-hi']=zlims[1]
                    pm.addline(f'inside box {box_ll[0]:.3f} {box_ll[1]:.3f} {specs["within-z-lo"]:.3f} {box_ur[0]:.3f} {box_ur[1]:.3f} {specs["within-z-hi"]:.3f}',indents=1)
                    pm.addline(f'constrain_rotation x {constrain_rotation} 10.0',indents=1)
                    pm.addline(f'constrain_rotation y {constrain_rotation} 10.0',indents=1)
                pm.addline(f'nloop {nloop}',indents=1)
                pm.addline(f'end structure')
        pm.writefile()
        self.result=pm.runscript()

        if self.result!=0:
            return super().do()
        anion_qtot=cation_qtot=0
        if sg=='+':
            anion_qtot=int(global_charge)
        else:
            cation_qtot=int(np.abs(global_charge))

        anion_n=int(anion_qtot//np.abs(int(anion_q)))
        cation_n=int(cation_qtot//np.abs(int(cation_q)))
        ions={anion_name:{'n':anion_n,'claimed':0},cation_name:{'n':cation_n,'claimed':0}}
        logger.debug(f'ions {ions}')

        LC=slices['LOWER-CHAMBER']
        LC['AVAILABLE-VOLUME']=mem_area*(LC['z-hi']-LC['z-lo'])
        UC=slices['UPPER-CHAMBER']
        UC['AVAILABLE-VOLUME']=mem_area*(UC['z-hi']-UC['z-lo'])
        # computes the volumes in the upper and lower chambers that are occupied
        # by lipid and/or protein atoms and adjusts the apparent chamber volumes
        # accordingly
        cdf=coorddf_from_pdb(packmol_output_pdb)
        # 1. find the "edges" where density is 1/2 max density in each chamber
        h,z=np.histogram(cdf['z'],bins=100)
        hmax=np.max(h)
        zmean=np.mean(z)
        low_zedge=z[:-1][(h<0.5*hmax)*(z[:-1]<zmean)][-1]
        hi_zedge=z[:-1][(h<0.5*hmax)*(z[:-1]>zmean)][0]
        logger.debug(f'low zedge {low_zedge} hi_zedge {hi_zedge}')
        # 2. Get the coordinates for all atoms in the lower chamber
        cl=cdf[cdf['z']<LC['z-hi']][['x','y','z']].to_numpy()
        if len(cl)>0:
            # compute the volume excluded as the volume of the complex hull over all
            # coordinates plus a 2-angstrom deep exluded layer along the interface
            # this will most likely be an overestimate of the excluded volume, since
            # the protein is likely convex and there are likely convexities 
            # along the membrane
            logger.debug(f'{len(cl)} lipid/protein atoms in the lower chamber')
            cl_ch=ConvexHull(cl)
            exvol=cl_ch.volume*prot_vol_factor
            # exarea=cl_ch.area
            # exdepth=LC['z-hi']-low_zedge
            # exarea-=mem_area+2*(boxdim[0]+boxdim[1])*exdepth
            # assert exarea>mem_area
            # exvol+=exarea*2.0
            LC['AVAILABLE-VOLUME']-=exvol
        cu=cdf[cdf['z']>UC['z-lo']][['x','y','z']].to_numpy()
        if len(cu)>0:
            logger.debug(f'{len(cu)} lipid/protein protein atoms in the upper chamber')
            cu_ch=ConvexHull(cl)
            exvol=cu_ch.volume*prot_vol_factor
            # exarea=cu_ch.area
            # exdepth=hi_zedge-UC['z-lo']
            # exarea-=mem_area+2*(boxdim[0]+boxdim[1])*exdepth
            # assert exarea>mem_area
            # exvol+=exarea*2.0
            UC['AVAILABLE-VOLUME']-=exvol

        total_available_volume=LC['AVAILABLE-VOLUME']+UC['AVAILABLE-VOLUME']
        LC['VOLUME-FRAC']=LC['AVAILABLE-VOLUME']/total_available_volume
        UC['VOLUME-FRAC']=UC['AVAILABLE-VOLUME']/total_available_volume
        for chamber in [LC,UC]:
            chamber['comp']={}
            for sn,sspec in solvent.items():
                mass=sspec['mass']
                x=sspec['mol-frac']
                chamber['comp'][sn]=int(x*molec_n(mass,solution_gcc,chamber['AVAILABLE-VOLUME']))
            for ion,ispec in ions.items():
                chamber['comp'][ion]=int(chamber['VOLUME-FRAC']*ispec['n'])
                ions[ion]['claimed']+=chamber['comp'][ion]

        for ion in ions:
            if ions[ion]['n']>ions[ion]['claimed']:
                ions[ion]['makeup']=ions[ion]['n']-ions[ion]['claimed']
                nn=ions[ion]['makeup']//2
                for chamber in [LC,UC]:
                    chamber['comp'][ion]+=nn
                if ions[ion]['makeup']%2==1:
                    LC['comp'][ion]+=1

        logger.debug(f'slices {slices}')

        # second packmol run builds solvent chambers
        pm.newscript(f'{self.basename}')
        packmol_output_pdb=f'{self.basename}.pdb'
        pm.comment('packmol input automatically generated by pestifer')
        pm.addline(f'output {packmol_output_pdb}')
        pm.addline(f'filetype pdb')
        if self.specs.get('seed',None) is not None:
            pm.addline(f'seed {self.specs.get("seed")}')
        pm.addline(f'pbc {" ".join([f"{_:.3f}" for _ in box_ll])} {" ".join([f"{_:.3f}" for _ in box_ur])}')
        pm.addline(f'tolerance {tolerance}')
        pm.addline(f'nloop {nloop_all}')
        if 'embed' in self.specs:
            pm.addline(f'structure {self.basename}-1.pdb')
            pm.addline( 'number 1',indents=1)
            pm.addline( 'fixed 0. 0. 0. 0. 0. 0.',indents=1)
            pm.addline( 'end structure')

        for chamber in [LC,UC]:
            components=chamber['comp']
            for comp,n in components.items():
                if n>0:
                    pm.addline(f'structure {comp}.pdb')
                    pm.addline(f'number {n}',indents=1)
                    pm.addline(f'inside box {box_ll[0]:.3f} {box_ll[1]:.3f} {chamber["z-lo"]:.3f} {box_ur[0]:.3f} {box_ur[1]:.3f} {chamber["z-hi"]:.3f}',indents=1)
                    pm.addline(f'nloop {nloop}',indents=1)
                    pm.addline(f'end structure')
        pm.writefile()
        self.result=pm.runscript()
        if self.result!=0:
            return super().do()
        # process output pdb to get new psf and pdb
        self.save_state(exts=['pdb'])
        self.next_basename('psfgen')
        self.result=self.psfgen(psf=psf,pdb=pdb,addpdb=packmol_output_pdb,additional_topologies=addl_params)
        if self.result!=0:
            return super().do()
        self.save_state(exts=['psf','pdb'])
        self.log_message('complete')
        return super().do()

    def psfgen(self,psf,pdb,addpdb,additional_topologies=[]):
        pg=self.writers['psfgen']
        pg.newscript(self.basename,additional_topologies=additional_topologies)
        pg.usescript('memb')
        pg.writescript(self.basename,guesscoord=False,regenerate=True,force_exit=True)
        result=pg.runscript(psf=psf,pdb=pdb,addpdb=addpdb,o=self.basename)
        return result
    
    def membrane_embed(self,pdb,zgroups,zval,outpdb='embed.pdb'):
        logger.debug(f'zgroups {zgroups}')
        ztopg,zbotg,zrefg=zgroups
        vm=self.writers['vmd']
        vm.newscript(self.basename,packages=['Orient','PestiferUtil'])
        vm.addline(f'mol new {pdb}')
        vm.addline(f'set TMP [atomselect top all]')
        vm.addline(f'set HEAD [atomselect top "{ztopg}"]')
        vm.addline(f'set TAIL [atomselect top "{zbotg}"]')
        vm.addline(f'set REF  [atomselect top "{zrefg}"]')
        vm.addline(f'vmdcon -info "[measure center $HEAD]"')
        vm.addline(f'vmdcon -info "[measure center $TAIL]"')
        vm.addline(f'set VEC [vecsub [measure center $HEAD] [measure center $TAIL]]')
        vm.addline( '$TMP move [orient $TMP $VEC {0 0 1}]')
        vm.addline(f'set COM [measure center $TMP]')
        vm.addline(f'$TMP moveby [vecscale $COM -1]')
        vm.addline(f'set REFZ [lindex [measure center $REF] 2]')
        vm.addline(f'$TMP moveby [list 0 0 [expr {zval}-($REFZ)]]')
        vm.addline(f'$TMP writepdb {outpdb}')
        vm.addline(f'set minmax [measure minmax $TMP]')
        vm.addline(f'set SEL_R [expr max([MaxRad $HEAD],[MaxRad $TAIL],[MaxRad $REF])]')
        vm.addline(f'set COM [measure center $TMP]')
        vm.addline(f'set f [open "{self.basename}-results.yaml" "w"]')
        vm.addline(r'puts $f "maxrad: [format %.5f $SEL_R]"')
        vm.addline(r'puts $f "lower_left: \[ [join [lindex $minmax 0] ,\ ] \]"')
        vm.addline(r'puts $f "upper_right: \[ [join [lindex $minmax 1] ,\ ] \]"')
        vm.addline(r'puts $f "com: \[ [join $COM ,\ ] \]"')
        vm.addline(r'puts $f "head_z: [format %.5f [lindex [measure center $HEAD] 2]]"')
        vm.addline(r'puts $f "tail_z: [format %.5f [lindex [measure center $TAIL] 2]]"')
        vm.addline(f'close $f')
        vm.writescript()
        vm.runscript()
