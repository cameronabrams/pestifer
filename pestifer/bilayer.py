# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import os 

import numpy as np

from .config import Config
from .controller import Controller
from .util.units import _UNITS_, _SYMBOLS_, cuA_of_nmolec
from .util.util import cell_from_xsc

sA_ =_SYMBOLS_['ANGSTROM']
sA2_=_UNITS_['SQUARE-ANGSTROMS']
sA3_=_UNITS_['CUBIC-ANGSTROMS']

logger=logging.getLogger(__name__)

class BilayerSpecString:
    """ A class for handling bilayer specification strings in memgen format """
    def __init__(self,specstring='',fracstring='',leaflet_delimiter='//',species_delimiter=':'):
        self.specstring=specstring
        self.fracstring=fracstring
        self.leaflet_delimiter=leaflet_delimiter
        self.species_delimiter=species_delimiter
        self.extrastrings={}
        if specstring:
            Sleft,Sright=(specstring.split(leaflet_delimiter)+[specstring])[:2]
            Lleft,Lright=Sleft.split(species_delimiter),Sright.split(species_delimiter)
            nSpleft,nSright=len(Lleft),len(Lright)
            if fracstring:
                FSleft,FSright=(fracstring.split(leaflet_delimiter)+[fracstring])[:2]
                Fleft,Fright=np.array([float(x) for x in FSleft.split(species_delimiter)]),np.array([float(x) for x in FSright.split(species_delimiter)])
                Fleft/=np.sum(Fleft)
                Fright/=np.sum(Fright)
                assert len(Fleft)==nSpleft,f'Left layer has {nSpleft} species but {len(Fleft)} fractions specified'
                assert len(Fright)==nSright,f'Right layer has {nSright} species but {len(Fright)} fractions specified'
            else:
                Fleft=np.array([1.0/nSpleft]*nSpleft)
                Fright=np.array([1.0/nSright]*nSright)
            
            self.left=[dict(name=n,frac=x) for n,x in zip(Lleft,Fleft)]
            self.right=[dict(name=n,frac=x) for n,x in zip(Lright,Fright)]

    def add_specstring(self,attr_name,specstring='',attr_type=str):
        if specstring:
            self.extrastrings[attr_name]=specstring
            Sleft,Sright=(specstring.split(self.leaflet_delimiter)+[specstring])[:2]
            Lleft,Lright=[attr_type(x) for x in Sleft.split(self.species_delimiter)],[attr_type(x) for x in Sright.split(self.species_delimiter)]
            nSpleft,nSright=len(Lleft),len(Lright)
            
            if nSpleft==len(self.left):
                for i in range(nSpleft):
                    self.left[i][attr_name]=Lleft[i]
            elif nSpleft==1:
                for i in range(len(self.left)):
                    self.left[i][attr_name]=Lleft[0]
            if nSright==len(self.right):
                for i in range(nSright):
                    self.right[i][attr_name]=Lright[i]
            elif nSright==1:
                for i in range(len(self.right)):
                    self.right[i][attr_name]=Lright[0]

def specstrings_builddict(lipid_specstring='',lipid_ratio_specstring='',lipid_conformers_specstring='',    
                          solvent_specstring='TIP3',solvent_ratio_specstring=''):
    L=BilayerSpecString(specstring=lipid_specstring,fracstring=lipid_ratio_specstring)
    L.add_specstring('conf',lipid_conformers_specstring,int)
    C=BilayerSpecString(specstring=solvent_specstring,fracstring=solvent_ratio_specstring)
    return {
        'upper_leaflet':L.left,
        'lower_leaflet':L.right,
        'upper_chamber':C.left,
        'lower_chamber':C.right
    }

class Bilayer:
    def __init__(self,composition_dict={},leaflet_nlipids=100,solvent_to_key_lipid_ratio=32.0,
                neutralizing_salt=['POT','CLA'],pdb_collection=None,resi_database=None,solvent_specstring='TIP3',solvent_ratio_specstring='1.0'):

        self.statevars={}

        # leaflet_nlipids is the number of lipids per leaflet in a patch
        self.leaflet_nlipids=leaflet_nlipids

        if not composition_dict:
            logger.debug('Empty bilayer')
            return None

        # complete leaflet entries in composition dictionary with species counts, charges, and MWs
        for l in ['upper_leaflet','lower_leaflet']:
            L=composition_dict[l]
            for d in L:
                if not 'patn' in d:
                    d['patn']=int(d['frac']*leaflet_nlipids)
                if not 'charge' in d:
                    d['charge']=resi_database['lipid'][d['name']].charge
                if not 'MW' in d:
                    d['MW']=resi_database['lipid'][d['name']].mass()

        # user need not have specified the solvent composition in the upper and lower chambers
        if 'upper_chamber' not in composition_dict or 'lower_chamber' not in composition_dict:
            logger.debug(f'Provided composition dictionary does not include solvent chamber specifications')
            logger.debug(f'Using specstrings \'{solvent_specstring}\' and \'{solvent_ratio_specstring}\' to build solvent composition')
            C=BilayerSpecString(specstring=solvent_specstring,fracstring=solvent_ratio_specstring)
            if 'upper_chamber' not in composition_dict:
                composition_dict['upper_chamber']=C.left
            if 'lower_chamber' not in composition_dict:
                composition_dict['lower_chamber']=C.right

        # complete chamber entries in composition dictionary with species counts, charges, and MWs
        for c in ['upper_chamber','lower_chamber']:
            L=composition_dict[c]
            for d in L:
                if not 'patn' in d:
                    d['patn']=int(d['frac']*leaflet_nlipids*solvent_to_key_lipid_ratio)
                if not 'charge' in d:
                    d['charge']=resi_database['water_ions'][d['name']].charge
                if not 'MW' in d:
                    d['MW']=resi_database['water_ions'][d['name']].mass()

        # set up some short-cut object labes
        self.slices={'lower_chamber':{},'lower_leaflet':{},'upper_leaflet':{},'upper_chamber':{}}
        self.LC=self.slices['lower_chamber']
        self.LL=self.slices['lower_leaflet']
        self.UL=self.slices['upper_leaflet']
        self.UC=self.slices['upper_chamber']
        for layer,data in self.slices.items():
            data['composition']=composition_dict[layer]
        # if the bilayer is asymmetric (each leaflet has a unique composition), we cannot assume a priori that
        # each leaflet in a patch has the same number of lipids.  We set a flag to indicate that the patch is 
        # asymmetric and that the number of lipids in each leaflet may need to be adjusted after equilibration
        # and measurment of the pressure profile.
        ul_lx,ll_lx=[(x['name'],x['frac']) for x in self.slices['upper_leaflet']['composition']],[(x['name'],x['frac']) for x in self.slices['lower_leaflet']['composition']]
        self.asymmetric=set(ul_lx)!=set(ll_lx)

        lipid_names=[x['name'] for x in self.slices['upper_leaflet']['composition']]
        lipid_names+= [x['name'] for x in self.slices['lower_leaflet']['composition']]
        self.lipid_names=list(set(lipid_names))
        
        solvent_names=[x['name'] for x in self.slices['upper_chamber']['composition']]
        solvent_names+= [x['name'] for x in self.slices['lower_chamber']['composition']]
        self.solvent_names=list(set(solvent_names))
        self.species_names=self.lipid_names+self.solvent_names

        self.species_data={}
        self.addl_streamfiles=[]
        if pdb_collection is not None:
            for l in self.species_names:
                logger.debug(f'Getting pdb for {l}')
                pdbstruct=pdb_collection.get_pdb(l)
                self.species_data[l]=pdbstruct
                for p in self.species_data[l].get_parameters():
                    if p.endswith('.str') and not p in self.addl_streamfiles:
                        self.addl_streamfiles.append(p)
        logger.debug(f'Additional stream files: {self.addl_streamfiles}')

        self.total_charge=0.0
        for layer,data in self.slices.items():
            data['charge']=0.0
            data['maxlength']=0.0
            data['composition']=composition_dict[layer]
            data['avgMW']=0.0
            data['patn']=0
            for species in data['composition']:
                data['patn']+=species['patn']
                data['charge']+=species['charge']*species['patn']
                if 'chamber' in layer:
                    data['avgMW']+=species['MW']*species['frac']
                elif 'leaflet' in layer:
                    for lipid in data['composition']:
                        lipid['reference_length']=self.species_data[lipid['name']].get_ref_length(index=lipid['conf'])
                        if lipid['reference_length']>data['maxlength']:
                            data['maxlength']=lipid['reference_length']
                    logger.debug(f'{layer} maxlength {data["maxlength"]:.3f}')
            self.total_charge+=data['charge']

        if self.total_charge!=0.0:
            logger.debug(f'Total charge of bilayer is {self.total_charge:.3f} e')
            cation_name,anion_name=neutralizing_salt
            if self.total_charge>0.0: # need to include anion
                ion_name=anion_name
            else:
                ion_name=cation_name
            if ion_name not in self.species_names:
                self.species_names.append(ion_name)
                self.species_data[ion_name]=pdb_collection.get_pdb(ion_name)
            ion_q=resi_database['water_ions'][ion_name].charge
            ion_patn=int(np.round(np.abs(self.total_charge),0)/np.abs(ion_q))
            if ion_patn%2==0:
                lc_ion_patn=ion_patn//2
                uc_ion_patn=ion_patn//2
            else:
                lc_ion_patn=ion_patn//2
                uc_ion_patn=ion_patn//2+1
            for chamber,nions in zip(['upper_chamber','lower_chamber'],[uc_ion_patn,lc_ion_patn]):
                chamber_names=[x['name'] for x in self.slices[chamber]['composition']]
                if ion_name not in chamber_names:
                    self.slices[chamber]['composition'].append({'name':ion_name,'patn':nions,'charge':ion_q,'MW':resi_database['water_ions'][ion_name].mass()})
                else:
                    # if the ion is already in the chamber, just add to the number of ions
                    for species in self.slices[chamber]['composition']:
                        if species['name']==ion_name:
                            species['patn']+=nions
                            break

        # finally, check out all required PDB input files for packmol
        for layer,data in self.slices.items():
            for species in data['composition']:
                species_name=species['name']
                species['local_name']=self.species_data[species_name].checkout()
                # logger.debug(f'Checked out {species_name} as {species["local_name"]}')

    def build_patch(self,SAPL=60.0,xy_aspect_ratio=1.0,midplane_z=0.0,half_mid_zgap=1.0,solution_gcc=1.0,rotation_pm=10.0):
        patch_area=SAPL*self.leaflet_nlipids
        Lx=np.sqrt(patch_area/xy_aspect_ratio)
        Ly=xy_aspect_ratio*Lx
        self.patch_area=patch_area
        lcvol=cuA_of_nmolec(self.LC['avgMW'],solution_gcc,self.LC['patn'])
        ucvol=cuA_of_nmolec(self.UC['avgMW'],solution_gcc,self.UC['patn'])
        lcdepth=lcvol/self.patch_area
        ucdepth=ucvol/self.patch_area
        ll_maxlength=self.LL['maxlength']
        ul_maxlength=self.UL['maxlength']
        ll_actlength=np.cos(np.deg2rad(rotation_pm))*ll_maxlength
        ul_actlength=np.cos(np.deg2rad(rotation_pm))*ul_maxlength
        # guarantees that the longest lipid is longer than the width of the leaflet
        zmin=midplane_z-ll_actlength-lcdepth-half_mid_zgap
        zmax=midplane_z+ul_actlength+ucdepth+half_mid_zgap
        self.UC['z-hi']=zmax
        self.UC['z-lo']=zmax-ucdepth
        self.UL['z-hi']=self.UC['z-lo']
        self.UL['z-lo']=midplane_z+half_mid_zgap
        self.LL['z-hi']=midplane_z-half_mid_zgap
        self.LL['z-lo']=zmin+lcdepth
        self.LC['z-hi']=self.LL['z-lo']
        self.LC['z-lo']=zmin
        logger.debug(f'++++++++++++++ {zmax:8.3f} ++++++++++++++')
        logger.debug(f'                                    ^  ')
        logger.debug(f'    U P P E R   C H A M B E R {zmax-self.UC["z-lo"]:8.3f}')
        logger.debug(f'                                    v  ')
        logger.debug(f'-------------- {self.UC["z-lo"]:8.3f} --------------')
        logger.debug(f'                                    ^  ')
        logger.debug(f'    U P P E R   L E A F L E T {self.UC["z-lo"]-self.UL["z-lo"]:8.3f}    ')
        logger.debug(f'                                    v  ')
        logger.debug(f'-------------- {self.UL["z-lo"]:8.3f} --------------')
        logger.debug(f'    M I D P L A N E  Z - G A P      ')
        logger.debug(f'-------------- {self.LL["z-hi"]:8.3f} --------------')
        logger.debug(f'                                    ^  ')
        logger.debug(f'    L O W E R   L E A F L E T {self.LL["z-hi"]-self.LL["z-lo"]:8.3f}    ')
        logger.debug(f'                                    v  ')
        logger.debug(f'-------------- {self.LC["z-hi"]:8.3f} --------------')
        logger.debug(f'                                    ^  ')
        logger.debug(f'    L O W E R   C H A M B E R {self.LC["z-hi"]-zmin:8.3f}')
        logger.debug(f'                                    v  ')
        logger.debug(f'++++++++++++++ {zmin:8.3f} ++++++++++++++')
        self.patch_ll_corner=np.array([0,0,zmin])
        self.patch_ur_corner=np.array([Lx,Ly,zmax])
        # box and origin
        self.box=np.array([[Lx,0,0],[0,Ly,0],[0,0,zmax-zmin]])
        self.origin=np.array([0,0,zmin])

    def write_packmol(self,pm,half_mid_zgap=2.0,rotation_pm=0.0,nloop=100):
        # first patch-specific packmol directives
        pm.addline(f'pbc {" ".join([f"{_:.3f}" for _ in self.patch_ll_corner])} {" ".join([f"{_:.3f}" for _ in self.patch_ur_corner])}')

        for leaflet in [self.LL,self.UL]:
            logger.debug(f'Leaflet species to pack :{leaflet["composition"]}')
            for specs in leaflet['composition']:
                name=specs['name']
                ref_atoms=self.species_data[name].get_ref_atoms()
                hs=' '.join([f"{x['serial']}" for x in ref_atoms['heads']])
                ts=' '.join([f"{x['serial']}" for x in ref_atoms['tails']])
                n=specs['patn']
                pm.addline(f'structure {specs["local_name"]}')
                pm.addline(f'number {n}',indents=1)
                lipid_length=self.species_data[name].get_ref_length(index=specs['conf'])

                leaflet_thickness=leaflet['z-hi']-leaflet['z-lo']

                if lipid_length>leaflet_thickness:
                    pm.addline(f'inside box {self.patch_ll_corner[0]:.3f} {self.patch_ll_corner[1]:.3f} {self.patch_ll_corner[2]:.3f} {self.patch_ur_corner[0]:.3f} {self.patch_ur_corner[1]:.3f} {self.patch_ur_corner[2]:.3f}',indents=1)
                    if leaflet is self.LL:
                        below_plane_z=leaflet['z-lo']
                        above_plane_z=leaflet['z-hi']-half_mid_zgap
                        pm.addline(f'atoms {hs}',indents=1)
                        pm.addline(f'below plane 0. 0. 1. {below_plane_z:.3f}',indents=2)
                        pm.addline( 'end atoms',indents=1)
                        pm.addline(f'atoms {ts}',indents=1)
                        pm.addline(f'above plane 0. 0. 1. {above_plane_z:.3f}',indents=2)
                        pm.addline(f'below plane 0. 0. 1. 0.0',indents=2)
                        pm.addline( 'end atoms',indents=1)
                    elif leaflet is self.UL:
                        below_plane_z=leaflet['z-lo']+half_mid_zgap
                        above_plane_z=leaflet['z-hi']   
                        pm.addline(f'atoms {ts}',indents=1)
                        pm.addline(f'below plane 0. 0. 1. {below_plane_z:.3f}',indents=2)
                        pm.addline(f'above plane 0. 0. 1. 0.0',indents=2)
                        pm.addline( 'end atoms',indents=1)
                        pm.addline(f'atoms {hs}',indents=1)
                        pm.addline(f'above plane 0. 0. 1. {above_plane_z:.3f}',indents=2)
                        pm.addline( 'end atoms',indents=1)
                else:
                    if leaflet is self.LL:
                        constrain_rotation=180.0
                    elif leaflet is self.UL:
                        constrain_rotation=0.0
                    inside_z_lo=leaflet['z-lo']
                    inside_z_hi=leaflet['z-hi']
                    pm.addline(f'inside box {self.patch_ll_corner[0]:.3f} {self.patch_ll_corner[1]:.3f} {inside_z_lo:.3f} {self.patch_ur_corner[0]:.3f} {self.patch_ur_corner[1]:.3f} {inside_z_hi:.3f}',indents=1)
                    pm.addline(f'constrain_rotation x {constrain_rotation} {rotation_pm}',indents=1)
                    pm.addline(f'constrain_rotation y {constrain_rotation} {rotation_pm}',indents=1)
                pm.addline(f'nloop {nloop}',indents=1)
                pm.addline(f'end structure')
        for chamber in [self.LC,self.UC]:
            logger.debug(f'Chamber species to pack :{chamber["composition"]}')
            for specs in chamber['composition']:
                name=specs['name']
                n=specs['patn']
                pm.addline(f'structure {specs["local_name"]}')
                pm.addline(f'number {n}',indents=1)
                inside_z_lo=chamber['z-lo']
                inside_z_hi=chamber['z-hi']
                pm.addline(f'inside box {self.patch_ll_corner[0]:.3f} {self.patch_ll_corner[1]:.3f} {inside_z_lo:.3f} {self.patch_ur_corner[0]:.3f} {self.patch_ur_corner[1]:.3f} {inside_z_hi:.3f}',indents=1)

                pm.addline(f'nloop {nloop}',indents=1)
                pm.addline(f'end structure')
        pm.writefile()

    def pack_patch(self,pm,specname,seed=None,tolerance=None,nloop_all=200,nloop=200,half_mid_zgap=1.0,rotation_pm=20):
        pm.newscript(specname)
        packmol_output_pdb=f'{specname}.pdb'
        pm.comment('packmol input automatically generated by pestifer')
        pm.addline(f'output {packmol_output_pdb}')
        pm.addline(f'filetype pdb')
        if seed is not None:
            pm.addline(f'seed {seed}')
        pm.addline(f'tolerance {tolerance}')
        pm.addline(f'nloop {nloop_all}')
        self.write_packmol(pm,half_mid_zgap=half_mid_zgap,rotation_pm=rotation_pm,nloop=nloop)
        pm.writefile()
        result=pm.runscript()
        logger.debug(f'{specname} packmol result {result}')
        if result!=0:
            raise Exception(f'Packmol failed with result {result}')
        return packmol_output_pdb

    def equilibrate(self,user_dict={},basename='equilibrate',index=0,spec=''):
        if user_dict=={}:
            return
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        xsc=self.statevars['xsc']
        logger.debug(f'Bilayer area before equilibration: {self.area:.3f} {sA2_}')
        user_dict['tasks']=[
            {'restart':dict(psf=psf,pdb=pdb,xsc=xsc,index=index)},
            {'md':dict(ensemble='minimize',minimize=1000,addl_paramfiles=self.addl_streamfiles)},
            {'md':dict(ensemble='NVT',nsteps=1000,addl_paramfiles=self.addl_streamfiles)},
            {'md':dict(ensemble='NPT',nsteps=200,addl_paramfiles=self.addl_streamfiles,
                       other_parameters=dict(useflexiblecell=True,useconstantratio=True))},
            {'md':dict(ensemble='NPT',nsteps=400,addl_paramfiles=self.addl_streamfiles,
                       other_parameters=dict(useflexiblecell=True,useconstantratio=True))},
            {'md':dict(ensemble='NPT',nsteps=800,addl_paramfiles=self.addl_streamfiles,
                       other_parameters=dict(useflexiblecell=True,useconstantratio=True))},
            {'md':dict(ensemble='NPT',nsteps=1600,addl_paramfiles=self.addl_streamfiles,
                       other_parameters=dict(useflexiblecell=True,useconstantratio=True))},
            {'md':dict(ensemble='NPT',nsteps=3200,addl_paramfiles=self.addl_streamfiles,
                       other_parameters=dict(useflexiblecell=True,useconstantratio=True))},
            {'md':dict(ensemble='NPT',nsteps=6400,addl_paramfiles=self.addl_streamfiles,
                       other_parameters=dict(useflexiblecell=True,useconstantratio=True))},
            {'md':dict(ensemble='NPT',nsteps=12800,addl_paramfiles=self.addl_streamfiles,
                       other_parameters=dict(useflexiblecell=True,useconstantratio=True))},
            {'md':dict(ensemble='NPT',nsteps=25600,addl_paramfiles=self.addl_streamfiles,
                       other_parameters=dict(useflexiblecell=True,useconstantratio=True))},
            {'mdplot':dict(traces=['density',['a_x','b_y','c_z']],legend=True,grid=True,savedata=f'{basename}-traces.csv',basename=basename)},
            {'terminate':dict(basename=basename,chainmapfile=f'{basename}-chainmap.yaml',statefile=f'{basename}-state.yaml')}                 
        ]
        subconfig=Config(userdict=user_dict)
        subcontroller=Controller(subconfig)
        for task in subcontroller.tasks:
            task_key=task.taskname
            task.override_taskname(f'bilayer{spec}-'+task_key)
        
        subcontroller.do_tasks()
        self.statevars=subcontroller.tasks[-1].statevars.copy()
        self.box,self.origin=cell_from_xsc(self.statevars['xsc'])
        self.area=self.box[0][0]*self.box[1][1]
        logger.debug(f'Bilayer area after equilibration: {self.area:.3f} {sA2_}')         

    def delete_lipid(self,vm,count,leaflet='upper'):
        selstr="and z>0.0" if leaflet=='upper' else "and z<0.0"
        pdb=self.statevars['pdb']
        psf=self.statevars['psf']
        basenamepdb,dum=os.path.splitext(pdb)
        basenamepsf,dum=os.path.splitext(psf)
        vm.newscript(f'delete-lipid')
        vm.addline(f'package require PestiferEnviron')
        vm.addline(f'namespace import ::PestiferEnviron::*')
        vm.addline(f'package require psfgen')
        vm.addline(f'readpsf {psf} pdb {pdb}')
        vm.addline(f'mol new {psf}')
        vm.addline(f'mol addfile {pdb}')
        vm.addline(f'set molid [molinfo top get id]')
        vm.addline(f'set residues [leaflet_apportionment $molid]')

        vm.addline(f'set sel [atomselect top "lipid and residue $residues($leaflet)"]')
        vm.addline(f'set lres [lsort -unique [$sel get residue]]')
        vm.addline(f'set nres [llength $lres]')
        vm.addline(r'set shuffled_resnums [lsort -integer -command {expr {rand() > 0.5 ? 1 : -1}} $lres]')
        vm.addline(f'set count {count}')
        vm.addline(r'set remove_resnums [lrange $shuffled_resnums 0 [expr {$count - 1}]]')
        vm.addline(r'set remove_sel [atomselect top "residue $remove_resnums"]')
        vm.addline(r'foreach segname [$remove_sel get segname] resid [$remove_sel get resid] {')
        vm.addline(r'   delatom $segname $resid')
        vm.addline(r'}')
        trim_pdb=basenamepdb+'-trim.pdb'
        trim_psf=basenamepsf+'-trim.psf'
        vm.addline(f'writepdb {trim_pdb}')
        vm.addline(f'writepsf {trim_psf}')
        vm.addline(f'close')
        vm.writescript()
        vm.runscript()
        self.statevars['pdb']=trim_pdb
        self.statevars['psf']=trim_psf