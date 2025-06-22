# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging

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

def specstrings_builddict(lipid_specstring='',lipid_ratio_specstring='',lipid_conformers_specstring='0',    
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
    def __init__(self,composition_dict={},leaflet_nlipids=dict(upper=100,lower=100),solvent_to_key_lipid_ratio=32.0,
                neutralizing_salt=['POT','CLA'],salt_concentration=0.0,solution_gcc=1.0,pdbrepository=None,resi_database=None,solvent_specstring='TIP3',solvent_ratio_specstring='1.0'):

        self.statevars={}

        # leaflet_nlipids is the number of lipids per leaflet in a patch
        self.leaflet_nlipids=leaflet_nlipids

        if not composition_dict:
            logger.debug('Empty bilayer')
            return None

        # complete leaflet entries in composition dictionary with species counts, charges, and MWs
        for l in ['upper_leaflet','lower_leaflet']:
            adjective='upper' if l=='upper_leaflet' else 'lower'
            L=composition_dict[l]
            logger.debug(f'Leaflet {l} composition: {L}')
            for d in L:
                resi=resi_database.get_resi(d['name'])
                if not 'patn' in d:
                    d['patn']=int(d['frac']*leaflet_nlipids[adjective])
                if not 'charge' in d:
                    d['charge']=resi.charge
                if not 'MW' in d:
                    d['MW']=resi.mass()

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
            adjective='upper' if c=='upper_chamber' else 'lower'
            L=composition_dict[c]
            Nsol=0
            AMW=0.0
            for d in L:
                resi=resi_database.get_resi(d['name'])
                if not 'patn' in d:
                    d['patn']=int(d['frac']*leaflet_nlipids[adjective]*solvent_to_key_lipid_ratio)
                if not 'charge' in d:
                    d['charge']=resi.charge
                if not 'MW' in d:
                    d['MW']=resi.mass()
                Nsol+=d['patn']
                AMW+=d['MW']*d['patn']
            if Nsol>0:
                AMW/=Nsol
            else:
                AMW=0.0
            if salt_concentration>0.0:
                Npm=int(AMW/1000*salt_concentration/solution_gcc*Nsol)
                if Npm>0:
                    cation_name,anion_name=neutralizing_salt
                    logger.debug(f'Salting at {salt_concentration} M, soln density {solution_gcc} gcc')
                    logger.debug(f'-> adding {Npm} {cation_name} and {Npm} {anion_name} to {c}')
                    cation,anion=resi_database.get_resi(cation_name),resi_database.get_resi(anion_name)
                    n_cation=int(np.round(Npm/np.abs(cation.charge),0))
                    n_anion=int(np.round(Npm/np.abs(anion.charge),0))
                    composition_dict[c].append({'name':cation_name,'patn':n_cation,'charge':cation.charge,'MW':cation.mass()})
                    composition_dict[c].append({'name':anion_name,'patn':n_anion,'charge':anion.charge,'MW':anion.mass()})
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
        if pdbrepository is not None:
            for l in self.species_names:
                logger.debug(f'Getting pdb for {l}')
                if not l in pdbrepository:
                    raise Exception(f'Cannot find {l} in PDB repository')
                pdbstruct=pdbrepository.checkout(l)
                self.species_data[l]=pdbstruct
                for p in self.species_data[l].get_parameters():
                    if p.endswith('.str') and not p in self.addl_streamfiles:
                        self.addl_streamfiles.append(p)
        logger.debug(f'Additional stream files: {self.addl_streamfiles}')

        self.total_charge=0.0
        for layer,data in self.slices.items():
            data['charge']=0.0
            data['maxthickness']=0.0
            data['composition']=composition_dict[layer]
            data['avgMW']=0.0
            data['patn']=0
            for species in data['composition']:
                logger.debug(f'Layer {layer} species {species["name"]} patn {species["patn"]} charge {species["charge"]} MW {species["MW"]}')
                data['patn']+=species['patn']
                data['charge']+=species['charge']*species['patn']
                if 'chamber' in layer:
                    data['avgMW']+=species['MW']*species['patn']
                elif 'leaflet' in layer:
                    for lipid in data['composition']:
                        lipid['reference_length']=self.species_data[lipid['name']].get_head_tail_length(conformerID=lipid.get('conf',0))
                        if lipid['reference_length']>data['maxthickness']:
                            data['maxthickness']=lipid['reference_length']
                    logger.debug(f'{layer} maxthickness {data["maxthickness"]:.3f}')
            logger.debug(f'Layer {layer} patn {data["patn"]} charge {data["charge"]:.3f} e avgMW*N {data["avgMW"]:.3f} g/mol')
            self.total_charge+=data['charge']
            data['avgMW']/=data['patn']

        if self.total_charge!=0.0:
            logger.debug(f'Total charge of bilayer is {self.total_charge:.3f} e')
            cation_name,anion_name=neutralizing_salt
            if self.total_charge>0.0: # need to include anion
                ion_name=anion_name
            else:
                ion_name=cation_name
            if ion_name not in self.species_names:
                self.species_names.append(ion_name)
                self.species_data[ion_name]=pdbrepository.checkout(ion_name)
            ion_resi=resi_database.get_resi(ion_name)
            ion_q=ion_resi.charge
            logger.debug(f'Adding {ion_name} with charge {ion_q:.3f} e')
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
                    self.slices[chamber]['composition'].append({'name':ion_name,'patn':nions,'charge':ion_q,'MW':ion_resi.mass()})
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
                conformerID=species.get('conf',0)
                noh=species.get('noh',False)
                species['local_name']=self.species_data[species_name].get_pdb(conformerID=conformerID,noh=noh)
                # logger.debug(f'Checked out {species_name} as {species["local_name"]}')

    def build_patch(self,SAPL=75.0,xy_aspect_ratio=1.0,half_mid_zgap=1.0,solution_gcc=1.0,rotation_pm=10.0):
        patch_area=SAPL*self.leaflet_nlipids['upper'] # assume symmetric
        Lx=np.sqrt(patch_area/xy_aspect_ratio)
        Ly=xy_aspect_ratio*Lx
        self.patch_area=patch_area
        lc_vol=cuA_of_nmolec(self.LC['avgMW'],solution_gcc,self.LC['patn'])
        uc_vol=cuA_of_nmolec(self.UC['avgMW'],solution_gcc,self.UC['patn'])
        lc_thickness=lc_vol/self.patch_area
        uc_thickness=uc_vol/self.patch_area
        ll_maxthickness=self.LL['maxthickness']
        ul_maxthickness=self.UL['maxthickness']
        ll_actthickness=np.cos(np.deg2rad(rotation_pm))*ll_maxthickness
        ul_actthickness=np.cos(np.deg2rad(rotation_pm))*ul_maxthickness
        # guarantees that the longest lipid is longer than the width of the leaflet
        zmin=0.0
        zmax=lc_thickness+ll_actthickness+2*half_mid_zgap+ul_actthickness+uc_thickness
        self.LC['z-lo']=zmin
        self.LC['z-hi']=self.LC['z-lo']+lc_thickness
        self.LL['z-lo']=self.LC['z-hi']
        self.LL['z-hi']=self.LL['z-lo']+ll_actthickness
        self.UL['z-lo']=self.LL['z-hi']+2*half_mid_zgap
        self.midplane_z=self.LL['z-hi']+half_mid_zgap
        self.UL['z-hi']=self.UL['z-lo']+ul_actthickness
        self.UC['z-lo']=self.UL['z-hi']
        self.UC['z-hi']=zmax
        logger.debug(f'++++++++++++++ {zmax:8.3f} ++++++++++++++')
        logger.debug(f'                                    ^  ')
        logger.debug(f'    U P P E R   C H A M B E R {zmax-self.UC["z-lo"]:8.3f}')
        logger.debug(f'                                    v  ')
        logger.debug(f'-------------- {self.UC["z-lo"]:8.3f} --------------')
        logger.debug(f'                                    ^  ')
        logger.debug(f'    U P P E R   L E A F L E T {self.UC["z-lo"]-self.UL["z-lo"]:8.3f}    ')
        logger.debug(f'                                    v  ')
        logger.debug(f'-------------- {self.UL["z-lo"]:8.3f} --------------')
        logger.debug(f'    M I D P L A N E  H A L F  G A P    ')
        logger.debug(f'-------------- {self.midplane_z:8.3f} --------------')
        logger.debug(f'    M I D P L A N E  H A L F  G A P    ')
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
        self.origin=np.array([self.box[i][i]/2 for i in range(3)])

    def write_packmol(self,pm,half_mid_zgap=2.0,rotation_pm=0.0,nloop=100):
        # first patch-specific packmol directives
        pm.addline(f'pbc {" ".join([f"{_:.3f}" for _ in self.patch_ll_corner])} {" ".join([f"{_:.3f}" for _ in self.patch_ur_corner])}')
        ll=self.patch_ll_corner
        ur=self.patch_ur_corner
        for leaflet in [self.LL,self.UL]:
            logger.debug(f'Leaflet species to pack :{leaflet["composition"]}')
            for specs in leaflet['composition']:
                name=specs['name']
                logger.debug(f'Packing {name}')
                lipid_max_length=self.species_data[name].get_max_internal_length(conformerID=specs.get('conf',0))
                lipid_headtail_length=self.species_data[name].get_head_tail_length(conformerID=specs.get('conf',0))
                lipid_overhang=lipid_max_length-lipid_headtail_length
                leaflet_thickness=leaflet['z-hi']-leaflet['z-lo']
                ref_atoms=self.species_data[name].get_ref_atoms()
                hs=' '.join([f"{x['serial']}" for x in ref_atoms['heads']])
                ts=' '.join([f"{x['serial']}" for x in ref_atoms['tails']])
                n=specs['patn']
                pm.addline(f'structure {specs["local_name"]}')
                pm.comment(f'  max int length {lipid_max_length:.3f}, head-tail length {lipid_headtail_length:.3f}, overhang {lipid_overhang:.3f}')
                pm.addline(f'number {n}',indents=1)
                # if the maximum length of the lipid is less than the desired leaflet thickness minus a margin,
                # we can pack directly into the leaflet slab using contstrained rotation to orient
                if lipid_max_length<leaflet_thickness:
                    logger.debug(f' -> lipid length {lipid_max_length:.3f} < leaflet thickness {leaflet_thickness:.3f}')
                    if leaflet is self.LL:
                        constrain_rotation=180.0
                    elif leaflet is self.UL:
                        constrain_rotation=0.0
                    inside_z_lo=leaflet['z-lo']
                    inside_z_hi=leaflet['z-hi']
                    pm.addline(f'inside box {ll[0]:.3f} {ll[1]:.3f} {inside_z_lo:.3f} {ur[0]:.3f} {ur[1]:.3f} {inside_z_hi:.3f}',indents=1)
                    pm.addline(f'constrain_rotation x {constrain_rotation} {rotation_pm}',indents=1)
                    pm.addline(f'constrain_rotation y {constrain_rotation} {rotation_pm}',indents=1)
                else:
                    # we need to pack by specifying some atoms above a plane and other atoms below a different
                    # plane, and the "inside box" should refer to the whole cell
                    # if the head-tail length is GREATER than the leaflet thickness, we can use the
                    # explicit leafleat boundaries as the planes
                    if lipid_headtail_length>leaflet_thickness:
                        below_plane_z=leaflet['z-lo']
                        above_plane_z=leaflet['z-hi']-half_mid_zgap
                    else:
                        # if the head-tail length is less than the leaflet thickness, while the max internal
                        # length is still greater, then we can't use the leaflet boundaries as the planes
                        span=leaflet_thickness-lipid_headtail_length
                        below_plane_z=leaflet['z-lo']+span/2.0
                        above_plane_z=leaflet['z-hi']-span/2.0
                    pm.addline(f'inside box {ll[0]:.3f} {ll[1]:.3f} {ll[2]:.3f} {ur[0]:.3f} {ur[1]:.3f} {ur[2]:.3f}',indents=1)
                    pm.addline(f'atoms {hs}',indents=1)
                    if leaflet is self.LL: # heads are low
                        pm.addline(f'below plane 0. 0. 1. {below_plane_z:.3f}',indents=2)
                    elif leaflet is self.UL: # heads are high
                        pm.addline(f'above plane 0. 0. 1. {above_plane_z:.3f}',indents=2)
                    pm.addline( 'end atoms',indents=1)
                    pm.addline(f'atoms {ts}',indents=1)
                    if leaflet is self.LL: # tails are high
                        pm.addline(f'above plane 0. 0. 1. {above_plane_z:.3f}',indents=2)
                    elif leaflet is self.UL: # tails are low
                        pm.addline(f'below plane 0. 0. 1. {below_plane_z:.3f}',indents=2)
                    pm.addline( 'end atoms',indents=1)

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
                pm.addline(f'inside box {ll[0]:.3f} {ll[1]:.3f} {inside_z_lo:.3f} {ur[0]:.3f} {ur[1]:.3f} {inside_z_hi:.3f}',indents=1)
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

    def equilibrate(self,user_dict={},
                    basename='equilibrate',index=0,
                    relaxation_protocol=None,parent_controller_index=0):
        if user_dict=={}:
            return
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        xsc=self.statevars['xsc']
        logger.debug(f'Bilayer area before equilibration: {self.area:.3f} {sA2_}')
        if not relaxation_protocol:
            logger.debug(f'Using hard-coded relaxation protocol for {basename}!!')
            relaxation_protocol=[            
                {'md':dict(ensemble='minimize',minimize=1000)},
                {'md':dict(ensemble='NVT',nsteps=1000)},
                {'md':dict(ensemble='NPT',nsteps=200)},
                {'md':dict(ensemble='NPT',nsteps=400)},
                {'md':dict(ensemble='NPT',nsteps=800)},
                {'md':dict(ensemble='NPT',nsteps=1600)},
                {'md':dict(ensemble='NPT',nsteps=3200)},
                {'md':dict(ensemble='NPT',nsteps=6400)},
                {'md':dict(ensemble='NPT',nsteps=12800)},
                {'md':dict(ensemble='NPT',nsteps=25600)}]
        else:
            logger.debug(f'Using user-specified relaxation protocol: {relaxation_protocol}')
        for stage in relaxation_protocol:
            specs=stage['md']
            specs['addl_paramfiles']=self.addl_streamfiles
            if specs.get('ensemble',None) in ['NPT','npt']:
                if not 'other_parameters' in specs: # never true due to ycleptic base.yaml
                    specs['other_parameters']={'useflexiblecell':True,'useconstantratio':True,
                                               'pressureProfile':'on','pressureProfileSlabs':30,
                                               'pressureProfileFreq':100}
                else:
                    if not 'useflexiblecell' in specs['other_parameters']:
                        specs['other_parameters']['useflexiblecell']=True
                    if not 'useconstantratio' in specs['other_parameters']:
                        specs['other_parameters']['useconstantratio']=True
                    if user_dict['namd']['processor-type']!='gpu': # GPU NAMD 3.0.1 does not support pressure profiles
                        if not 'pressureProfile' in specs['other_parameters']:
                            specs['other_parameters']['pressureProfile']='on'
                        if not 'pressureProfileSlabs' in specs['other_parameters']:
                            specs['other_parameters']['pressureProfileSlabs']=30
                        if not 'pressureProfileFreq' in specs['other_parameters']:
                            specs['other_parameters']['pressureProfileFreq']=100
        traces=['density',['a_x','b_y','c_z']]
        profiles=['pressure']
        if user_dict['namd']['processor-type']!='gpu':
            traces.append('pressure') # To do: change this to pressureProfile plotting
        user_dict['tasks']=[
            {'restart':dict(psf=psf,pdb=pdb,xsc=xsc,index=index)}
            ]+relaxation_protocol+[
            {'mdplot':dict(traces=traces,profiles=profiles,legend=True,grid=True,basename=basename)},
            {'terminate':dict(basename=basename,chainmapfile=f'{basename}-chainmap.yaml',statefile=f'{basename}-state.yaml')}                 
        ]
        user_dict['title']=f'Bilayer equilibration from {basename}'
        subconfig=Config(userdict=user_dict,quiet=True)
        subcontroller=Controller(subconfig,index=parent_controller_index+1)
        for task in subcontroller.tasks:
            task_key=task.taskname
            task.override_taskname(f'{basename}-'+task_key)
        
        subcontroller.do_tasks()
        self.statevars=subcontroller.tasks[-1].statevars.copy()
        self.box,self.origin=cell_from_xsc(self.statevars['xsc'])
        self.area=self.box[0][0]*self.box[1][1]
        logger.debug(f'{basename} area after equilibration: {self.area:.3f} {sA2_}')         
