# Author: Cameron F. Abrams, <cfa2@drexel.edu>
import logging
import numpy as np
import yaml

from ..basetask import BaseTask
from ..charmmtop import CharmmResiDatabase
from ..config import Config
from ..util.coord import coorddf_from_pdb
from ..packmol import PackmolInputWriter
from ..psfutil.psfcontents import PSFContents
from ..util.units import _UNITS_, _SYMBOLS_
from ..stringthings import my_logger
from ..util.util import cell_to_xsc
from ..util.units import cuA_of_nmolec, nmolec_in_cuA

sA_ =_SYMBOLS_['ANGSTROM']
sA2_=_UNITS_['SQUARE-ANGSTROMS']
sA3_=_UNITS_['CUBIC-ANGSTROMS']

logger=logging.getLogger(__name__)

def bilayer_stringsplit(input_string,delimiter0='//',delimiter1=':',return_type=str):
    """Splits a string into two lists of strings, separated by the first delimiter.
    The two lists are separated by the second delimiter.
    The first delimiter is used to split the string into two parts, and the second delimiter
    is used to split each part into a list of strings.
    The function returns two lists of strings.
    If the input string is empty, the function returns two empty lists. """
    uls,lls=(input_string.split(delimiter0)+[input_string])[:2]
    ul,ll=uls.split(delimiter1),lls.split(delimiter1)
    if return_type==str:
        return ul,ll
    ul,ll=list(map(return_type,ul)),list(map(return_type,ll))
    return np.array(ul),np.array(ll)

class Bilayer:
    def __init__(self,composition_dict={},lipid_specstring='',lipid_ratio_specstring='',lipid_conformers_specstring='',    
                leaflet_patch_nlipids=100,solvent_specstring='',solvent_ratio_specstring='',solvent_to_lipid_ratio_specstring='',
                neutralizing_salt=['POT','CLA'],pdb_collection=None,resi_database=None):
        # leaflet_patch_nlipids is the number of lipids per leaflet in a patch
        self.leaflet_patch_nlipids=leaflet_patch_nlipids
        self.slices={'lower_chamber':{},'lower_leaflet':{},'upper_leaflet':{},'upper_chamber':{}}
        self.LC=self.slices['lower_chamber']
        self.LL=self.slices['lower_leaflet']
        self.UL=self.slices['upper_leaflet']
        self.UC=self.slices['upper_chamber']
        if not composition_dict and not lipid_specstring and not lipid_ratio_specstring and not lipid_conformers_specstring:
            logger.debug('Empty bilayer')
            return
        if not composition_dict:
            # old-style bilayer composition specification with memgen-style specstrings
            ul_lip,ll_lip=bilayer_stringsplit(lipid_specstring)
            ul_x,ll_x=bilayer_stringsplit(lipid_ratio_specstring,return_type=float)
            ul_xs,ll_xs=np.sum(ul_x),np.sum(ll_x)
            ul_x/=ul_xs
            ll_x/=ll_xs
            ul_c,ll_c=bilayer_stringsplit(lipid_conformers_specstring,return_type=int)
            uc_s,lc_s=bilayer_stringsplit(solvent_specstring)
            uc_x,lc_x=bilayer_stringsplit(solvent_ratio_specstring,return_type=float)
            uc_xs,lc_xs=np.sum(uc_x),np.sum(lc_x)
            uc_x/=uc_xs
            lc_x/=lc_xs
            uc_slr,lc_slr=bilayer_stringsplit(solvent_to_lipid_ratio_specstring,return_type=float)
            assert len(ul_lip)==len(ul_x),f'Upper leaflet has {len(ul_lip)} lipids but {len(ul_x)} mole fractions specified'
            assert len(ul_lip)==len(ul_c),f'Upper leaflet has {len(ul_lip)} lipids but {len(ul_c)} conformers specified'
            assert len(ul_x)==len(ul_c),f'Upper leaflet has {len(ul_x)} mole fractions but {len(ul_c)} conformers specified'
            assert len(ll_lip)==len(ll_x),f'Upper leaflet has {len(ul_lip)} lipids but {len(ul_x)} mole fractions specified'
            assert len(ll_lip)==len(ll_c),f'Upper leaflet has {len(ul_lip)} lipids but {len(ul_c)} conformers specified'
            assert len(ll_x)==len(ll_c),f'Upper leaflet has {len(ul_x)} mole fractions but {len(ul_c)} conformers specified'
            # build the equivalent new-style composition dictionary; 'patn' is the number of molecules for a minimal patch
            composition_dict={
                'upper_chamber':[{'name':n,'frac':x,'patn':int(x*leaflet_patch_nlipids*s),'MW':resi_database['water_ions'][n].mass()} for n,x,s in zip(uc_s,uc_x,uc_slr)],
                'lower_chamber':[{'name':n,'frac':x,'patn':int(x*leaflet_patch_nlipids*s),'MW':resi_database['water_ions'][n].mass()} for n,x,s in zip(lc_s,lc_x,lc_slr)],
                'upper_leaflet':[{'name':n,'frac':x,'conf':c,'patn':int(x*leaflet_patch_nlipids)} for n,x,c in zip(ul_lip,ul_x,ul_c)],
                'lower_leaflet':[{'name':n,'frac':x,'conf':c,'patn':int(x*leaflet_patch_nlipids)} for n,x,c in zip(ll_lip,ll_x,ll_c)],
            }
        else:
            ul_lip,ll_lip=[x['name'] for x in composition_dict['upper_leaflet']],[x['name'] for x in composition_dict['lower_leaflet']]
            ul_x,ll_x=[x['frac'] for x in composition_dict['upper_leaflet']],[x['frac'] for x in composition_dict['lower_leaflet']]
            # new-style bilayer composition specification with composition dictionary
            for l in ['upper_leaflet','lower_leaflet']:
                L=composition_dict[l]
                for d in L:
                    d['patn']=int(d['frac']*leaflet_patch_nlipids)
            # user need not have specified the solvent composition in the upper and lower chambers
            if 'upper_chamber' not in composition_dict or 'lower_chamber' not in composition_dict:
                uc_s,lc_s=bilayer_stringsplit(solvent_specstring)
                uc_x,lc_x=bilayer_stringsplit(solvent_ratio_specstring,return_type=float)
                uc_xs,lc_xs=np.sum(uc_x),np.sum(lc_x)
                uc_x/=uc_xs
                lc_x/=lc_xs
                if 'upper_chamber' not in composition_dict:
                    composition_dict['upper_chamber']=[{'name':n,'frac':x,'MW':resi_database['water_ions'][n].mass()} for n,x in zip(uc_s,uc_x)]
                if 'lower_chamber' not in composition_dict:
                    composition_dict['lower_chamber']=[{'name':n,'frac':x,'MW':resi_database['water_ions'][n].mass()} for n,x in zip(lc_s,lc_x)]

            uc_slr,lc_slr=bilayer_stringsplit(solvent_to_lipid_ratio_specstring,return_type=float)
            for c,slr in zip(['upper_chamber','lower_chamber'],[uc_slr,lc_slr]):
                L=composition_dict[c]
                for d in L:
                    d['patn']=int(d['frac']*leaflet_patch_nlipids*slr)

        # if the bilayer is asymmetric (each leaflet has a unique composition), we cannot assume a priori that
        # each leaflet in a patch has the same number of lipids.  We set a flag to indicate that the patch is 
        # asymmetric and that the number of lipids in each leaflet may need to be adjusted after equilibration
        # and measurment of the pressure profile.
        ul_lx,ll_lx=[(x['name'],x['frac']) for x in composition_dict['upper_leaflet']],[(x['name'],x['frac']) for x in composition_dict['lower_leaflet']]
        self.asymmetric=set(ul_lx)!=set(ll_lx)
        
        lipid_names=[x['name'] for x in composition_dict['upper_leaflet']]
        lipid_names+= [x['name'] for x in composition_dict['lower_leaflet']]
        self.lipid_names=list(set(lipid_names))
        
        solvent_names=[x['name'] for x in composition_dict['upper_chamber']]
        solvent_names+= [x['name'] for x in composition_dict['lower_chamber']]
        self.solvent_names=list(set(solvent_names))
        
        if pdb_collection is not None:
            self.species_data={}
            self.addl_streamfiles=[]
            for l in self.lipid_names+self.solvent_names:
                pdbstruct=pdb_collection.get_pdb(l)
                self.species_data[l]=pdbstruct
                for p in self.species_data[l].get_parameters():
                    if p.endswith('.str') and not p in self.addl_streamfiles:
                        self.addl_streamfiles.append(p)

        self.total_charge=0.0
        for layer,data in self.slices.items():
            data['charge']=0.0
            data['maxlength']=0.0
            data['composition']=composition_dict[layer]
            data['avgMW']=0.0
            data['patn']=0
            for species in data['composition']:
                data['patn']+=species['patn']
                if 'chamber' in layer:
                    data['charge']+=self.species_data[species['name']].get_charge()
                    data['avgMW']+=species['MW']*species['frac']
                elif 'leaflet' in layer:
                    data['charge']+=self.species_data[species['name']].get_charge()
                    for lipid in data:
                        lipid['reference_length']=self.species_data[lipid['name']].get_ref_length(index=lipid['conf'])
                        if lipid['reference_length']>data['maxlength']:
                            data['maxlength']=lipid['reference_length']
            self.total_charge+=self.layer_charge[layer]

        if self.total_charge!=0.0:
            logger.debug(f'Total charge of bilayer is {self.total_charge:.3f} e')
            cation_name,anion_name=neutralizing_salt
            if self.total_charge>0.0: # need to include anion
                ion_name=anion_name
            else:
                ion_name=cation_name
            self.species_data[ion_name]=self.pdb_collection.get_pdb(ion_name)
            ion_q=self.solvent_data[ion_name].get_charge()
            ion_patn=int(np.round(np.abs(self.total_charge),0)/np.abs(ion_q))
            if ion_patn%2==0:
                lc_ion_patn=ion_patn//2
                uc_ion_patn=ion_patn//2
            else:
                lc_ion_patn=ion_patn//2
                uc_ion_patn=ion_patn//2+1
            self.UC['composition'].append({'name':ion_name,'patn':uc_ion_patn})
            self.LC['composition'].append({'name':ion_name,'patn':lc_ion_patn}) 

        for layer,data in self.slices.items():
            composition=data['composition']
            for species in composition:
                species_name=species['name']

    # def set_slice_bounds(self,chamber_thickness=[5.0,5.0],midplane_z=0.0,leaflet_thickness=[22.0,22.0]):
    #     global_z_min=midplane_z-chamber_thickness[0]-leaflet_thickness[0]
    #     global_z_max=midplane_z+chamber_thickness[1]+leaflet_thickness[1]
    #     self.LC['z-lo']= global_z_min
    #     self.LC['z-hi']= midplane_z-leaflet_thickness[0]
    #     self.LL['z-lo']= midplane_z-leaflet_thickness[0]
    #     self.LL['z-hi']= midplane_z
    #     self.UL['z-lo']= midplane_z
    #     self.UL['z-hi']= midplane_z+leaflet_thickness[1]
    #     self.UC['z-lo']= midplane_z+leaflet_thickness[1]
    #     self.UC['z-hi']= global_z_max
    #     for k,v in self.slices.items():
    #         v['THICKNESS']=v['z-hi']-v['z-lo']

    # def volumizer(self,lateral_area):
    #     for k,v in self.slices.items():
    #         v['INIT-VOLUME']=v['THICKNESS']*lateral_area
    #         v['INIT-NWATEREQUIV']=nmolec_in_cuA(18.0,1.0,v['INIT-VOLUME'])

    def build_patch(self,SAPL=60.0,xy_aspect_ratio=1.0,midplane_z=0.0):
        patch_area=SAPL*self.leaflet_patch_nlipids
        Lx=np.sqrt(patch_area/xy_aspect_ratio)
        Ly=xy_aspect_ratio*Lx
        self.patch_area=patch_area
        lcvol=cuA_of_nmolec(self.LC['avgMW'],1.0,self.LC['patn'])
        ucvol=cuA_of_nmolec(self.UC['avgMW'],1.0,self.UC['patn'])
        lcdepth=lcvol/self.patch_area
        ucdepth=ucvol/self.patch_area
        zmin=midplane_z-self.LL['maxlength']-lcdepth
        zmax=midplane_z+self.UL['maxlength']+ucdepth
        self.patch_ll_corner=np.array([0,0,zmin])
        self.patch_ur_corner=np.array([Lx,Ly,zmax])
        # self.LL_charge=sum([x['patn']*self.lipid_data[x['name']].get_charge() for x in self.LL['lipids']])
        # self.UL_charge=sum([x['patn']*self.lipid_data[x['name']].get_charge() for x in self.UL['lipids']])
        # patch_charge=self.LL_charge+self.UL_charge
        # self.patch_charge=patch_charge
        sg='+' if self.patch_charge>0 else ''
        logger.debug(f'Total patch charge: {sg}{self.patch_charge:.3f}')
        anion_qtot=cation_qtot=0
        if sg=='+':
            anion_qtot=int(self.patch_charge)
        else:
            cation_qtot=int(np.abs(np.round(self.patch_charge)))
        for leaflet in [self.LL,self.UL]:
            for component in leaflet['lipids']:
                component['local_name']=self.lipid_data[component['name']].checkout(index=component['conf'])
        self.anion_n=int(anion_qtot//np.abs(int(anion_q)))
        self.cation_n=int(cation_qtot//np.abs(int(cation_q)))
        ions={anion_name:{'n':anion_n,'claimed':0},cation_name:{'n':cation_n,'claimed':0}}

    def solvate_patch(self,solvent_specstring,solv_molfrac_specstring,cation_name='POT',cation_q=1,anion_name='CLA',anion_q=-1):
        solvent_names=solvent_specstring.split(':')
        solvent_molfracs=np.array(list(map(float,solv_molfrac_specstring.split(':'))))
        assert len(solvent_names)==len(solvent_molfracs),f'You have specified {len(solvent_names)} solvent names but {len(solvent_molfracs)} mole fractions'
        self.solvent={}
        for sn,sm in zip(solvent_names,solvent_molfracs):
            solvent[sn]={}
            solvent[sn]['mol-frac']=sm
            solv_pdbstruct=self.pdb_collection.get_pdb(sn)
            solv_pdbstruct.checkout()
            solvent[sn]['mass']=self.RDB['water_ions'][sn].mass()

        anion_qtot=cation_qtot=0
        if sg=='+':
            anion_qtot=int(self.patch_charge)
        else:
            cation_qtot=int(np.abs(np.round(self.patch_charge)))

        anion_n=int(anion_qtot//np.abs(int(anion_q)))
        cation_n=int(cation_qtot//np.abs(int(cation_q)))
        ions={anion_name:{'n':anion_n,'claimed':0},cation_name:{'n':cation_n,'claimed':0}}


class BilayerEmbedTask(BaseTask):
    """ A class for handling embedding proteins into bilayers
    
    Attributes
    ----------
    yaml_header(str) 

    Methods
    -------
    do(): 
        Based on specs, writes packmol input files to generate a membrane-embedded protein and
        then runs packmol

    """
    yaml_header='bilayer'
    def __init__(self,input_dict,taskname,config:Config,writers,prior):
        super().__init__(input_dict,taskname,config,writers,prior)
        self.progress=config.progress
        self.pdb_collection=config.RM.pdb_collection
        self.RDB=CharmmResiDatabase()
        self.RDB.add_stream('lipid')
        self.RDB.add_topology('toppar_all36_moreions.str',streamnameoverride='water_ions')
        lipid_specstring=self.specs.get('lipids','')
        ratio_specstring=self.specs.get('mole_fractions','')
        conformers_specstring=self.specs.get('conformers','')
        solvent_specstring=self.specs.get('solvents','')
        solvent_ratio_specstring=self.specs.get('solvent_mole_fractions','')
        solvent_to_lipid_ratio_specstring=self.specs.get('solvent_to_lipid_ratio','')
        composition_dict=self.specs.get('composition',{})
        self.patch=Bilayer(composition_dict,
                            lipid_specstring=lipid_specstring,
                            ratio_specstring=ratio_specstring,
                            conformers_specstring=conformers_specstring,
                            solvent_specstring=solvent_specstring,
                            solvent_ratio_specstring=solvent_ratio_specstring,
                            solvent_to_lipid_ratio_specstring=solvent_to_lipid_ratio_specstring,
                            pdb_collection=self.pdb_collection)

    def do(self):
        # TODO: switch from full packmol to the following
        # 1. use packmol to build a minimal patch with the desired composition
        # 2. equilibrate the hell out of that patch
        # 3. use a modified version of the membrane plugin to generate a full bilayer
        # 4. minimize and thermalize that bilayer
        # 5. use psfgen to embed the protein, deleting the bad lipids
        # 6. use packmol to solvate and ionize (more flexible than vmd/solvate)
        # 7. profit!
        self.log_message('initiated')
        self.inherit_state()
        self.pro_psf=self.statevars.get('psf',None)
        if self.pro_psf is not None:
            self.pro_psc=PSFContents(self.pro_psf)
            self.pro_charge=self.pro_psc.get_charge()
            self.pro_pdb=self.statevars.get('pdb',None)
            if self.pro_psf is not None and self.pro_pdb is not None:
                logger.debug(f'BilayerEmbedTask will use psf {self.pro_psf} and pdb {self.pro_pdb} as inputs')

        self.build_patch()
        # self.relax_patch(patch)
        # self.make_bilayer_from_patch(patch)
        # self.embed_protein()
        # self.solvate()
        # self.log_message('complete')
        return super().do()

    def build_patch(self):
        self.next_basename('build_patch')

        rotation_pm=self.specs.get('rotation_pm',20.)
        # fuzz_factor=self.specs.get('fuzz_factor',0.5)
        half_mid_zgap=self.specs.get('half_mid_zgap',1.0)

        SAPL=self.specs.get('SAPL',60.0)
        # scale_excluded_volume=self.specs.get('scale_excluded_volume',1.0)
        seed=self.specs.get('seed',27021972)
        tolerance=self.specs.get('tolerance',2.0)
        # solution_gcc=self.specs.get('solution_gcc',1.0)
        xy_aspect_ratio=self.specs.get('xy_aspect_ratio',1.0)
        # leaflet_thickness=self.specs.get('leaflet_thickness',27.0) # initial value
        # chamber_thickness=self.specs.get('chamber_thickness',10.0)
        nloop=self.specs.get('nloop',200)
        nloop_all=self.specs.get('nloop_all',200)

        self.patch.set_slice_bounds(chamber_thickness=[chamber_thickness,chamber_thickness],midplane_z=0.0,leaflet_thickness=[leaflet_thickness,leaflet_thickness])
        self.patch.build_patch(SAPL=SAPL,xy_aspect_ratio=xy_aspect_ratio)
        pm=PackmolInputWriter(self.config)
        pm.newscript(f'{self.basename}-patch')
        packmol_output_pdb=f'{self.basename}-patch.pdb'
        pm.comment('packmol input automatically generated by pestifer')
        pm.addline(f'output {packmol_output_pdb}')
        pm.addline(f'filetype pdb')
        if seed is not None:
            pm.addline(f'seed {seed}')
        pm.addline(f'pbc {" ".join([f"{_:.3f}" for _ in self.patch.patch_ll_corner])} {" ".join([f"{_:.3f}" for _ in self.patch.patch_ur_corner])}')
        pm.addline(f'tolerance {tolerance}')
        pm.addline(f'nloop {nloop_all}')

        for leaflet in [self.patch.LL,self.patch.UL]:
            for specs in leaflet['lipids']:
                name=specs['name']
                ref_atoms=self.patch.lipid_data[name].get_ref_atoms()
                hs=' '.join([f"{x['serial']}" for x in ref_atoms['heads']])
                ts=' '.join([f"{x['serial']}" for x in ref_atoms['tails']])
                n=specs['patn']
                pm.addline(f'structure {specs["local_name"]}')
                pm.addline(f'number {n}',indents=1)
                lipid_length=self.patch.lipid_data[name].get_ref_length(index=specs['conf'])
                Dz=np.cos(np.deg2rad(rotation_pm))*lipid_length
                fuzz=Dz-lipid_length # negative
                lipid_interior_length=lipid_length+2*fuzz-2.0
                # fuzz_out,fuzz_in=fuzz*fuzz_factor,fuzz*(1-fuzz_factor)
                leaflet_thickness=leaflet['z-hi']-leaflet['z-lo']
                logger.debug(f'{name}: {specs["local_name"]} length {lipid_length} lipid interior length {lipid_interior_length:.2f}')
                if lipid_length>leaflet_thickness:
                    if leaflet is self.patch.LL:
                        below_plane_z=leaflet['z-hi']-lipid_interior_length-half_mid_zgap
                        above_plane_z=leaflet['z-hi']-half_mid_zgap
                        below_plane_atoms=hs
                        above_plane_atoms=ts
                    elif leaflet is self.patch.UL:
                        below_plane_z=leaflet['z-lo']+half_mid_zgap
                        above_plane_z=leaflet['z-lo']+lipid_interior_length+half_mid_zgap   
                        below_plane_atoms=ts
                        above_plane_atoms=hs
                    pm.addline(f'inside box {self.patch.patch_ll_corner[0]:.3f} {self.patch.patch_ll_corner[1]:.3f} {self.patch.patch_ll_corner[2]:.3f} {self.patch.patch_ur_corner[0]:.3f} {self.patch.patch_ur_corner[1]:.3f} {self.patch.patch_ur_corner[2]:.3f}',indents=1)
                    pm.addline(f'atoms {below_plane_atoms}',indents=1)
                    pm.addline(f'below plane 0. 0. 1. {below_plane_z:.3f}',indents=2)
                    pm.addline( 'end atoms',indents=1)
                    pm.addline(f'atoms {above_plane_atoms}',indents=1)
                    pm.addline(f'above plane 0. 0. 1. {above_plane_z:.3f}',indents=2)
                    pm.addline( 'end atoms',indents=1)
                else:
                    if leaflet is self.patch.LL:
                        constrain_rotation=180.0
                    elif leaflet is self.patch.UL:
                        constrain_rotation=0.0
                    inside_z_lo=leaflet['z-lo']
                    inside_z_hi=leaflet['z-hi']
                    pm.addline(f'inside box {self.patch.patch_ll_corner[0]:.3f} {self.patch.patch_ll_corner[1]:.3f} {inside_z_lo:.3f} {self.patch.patch_ur_corner[0]:.3f} {self.patch.patch_ur_corner[1]:.3f} {inside_z_hi:.3f}',indents=1)
                    pm.addline(f'constrain_rotation x {constrain_rotation} {rotation_pm}',indents=1)
                    pm.addline(f'constrain_rotation y {constrain_rotation} {rotation_pm}',indents=1)
                pm.addline(f'nloop {nloop}',indents=1)
                pm.addline(f'end structure')
        pm.writefile()
        self.result=pm.runscript()
        logger.debug(f'patch packmol result {self.result}')
        if self.result!=0:
            return super().do()
    
    # def solvate(self):
    #     solvent_specstring=self.specs.get('solvents','TIP3')
    #     solv_molfrac_specstring=self.specs.get('solvent_mole_fractions','1.0')
    #     cation_name=self.specs.get('cation','POT')
    #     cation_q=self.RDB['water_ions'][cation_name].charge
    #     anion_name=self.specs.get('anion','CLA')
    #     anion_q=self.RDB['water_ions'][anion_name].charge
    #     assert cation_name in self.RDB['water_ions'],f'Cation {cation_name} not found in available CHARMM topologies'
    #     assert anion_name in self.RDB['water_ions'],f'Anion {anion_name} not found in available CHARMM topologies'
    #     cation_pdbstruct=self.pdb_collection.get_pdb(cation_name)
    #     cation_pdbstruct.checkout()
    #     anion_pdbstruct=self.pdb_collection.get_pdb(anion_name)
    #     anion_pdbstruct.checkout()
    #     salt_con_M=self.specs.get('salt_con',0.0)
    #     solvent_names=solvent_specstring.split(':')
    #     solvent_molfracs=np.array(list(map(float,solv_molfrac_specstring.split(':'))))
    #     assert len(solvent_names)==len(solvent_molfracs),f'You have specified {len(solvent_names)} solvent names but {len(solvent_molfracs)} mole fractions'
    #     solvent={}
    #     for sn,sm in zip(solvent_names,solvent_molfracs):
    #         assert sn in self.RDB['water_ions'],f'solvent {sn} is not found in the available CHARMM topologies'
    #         solvent[sn]={}
    #         solvent[sn]['mol-frac']=sm
    #         solv_pdbstruct=self.pdb_collection.get_pdb(sn)
    #         solv_pdbstruct.checkout()
    #         solvent[sn]['mass']=self.RDB['water_ions'][sn].mass()

    #     my_logger(solvent,logger.debug)

    #     # boxdim=np.array(list(map(float,self.specs.get('dims',[]))))

    # def embed_protein(self):
    #     if 'embed' in self.specs:
    #         assert pdb is not None
    #         # use the protein and embedding specs to size the box
    #         embed_specs=self.specs['embed']
    #         prot_rad_scal=embed_specs.get('protein_radius_scaling',1.0)
    #         xydist=embed_specs.get('xydist',0.0)
    #         zdist=embed_specs.get('zdist',0.0)
    #         self.next_basename('embed')
    #         self.membrane_embed(pdb,[embed_specs['z_head_group'],embed_specs['z_tail_group'],embed_specs['z_ref_group']['text']],embed_specs['z_ref_group']['z_value'],outpdb=f'{self.basename}.pdb')
    #         self.save_state(exts=['pdb'])
    #         pdb=self.statevars.get('pdb',None)
    #         with open(f'{self.basename}-results.yaml','r') as f:
    #             embed_results=yaml.safe_load(f)
    #         protein_rad=float(embed_results["maxrad"])
    #         protein_ll=np.array(list(map(float,embed_results["lower_left"])))
    #         protein_ur=np.array(list(map(float,embed_results["upper_right"])))
    #         # protein_com=np.array(list(map(float,embed_results["com"])))
    #         midplane=0.0
    #         # logger.debug(f'protein ll {protein_ll} ur {protein_ur} com {protein_com}')
    #         box_ll=protein_ll-np.array([xydist,xydist,zdist])
    #         box_ur=protein_ur+np.array([xydist,xydist,zdist])
            
    #         shift=-box_ll
    #         box_ur+=shift
    #         box_ll+=shift
    #         midplane+=shift[2]

    #         origin=0.5*(box_ur+box_ll)
    #         boxdim=box_ur-box_ll
    #         box_area=boxdim[0]*boxdim[1]
    #         protein_xyarea=np.pi*(protein_rad*prot_rad_scal)**2
    #         logger.debug(f'box area {box_area:.3f} {sA2_}')
    #         logger.debug(f'protein radius {protein_rad:.3f} {sA_} xyarea {protein_xyarea:.3f} {sA2_}')
    #         mem_area=box_area-protein_xyarea
    #         logger.debug(f'protein excludes approx. {np.floor(protein_xyarea/SAPL)} lipids per leaflet')
    #         logger.debug(f'actual membrane area {mem_area:.3f} {sA2_}')
    #         lower_chamber_thickness=midplane-leaflet_thickness-box_ll[2]
    #         upper_chamber_thickness= box_ur[2]-(midplane+leaflet_thickness)
    #         bilayer.set_slice_bounds(chamber_thickness=[lower_chamber_thickness,upper_chamber_thickness],
    #                                  midplane_z=midplane,
    #                                  leaflet_thickness=[leaflet_thickness,leaflet_thickness])
    #     else:
    #         box_ll=np.zeros(3,dtype=float)
    #         box_ur=boxdim
    #         origin=0.5*(box_ur+box_ll)
    #         box_area=boxdim[0]*boxdim[1]
    #         mem_area=box_area
    #         lower_chamber_thickness=midplane-leaflet_thickness-box_ll[2]
    #         upper_chamber_thickness= box_ur[2]-(midplane+leaflet_thickness)
    #         bilayer.set_slice_bounds(chamber_thickness=[lower_chamber_thickness,upper_chamber_thickness],
    #                                  midplane_z=midplane,
    #                                  leaflet_thickness=[leaflet_thickness,leaflet_thickness])

    #     bilayer.volumizer(box_area)
        
    #     my_logger(self.slices,logger.debug)

    #     boxV=boxdim.prod()
    #     logger.debug(f'box corners ll {box_ll} ur {box_ur}')
    #     logger.debug(f'box dimensions {boxdim}')
    #     logger.debug(f'box area {box_area} {sA2_}; membrane_area {mem_area} {sA2_}')
    #     logger.debug(f'box volume {boxV:.3f} {sA3_} (check sum-over-slices: {np.sum(np.array([x["INIT-VOLUME"] for x in bilayer.slices.values()]))} {sA3_}) (or {nmolec_in_cuA(18.0,1.0,boxV)} waters at 1 g/cc)')
    #     logger.debug(f'membrane volume {box_area*2*leaflet_thickness:.3f} {sA3_}')
    #     logger.debug(f'   lower chamber: {bilayer.LC["INIT-VOLUME"]} {sA3_}; upper chamber: {bilayer.UC["INIT-VOLUME"]} {sA3_}')
    #     logger.debug(f'due to membrane, expecting {nmolec_in_cuA(18.0,1.0,boxV-(box_area*2*leaflet_thickness))} waters')
    #     logger.debug(f'   lower chamber: {bilayer.LC["INIT-NWATEREQUIV"]}; upper chamber: {bilayer.UC["INIT-NWATEREQUIV"]}')

    #     cell_to_xsc(np.diag(boxdim),origin,f'{self.basename}.xsc')
    #     self.save_state(exts=['xsc'])


    #     logger.debug(f'leaflets')
    #     my_logger(bilayer.LL,logger.debug)
    #     my_logger(bilayer.UL,logger.debug)



    #     # TODO: build the minimal patch!

    #     self.next_basename('packmol')
    #     # first packmol run builds the bilayer patch
  
        
    #     """ finish building bilayer with embedded protein, no solvent """
        
    #     anion_qtot=cation_qtot=0
    #     if sg=='+':
    #         anion_qtot=int(global_charge)
    #     else:
    #         cation_qtot=int(np.abs(np.round(global_charge)))

    #     anion_n=int(anion_qtot//np.abs(int(anion_q)))
    #     cation_n=int(cation_qtot//np.abs(int(cation_q)))
    #     ions={anion_name:{'n':anion_n,'claimed':0},cation_name:{'n':cation_n,'claimed':0}}
    #     my_logger(ions,logger.debug)

    #     # computes the volumes in the upper and lower chambers that are occupied
    #     # by lipid and/or protein atoms and adjusts the apparent chamber volumes
    #     # accordingly
    #     cdf=coorddf_from_pdb(packmol_output_pdb)
    #     density_calc_resolution=2. # Angstrom
    #     nxbins=int(boxdim[0]/density_calc_resolution)
    #     nybins=int(boxdim[1]/density_calc_resolution)
    #     for S in [LC,UC]:
    #         slice_atoms=cdf[(cdf['z']<S['z-hi'])&(cdf['z']>S['z-lo'])]
    #         my_logger(slice_atoms.head(),logger.debug)
    #         logger.debug(f'{slice_atoms.shape[0]} atoms')
    #         S['AVAILABLE-VOLUME']=S['INIT-VOLUME']
    #         if slice_atoms.shape[0]>0:
    #             nzbins=int((S['THICKNESS'])/density_calc_resolution)
    #             h,edges=np.histogramdd(slice_atoms[['x','y','z']].to_numpy(),bins=(nxbins,nybins,nzbins))
    #             binV=boxdim[0]*boxdim[1]*(S['THICKNESS'])/(nxbins*nybins*nzbins)
    #             logger.debug(f'bin volume {binV:.3f} {sA3_}')
    #             hit_bin_count=h[h>0.0].size
    #             logger.debug(f'{hit_bin_count}/{h.size} bins have some atom density')
    #             excluded_volume=hit_bin_count*binV
    #             logger.debug(f'excluded volume {excluded_volume} {sA3_}')
    #             S['EXCLUDED-VOLUME']=excluded_volume
    #             S['AVAILABLE-VOLUME']-=S['EXCLUDED-VOLUME']*scale_excluded_volume
    #             logger.debug(f'available volume {S["AVAILABLE-VOLUME"]} {sA3_} ({nmolec_in_cuA(18.0,1.0,S["AVAILABLE-VOLUME"])} water-equivs)')

    #     my_logger(self.slices,logger.debug)

    #     total_available_volume=LC['AVAILABLE-VOLUME']+UC['AVAILABLE-VOLUME']
    #     logger.debug(f'total volume available for solvent {total_available_volume:.3f} {sA3_} ({total_available_volume/boxV*100:2f} of box)')
    #     for chamber in [LC,UC]:
    #         chamber['VOLUME-FRAC']=chamber['AVAILABLE-VOLUME']/total_available_volume
    #         chamber['comp']={}
    #         for sn,sspec in solvent.items():
    #             mass=sspec['mass']
    #             x=sspec['mol-frac']
    #             chamber['comp'][sn]=int(x*nmolec_in_cuA(mass,solution_gcc,chamber['AVAILABLE-VOLUME']))
    #             logger.debug(f'solvent component {sn}: mass {mass} x {x} n {chamber["comp"][sn]}')
    #         for ion,ispec in ions.items():
    #             chamber['comp'][ion]=int(chamber['VOLUME-FRAC']*ispec['n'])
    #             ions[ion]['claimed']+=chamber['comp'][ion]

    #     for ion in ions:
    #         if ions[ion]['n']>ions[ion]['claimed']:
    #             ions[ion]['makeup']=ions[ion]['n']-ions[ion]['claimed']
    #             nn=ions[ion]['makeup']//2
    #             for chamber in [LC,UC]:
    #                 chamber['comp'][ion]+=nn
    #             if ions[ion]['makeup']%2==1:
    #                 LC['comp'][ion]+=1

    #     my_logger(self.slices,logger.debug)

    #     # second packmol run builds solvent chambers
    #     pm.newscript(f'{self.basename}')
    #     packmol_output_pdb=f'{self.basename}.pdb'
    #     pm.comment('packmol input automatically generated by pestifer')
    #     pm.addline(f'output {packmol_output_pdb}')
    #     pm.addline(f'filetype pdb')
    #     if self.specs.get('seed',None) is not None:
    #         pm.addline(f'seed {self.specs.get("seed")}')
    #     pm.addline(f'pbc {" ".join([f"{_:.3f}" for _ in box_ll])} {" ".join([f"{_:.3f}" for _ in box_ur])}')
    #     pm.addline(f'tolerance {tolerance}')
    #     pm.addline(f'nloop {nloop_all}')
    #     if 'embed' in self.specs:
    #         pm.addline(f'structure {self.basename}-1.pdb')
    #         pm.addline( 'number 1',indents=1)
    #         pm.addline( 'fixed 0. 0. 0. 0. 0. 0.',indents=1)
    #         pm.addline( 'end structure')

    #     for chamber in [LC,UC]:
    #         components=chamber['comp']
    #         for comp,n in components.items():
    #             if n>0:
    #                 pm.addline(f'structure {comp}.pdb')
    #                 pm.addline(f'number {n}',indents=1)
    #                 pm.addline(f'inside box {box_ll[0]:.3f} {box_ll[1]:.3f} {chamber["z-lo"]:.3f} {box_ur[0]:.3f} {box_ur[1]:.3f} {chamber["z-hi"]:.3f}',indents=1)
    #                 pm.addline(f'nloop {nloop}',indents=1)
    #                 pm.addline(f'end structure')
    #     pm.writefile()
    #     self.result=pm.runscript()
    #     if self.result!=0:
    #         return super().do()
    #     # process output pdb to get new psf and pdb
    #     self.save_state(exts=['pdb'])
    #     self.next_basename('psfgen')
    #     self.result=self.psfgen(psf=psf,pdb=pdb,addpdb=packmol_output_pdb,additional_topologies=addl_streamfiles)
    #     if self.result!=0:
    #         return super().do()
    #     self.save_state(exts=['psf','pdb'])
    #     self.log_message('complete')
    #     return super().do()

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
