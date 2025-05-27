# Author: Cameron F. Abrams, <cfa2@drexel.edu>
import logging
import yaml

import numpy as np

from copy import deepcopy

from ..bilayer import Bilayer, specstrings_builddict
from ..basetask import BaseTask
from ..charmmtop import CharmmResiDatabase
from ..config import Config
from ..scriptwriters import PackmolInputWriter
from ..psfutil.psfcontents import get_toppar_from_psf
from ..util.util import cell_to_xsc,cell_from_xsc, protect_str_arg
from ..util.units import _UNITS_

sA2_=_UNITS_['SQUARE-ANGSTROMS']

logger=logging.getLogger(__name__)

class MakeMembraneSystemTask(BaseTask):
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
    yaml_header='make_membrane_system'
    def __init__(self,config_specs={},controller_specs={}):
        super().__init__(config_specs,controller_specs)
        self.patchA=self.patchB=self.patch=None
        self.progress=self.config.progress
        self.pdb_collection=self.config.RM.pdb_collection
        self.RDB=CharmmResiDatabase()
        self.RDB.add_stream('lipid')
        self.RDB.add_topology('toppar_all36_moreions.str',streamnameoverride='water_ions')
        self.bilayer_specs=self.specs.get('bilayer',{})
        self.embed_specs=self.specs.get('embed',{})
        self.using_prebuilt_bilayer=False
        if 'prebuilt' in self.bilayer_specs and 'pdb' in self.bilayer_specs['prebuilt']:
            logger.debug('Using prebuilt bilayer')
            self.using_prebuilt_bilayer=True
            self.quilt=Bilayer()
            self.quilt.statevars['pdb']=self.bilayer_specs['prebuilt']['pdb']
            self.quilt.statevars['psf']=self.bilayer_specs['prebuilt']['psf']
            self.quilt.statevars['xsc']=self.bilayer_specs['prebuilt']['xsc']
            self.quilt.box,self.quilt.origin=cell_from_xsc(self.quilt.statevars['xsc'])
            self.quilt.area=self.quilt.box[0][0]*self.quilt.box[1][1]
            additional_topologies=get_toppar_from_psf(self.quilt.statevars['psf'])
            self.quilt.addl_streamfiles=additional_topologies
        else:
            self.initialize()

    def initialize(self):
        lipid_specstring=self.bilayer_specs.get('lipids','')
        ratio_specstring=self.bilayer_specs.get('mole_fractions','')
        conformers_specstring=self.bilayer_specs.get('conformers','')
        solvent_specstring=self.bilayer_specs.get('solvents','TIP3')
        solvent_ratio_specstring=self.bilayer_specs.get('solvent_mole_fractions','1.0')
        solvent_to_lipid_ratio=self.bilayer_specs.get('solvent_to_lipid_ratio',32.0)
        leaflet_nlipids=self.bilayer_specs.get('leaflet_nlipids',100)
        cation_name=self.bilayer_specs.get('cation','POT')
        anion_name=self.bilayer_specs.get('anion','CLA')
        neutralizing_salt=[cation_name,anion_name]
        composition_dict=self.bilayer_specs.get('composition',{})

        if not composition_dict:
            composition_dict=specstrings_builddict(lipid_specstring,
                                                  ratio_specstring,
                                                  conformers_specstring,
                                                  solvent_specstring,
                                                  solvent_ratio_specstring)
        logger.debug(f'Main composition dict {composition_dict}')
        self.patch=Bilayer(composition_dict,
                            neutralizing_salt=neutralizing_salt,
                            solvent_specstring=solvent_specstring,
                            solvent_ratio_specstring=solvent_ratio_specstring,
                            solvent_to_key_lipid_ratio=solvent_to_lipid_ratio,
                            leaflet_nlipids=leaflet_nlipids,
                            pdb_collection=self.pdb_collection,resi_database=self.RDB)
        logger.debug(f'Main composition dict after call {composition_dict}')
        if self.patch.asymmetric:
            logger.debug(f'Requested patch is asymmetric; generating two symmetric patches')
            logger.debug(f'Symmetrizing bilayer to upper leaflet')
            composition_dict['lower_leaflet_saved']=composition_dict['lower_leaflet']
            composition_dict['lower_chamber_saved']=composition_dict['lower_chamber']
            composition_dict['lower_leaflet']=composition_dict['upper_leaflet']
            composition_dict['lower_chamber']=composition_dict['upper_chamber']
            self.patchA=Bilayer(composition_dict,
                                neutralizing_salt=neutralizing_salt,
                                solvent_specstring=solvent_specstring,
                                solvent_ratio_specstring=solvent_ratio_specstring,
                                solvent_to_key_lipid_ratio=solvent_to_lipid_ratio,
                                leaflet_nlipids=leaflet_nlipids,
                                pdb_collection=self.pdb_collection,resi_database=self.RDB)
            logger.debug(f'Symmetrizing bilayer to lower leaflet')
            composition_dict['upper_leaflet_saved']=composition_dict['upper_leaflet']
            composition_dict['upper_chamber_saved']=composition_dict['upper_chamber']
            composition_dict['lower_leaflet']=composition_dict['lower_leaflet_saved']
            composition_dict['lower_chamber']=composition_dict['lower_chamber_saved']
            composition_dict['upper_leaflet']=composition_dict['lower_leaflet']
            composition_dict['upper_chamber']=composition_dict['lower_chamber']
            self.patchB=Bilayer(composition_dict,
                                neutralizing_salt=neutralizing_salt,
                                solvent_specstring=solvent_specstring,
                                solvent_ratio_specstring=solvent_ratio_specstring,
                                solvent_to_key_lipid_ratio=solvent_to_lipid_ratio,
                                leaflet_nlipids=leaflet_nlipids,
                                pdb_collection=self.pdb_collection,resi_database=self.RDB)
            composition_dict['upper_leaflet']=composition_dict['upper_leaflet_saved']
            composition_dict['upper_chamber']=composition_dict['upper_chamber_saved']
            self.patch=None


    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        # as part of a list of tasks, this task expects to be fed a protein system to embed
        self.pro_psf=self.statevars.get('psf',None)
        if self.pro_psf is not None:
            self.pro_pdb=self.statevars.get('pdb',None)
            if self.pro_psf is not None and self.pro_pdb is not None:
                logger.debug(f'will use psf {self.pro_psf} and pdb {self.pro_pdb} as inputs')

        if not self.using_prebuilt_bilayer:
            self.build_patch()
            self.make_bilayer_from_patch()
        self.embed_protein()
        self.log_message('complete')
        return super().do()

    def build_patch(self):
        logger.debug(f'Bilayer specs: {self.bilayer_specs}')
        rotation_pm=self.bilayer_specs.get('rotation_pm',10.)
        half_mid_zgap=self.bilayer_specs.get('half_mid_zgap',1.0)
        SAPL=self.bilayer_specs.get('SAPL',60.0)
        seed=self.bilayer_specs.get('seed',27021972)
        tolerance=self.bilayer_specs.get('tolerance',2.0)
        xy_aspect_ratio=self.bilayer_specs.get('xy_aspect_ratio',1.0)
        nloop=self.bilayer_specs.get('nloop',100)
        nloop_all=self.bilayer_specs.get('nloop_all',100)
        relaxation_protocols=self.bilayer_specs.get('relaxation_protocols',{})
        relaxation_protocol=relaxation_protocols.get('patch',{})
        logger.debug(f'relaxation protocols: {relaxation_protocols}')
        # we now build the patch, or if asymmetric, two patches
        for patch,spec in zip([self.patch,self.patchA,self.patchB],['','A','B']):
            if patch is None:
                continue
            self.next_basename(f'patch{spec}')
            specname=self.basename
            logger.debug(f'building {specname}')
            patch.build_patch(SAPL=SAPL,xy_aspect_ratio=xy_aspect_ratio,
                              rotation_pm=rotation_pm,
                              half_mid_zgap=half_mid_zgap)
            pm=PackmolInputWriter(self.config)
            packmol_output_pdb=patch.pack_patch(pm,specname,seed=seed,
                                                tolerance=tolerance,
                                                nloop_all=nloop_all,
                                                half_mid_zgap=half_mid_zgap,
                                                rotation_pm=rotation_pm,
                                                nloop=nloop)
            self.next_basename(f'patch{spec}-psfgen')
            pg=self.writers['psfgen']
            pg.newscript(self.basename,additional_topologies=patch.addl_streamfiles)
            pg.usescript('bilayer_patch')
            pg.writescript(self.basename,guesscoord=False,regenerate=True,force_exit=True)
            result=pg.runscript(pdb=packmol_output_pdb,o=self.basename)
            cell_to_xsc(patch.box,patch.origin,f'{self.basename}.xsc')
            patch.area=patch.box[0][0]*patch.box[1][1]
            patch.statevars['pdb']=f'{self.basename}.pdb'
            patch.statevars['psf']=f'{self.basename}.psf'
            patch.statevars['xsc']=f'{self.basename}.xsc'
            self.next_basename(f'patch{spec}-equilibrate')
            patch.equilibrate(user_dict=deepcopy(self.config['user']),
                              basename=self.basename,index=self.index,
                              spec=spec,relaxation_protocol=relaxation_protocol,
                              parent_controller_index=self.controller_index)

    def make_bilayer_from_patch(self):
        logger.debug(f'Creating bilayer from patch')
        self.next_basename('bilayer')
        additional_topologies=[]
        if self.patch is not None:
            pdb=self.patch.statevars.get('pdb',None)
            xsc=self.patch.statevars.get('xsc',None)
            psf=self.patch.statevars.get('psf',None)
            pdbA=pdbB=pdb
            psfA=psfB=psf
            xscA=xscB=xsc
        elif self.patchA is not None and self.patchB is not None:
            psfA=self.patchA.statevars.get('psf',None)
            pdbA=self.patchA.statevars.get('pdb',None)
            xscA=self.patchA.statevars.get('xsc',None)
            psfB=self.patchB.statevars.get('psf',None)
            pdbB=self.patchB.statevars.get('pdb',None)
            xscB=self.patchB.statevars.get('xsc',None)

        for patch in [self.patch,self.patchA,self.patchB]:
            if patch is None:
                continue
            additional_topologies+=patch.addl_streamfiles
        additional_topologies=list(set(additional_topologies))
        pg=self.writers['psfgen']
        pg.newscript(self.basename,additional_topologies=additional_topologies)
        pg.usescript('bilayer_quilt')
        pg.writescript(self.basename,guesscoord=False,regenerate=False,force_exit=True,writepsf=False,writepdb=False)
        margin=self.embed_specs.get('xydist',10.0)
        if hasattr(self,"pro_pdb"):
            # we will eventually embed a protein in here, so send its pdb along to help size the bilayer
            result=pg.runscript(propdb=self.pro_pdb,margin=margin,psfA=psfA,pdbA=pdbA,psfB=psfB,pdbB=pdbB,xscA=xscA,xscB=xscB,o=self.basename)
        else:
            dimx,dimy=self.bilayer_specs.get('dims',(0,0))
            npatchx,npatchy=self.bilayer_specs.get('npatch',(0,0))
            if npatchx!=0 and npatchy!=0:
                result=pg.runscript(nx=npatchx,ny=npatchy,psfA=psfA,pdbA=pdbA,
                                    psfB=psfB,pdbB=pdbB,xscA=xscA,xscB=xscB,o=self.basename)
            elif dimx!=0 and dimy!=0:
                result=pg.runscript(dimx=dimx,dimy=dimy,psfA=psfA,pdbA=pdbA,
                                    psfB=psfB,pdbB=pdbB,xscA=xscA,xscB=xscB,o=self.basename)
        self.quilt=Bilayer()
        self.quilt.addl_streamfiles=additional_topologies
        self.statevars['pdb']=f'{self.basename}.pdb'
        self.statevars['psf']=f'{self.basename}.psf'
        self.statevars['xsc']=f'{self.basename}.xsc'
        self.statevars['topologies']=additional_topologies
        self.quilt.statevars=self.statevars.copy()
        self.quilt.box,self.quilt.origin=cell_from_xsc(f'{self.basename}.xsc')
        self.quilt.area=self.quilt.box[0][0]*self.quilt.box[1][1]
        relaxation_protocol=self.bilayer_specs.get('relaxation_protocols',{}).get('bilayer',{})
        self.quilt.equilibrate(user_dict=deepcopy(self.config['user']),
                                spec='',
                                relaxation_protocol=relaxation_protocol,
                                parent_controller_index=self.controller_index)

    def embed_protein(self):
        if not self.embed_specs:
            logger.debug('No embed specs.')
            return
        z_head_group=self.embed_specs.get('z_head_group',None)
        z_tail_group=self.embed_specs.get('z_tail_group',None)
        z_ref_group=self.embed_specs.get('z_ref_group',{}).get('text',None)
        z_value=self.embed_specs.get('z_ref_group',{}).get('z_value',0.0)
        self.next_basename('embed')
        pg=self.writers['psfgen']
        pg.newscript(self.basename,additional_topologies=self.quilt.addl_streamfiles)
        pg.usescript('bilayer_embed')
        pg.writescript(self.basename,guesscoord=False,regenerate=True,force_exit=True,writepsf=False,writepdb=False)
        result=pg.runscript(psf=self.pro_psf,
                            pdb=self.pro_pdb,
                            bilayer_psf=self.quilt.statevars['psf'],
                            bilayer_pdb=self.quilt.statevars['pdb'],
                            bilayer_xsc=self.quilt.statevars['xsc'],
                            z_head_group=protect_str_arg(z_head_group),
                            z_tail_group=protect_str_arg(z_tail_group),
                            z_ref_group=protect_str_arg(z_ref_group),
                            z_value=z_value,
                            o=self.basename)
        self.statevars['pdb']=f'{self.basename}.pdb'
        self.statevars['psf']=f'{self.basename}.psf'
        self.statevars['coor']=f'{self.basename}.coor'
        self.statevars['xsc']=f'{self.basename}.xsc'
        del self.statevars['vel']
        if 'charmmff_paramfiles' not in self.statevars:
            self.statevars['charmmff_paramfiles']=[]
        self.statevars['charmmff_paramfiles']+=self.quilt.addl_streamfiles
        self.statevars['charmmff_paramfiles']=list(set(self.statevars['charmmff_paramfiles']))
        self.quilt.statevars.update(self.statevars)
        return result
    
    def fill_box(self):
        fill_top=False
        fill_bot=False
        zdist=self.embed_specs.get('zdist',10.0)
        q_tolerance=self.embed_specs.get('q_tolerance',1.e-4)
        embed_results_file=f'{self.basename}_embed_prefilling.yaml'
        with open(embed_results_file,'r') as f:
            self.embed_results=yaml.safe_load(f)
        net_charge=self.embed_results['net_charge']
        botz=self.embed_results['protein']['min_z']-zdist
        if botz<self.embed_results['bilayer']['min_z']:
            fillwidth_bot=self.embed_results['bilayer']['min_z']-botz
            if fillwidth_bot<3.0:
                fillwidth_bot=0.0
            else:
                bot_zupper=self.embed_results['bilayer']['min_z']
                bot_zlower=botz
                fill_bot=True
        topz=self.embed_results['protein']['max_z']+zdist
        if topz>self.embed_results['bilayer']['max_z']:
            fillwidth_top=topz-self.embed_results['bilayer']['max_z']
            if fillwidth_top<3.0:
                fillwidth_top=0.0
            else:
                top_zupper=self.embed_results['bilayer']['max_z']
                top_zlower=topz
                fill_top=True
        if not fill_top and not fill_bot:
            if np.abs(net_charge)>q_tolerance:
                logger.debug(f'Net charge {net_charge} exceeds tolerance {q_tolerance}')
                # do_autoionize()
            else:
                logger.debug('No filling/ionization needed')
            return
        
        
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

    def psfgen(self,psf='',pdb='',addpdb='',additional_topologies=[]):
        logger.debug(f'psfgen {self.basename} {psf} {pdb} {addpdb}')
        pg=self.writers['psfgen']
        pg.newscript(self.basename,additional_topologies=additional_topologies)
        pg.usescript('memb')
        pg.writescript(self.basename,guesscoord=False,regenerate=True,force_exit=True)
        result=pg.runscript(psf=psf,pdb=pdb,addpdb=addpdb,o=self.basename)
        return result
    
    # def membrane_embed(self,pdb,zgroups,zval,outpdb='embed.pdb',fixed_orientation=False):
    #     logger.debug(f'zgroups {zgroups}')
    #     ztopg,zbotg,zrefg=zgroups
    #     vm=self.writers['vmd']
    #     vm.newscript(self.basename,packages=['Orient','PestiferUtil'])
    #     vm.addline(f'mol new {pdb}')
    #     vm.addline(f'set TMP [atomselect top all]')
    #     vm.addline(f'set REF  [atomselect top "{zrefg}"]')
    #     if not fixed_orientation:
    #         vm.addline(f'set HEAD [atomselect top "{ztopg}"]')
    #         vm.addline(f'set TAIL [atomselect top "{zbotg}"]')
    #         vm.addline(f'vmdcon -info "[measure center $HEAD]"')
    #         vm.addline(f'vmdcon -info "[measure center $TAIL]"')
    #         vm.addline(f'set VEC [vecsub [measure center $HEAD] [measure center $TAIL]]')
    #         vm.addline( '$TMP move [orient $TMP $VEC {0 0 1}]')
    #     vm.addline(f'set COM [measure center $TMP]')
    #     vm.addline(f'$TMP moveby [vecscale $COM -1]')
    #     vm.addline(f'set REFZ [lindex [measure center $REF] 2]')
    #     vm.addline(f'$TMP moveby [list 0 0 [expr {zval}-($REFZ)]]')
    #     vm.addline(f'$TMP writepdb {outpdb}')
    #     vm.addline(f'set minmax [measure minmax $TMP]')
    #     vm.addline(f'set SEL_R [expr max([MaxRad $HEAD],[MaxRad $TAIL],[MaxRad $REF])]')
    #     vm.addline(f'set COM [measure center $TMP]')
    #     vm.addline(f'set f [open "{self.basename}-results.yaml" "w"]')
    #     vm.addline(r'puts $f "maxrad: [format %.5f $SEL_R]"')
    #     vm.addline(r'puts $f "lower_left: \[ [join [lindex $minmax 0] ,\ ] \]"')
    #     vm.addline(r'puts $f "upper_right: \[ [join [lindex $minmax 1] ,\ ] \]"')
    #     vm.addline(r'puts $f "com: \[ [join $COM ,\ ] \]"')
    #     vm.addline(r'puts $f "head_z: [format %.5f [lindex [measure center $HEAD] 2]]"')
    #     vm.addline(r'puts $f "tail_z: [format %.5f [lindex [measure center $TAIL] 2]]"')
    #     vm.addline(f'close $f')
    #     vm.writescript()
    #     vm.runscript()
