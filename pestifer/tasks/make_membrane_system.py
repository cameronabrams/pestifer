# Author: Cameron F. Abrams, <cfa2@drexel.edu>
import logging
# import yaml

# import numpy as np

from copy import deepcopy

from ..bilayer import Bilayer, specstrings_builddict
from ..basetask import BaseTask
from ..charmmffcontent import CHARMMFFResiDatabase
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
        self.pdbrepository=self.config.RM.charmmff_content.pdbrepository
        self.charmmff_content=self.config.RM.charmmff_content
        self.RDB=CHARMMFFResiDatabase(self.charmmff_content,streamIDs=[])
        self.RDB.add_stream('lipid')
        self.RDB.add_topology('toppar_all36_moreions.str',streamIDoverride='water_ions')
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
        patch_nlipids=self.bilayer_specs.get('patch_nlipids',dict(upper=100,lower=100))
        cation_name=self.bilayer_specs.get('cation','POT')
        anion_name=self.bilayer_specs.get('anion','CLA')
        neutralizing_salt=[cation_name,anion_name]
        salt_con=self.bilayer_specs.get('salt_con',0.0)  # Molar concentration
        composition_dict=self.bilayer_specs.get('composition',{})

        if not composition_dict['upper_leaflet'] or not composition_dict['lower_leaflet']:
            logger.debug('No upper or lower leaflet specified in composition; building from memgen-format specstrings')
            composition_dict=specstrings_builddict(lipid_specstring,
                                                  ratio_specstring,
                                                  conformers_specstring,
                                                  solvent_specstring,
                                                  solvent_ratio_specstring)
        logger.debug(f'Main composition dict {composition_dict}')
        self.patch=Bilayer(composition_dict,
                            neutralizing_salt=neutralizing_salt,
                            salt_concentration=salt_con,
                            solvent_specstring=solvent_specstring,
                            solvent_ratio_specstring=solvent_ratio_specstring,
                            solvent_to_key_lipid_ratio=solvent_to_lipid_ratio,
                            leaflet_nlipids=patch_nlipids,
                            pdbrepository=self.pdbrepository,resi_database=self.RDB)
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
                                salt_concentration=salt_con,
                                solvent_specstring=solvent_specstring,
                                solvent_ratio_specstring=solvent_ratio_specstring,
                                solvent_to_key_lipid_ratio=solvent_to_lipid_ratio,
                                leaflet_nlipids=patch_nlipids,
                                pdbrepository=self.pdbrepository,resi_database=self.RDB)
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
                                leaflet_nlipids=patch_nlipids,
                                pdbrepository=self.pdbrepository,resi_database=self.RDB)
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
            self.make_quilt_from_patch()
        self.embed_protein()
        self.log_message('complete')
        return super().do()

    def build_patch(self):
        logger.debug(f'Bilayer specs: {self.bilayer_specs}')
        solution_gcc=self.bilayer_specs.get('solution_gcc',1.0)
        rotation_pm=self.bilayer_specs.get('rotation_pm',10.)
        half_mid_zgap=self.bilayer_specs.get('half_mid_zgap',1.0)
        SAPL=self.bilayer_specs.get('SAPL',75.0)
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
                              rotation_pm=rotation_pm,solution_gcc=solution_gcc,
                              half_mid_zgap=half_mid_zgap)
            pm=PackmolInputWriter(self.config)
            packmol_output_pdb=patch.pack_patch(pm,specname,seed=seed,
                                                tolerance=tolerance,
                                                nloop_all=nloop_all,
                                                half_mid_zgap=half_mid_zgap,
                                                rotation_pm=rotation_pm,
                                                nloop=nloop)
            self.next_basename(f'patch{spec}-build')
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
            patch.equilibrate(user_dict=deepcopy(self.config['user']),
                              basename=f'patch{spec}',index=self.index,
                              relaxation_protocol=relaxation_protocol,
                              parent_controller_index=self.controller_index)

    def make_quilt_from_patch(self):
        logger.debug(f'Creating quilt from patch')
        self.next_basename('quilt')
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
        relaxation_protocol=self.bilayer_specs.get('relaxation_protocols',{}).get('quilt',{})
        self.quilt.equilibrate(user_dict=deepcopy(self.config['user']),
                                basename='quilt',
                                relaxation_protocol=relaxation_protocol,
                                parent_controller_index=self.controller_index)

    def embed_protein(self):
        if not self.embed_specs:
            logger.debug('No embed specs.')
            return
        no_orient=self.embed_specs.get('no_orient',False)
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
                            no_orient=no_orient,
                            o=self.basename)
        self.statevars['pdb']=f'{self.basename}.pdb'
        self.statevars['psf']=f'{self.basename}.psf'
        self.statevars['coor']=f'{self.basename}.coor'
        self.statevars['xsc']=f'{self.basename}.xsc'
        if 'vel' in self.statevars:
            del self.statevars['vel']
        if 'charmmff_paramfiles' not in self.statevars:
            self.statevars['charmmff_paramfiles']=[]
        self.statevars['charmmff_paramfiles']+=self.quilt.addl_streamfiles
        self.statevars['charmmff_paramfiles']=list(set(self.statevars['charmmff_paramfiles']))
        self.quilt.statevars.update(self.statevars)
        return result

    
