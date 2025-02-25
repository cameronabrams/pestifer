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
from ..util.util import cell_to_xsc, nmolec_in_cuA

sA_ =_SYMBOLS_['ANGSTROM']
sA2_=_UNITS_['SQUARE-ANGSTROMS']
sA3_=_UNITS_['CUBIC-ANGSTROMS']

logger=logging.getLogger(__name__)

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
        if psf is not None and pdb is not None:
            logger.debug(f'BilayerEmbedTask will use psf {psf} and pdb {pdb} as inputs')

        self.next_basename('bilayer')

        lipid_specstring=self.specs['lipids']
        ratio_specstring=self.specs['mole_fractions']
        conformers_specstring=self.specs['lipid_conformers']

        solvent_specstring=self.specs.get('solvents','TIP3')
        solv_molfrac_specstring=self.specs.get('solvent_mole_fractions','1.0')
        cation_name=self.specs.get('cation','POT')
        cation_q=RDB['water_ions'][cation_name].charge
        anion_name=self.specs.get('anion','CLA')
        anion_q=RDB['water_ions'][anion_name].charge
        assert cation_name in RDB['water_ions'],f'Cation {cation_name} not found in available CHARMM topologies'
        assert anion_name in RDB['water_ions'],f'Anion {anion_name} not found in available CHARMM topologies'
        cation_pdbstruct=self.pdb_collection.get_pdb(cation_name)
        cation_pdbstruct.checkout()
        anion_pdbstruct=self.pdb_collection.get_pdb(anion_name)
        anion_pdbstruct.checkout()
        salt_con_M=self.specs.get('salt_con',0.0)

        rotation_pm=self.specs.get('rotation_pm',10.)
        fuzz_factor=self.specs.get('fuzz_factor',0.5)

        SAPL=self.specs.get('SAPL',60.0)
        scale_excluded_volume=self.specs.get('scale_excluded_volume',1.0)
        seed=self.specs.get('seed',27021972)
        tolerance=self.specs.get('tolerance',2.0)
        solution_gcc=self.specs.get('solution_gcc',1.0)
        leaflet_thickness=self.specs.get('leaflet_thickness',22.0)
        nloop=self.specs.get('nloop',200)
        nloop_all=self.specs.get('nloop_all',200)
        boxdim=np.array(list(map(float,self.specs.get('dims',[]))))

        self.slices={'LOWER-CHAMBER':{},'LOWER-LEAFLET':{},'UPPER-LEAFLET':{},'UPPER-CHAMBER':{}}
        LC=self.slices['LOWER-CHAMBER']
        LL=self.slices['LOWER-LEAFLET']
        UL=self.slices['UPPER-LEAFLET']
        UC=self.slices['UPPER-CHAMBER']
        
        if 'embed' in self.specs:
            assert pdb is not None
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
            box_area=boxdim[0]*boxdim[1]
            protein_xyarea=np.pi*(protein_rad*prot_rad_scal)**2
            logger.debug(f'box area {box_area:.3f} {sA2_}')
            logger.debug(f'protein radius {protein_rad:.3f} {sA_} xyarea {protein_xyarea:.3f} {sA2_}')
            mem_area=box_area-protein_xyarea
            logger.debug(f'protein excludes approx. {np.floor(protein_xyarea/SAPL)} lipids per leaflet')
            logger.debug(f'actual membrane area {mem_area:.3f} {sA2_}')
            LC['z-lo']= box_ll[2]
            LC['z-hi']= midplane-leaflet_thickness
            LL['z-lo']= midplane-leaflet_thickness
            LL['z-hi']= midplane
            UL['z-lo']= midplane
            UL['z-hi']= midplane+leaflet_thickness
            UC['z-lo']= midplane+leaflet_thickness
            UC['z-hi']= box_ur[2]
        else:
            box_ll=np.zeros(3,dtype=float)
            box_ur=boxdim
            origin=0.5*(box_ur+box_ll)
            box_area=boxdim[0]*boxdim[1]
            mem_area=box_area
            LC['z-lo']= box_ll[2]
            LC['z-hi']= origin[2]-leaflet_thickness
            LL['z-lo']= origin[2]-leaflet_thickness
            LL['z-hi']= origin[2]
            UL['z-lo']= origin[2]
            UL['z-hi']= origin[2]+leaflet_thickness
            UC['z-lo']= origin[2]+leaflet_thickness
            UC['z-hi']= box_ur[2]

        for S in [LC,LL,UL,UC]:
            S['THICKNESS']=S['z-hi']-S['z-lo']
            S['INIT-VOLUME']=S['THICKNESS']*box_area
            S['INIT-NWATEREQUIV']=nmolec_in_cuA(18.0,1.0,S['INIT-VOLUME'])
        
        my_logger(self.slices,logger.debug)

        boxV=boxdim.prod()
        logger.debug(f'box corners ll {box_ll} ur {box_ur}')
        logger.debug(f'box dimensions {boxdim}')
        logger.debug(f'box area {box_area} {sA2_}; membrane_area {mem_area} {sA2_}')
        logger.debug(f'box volume {boxV:.3f} {sA3_} (check sum-over-slices: {np.sum(np.array([x["INIT-VOLUME"] for x in self.slices.values()]))} {sA3_}) (or {nmolec_in_cuA(18.0,1.0,boxV)} waters at 1 g/cc)')
        logger.debug(f'membrane volume {box_area*2*leaflet_thickness:.3f} {sA3_}')
        logger.debug(f'   lower chamber: {LC["INIT-VOLUME"]} {sA3_}; upper chamber: {UC["INIT-VOLUME"]} {sA3_}')
        logger.debug(f'due to membrane, expecting {nmolec_in_cuA(18.0,1.0,boxV-(box_area*2*leaflet_thickness))} waters')
        logger.debug(f'   lower chamber: {LC["INIT-NWATEREQUIV"]}; upper chamber: {UC["INIT-NWATEREQUIV"]}')

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
            solv_pdbstruct=self.pdb_collection.get_pdb(sn)
            solv_pdbstruct.checkout()
            solvent[sn]['mass']=RDB['water_ions'][sn].mass()

        my_logger(solvent,logger.debug)

        leaf_lipspec=lipid_specstring.split('//')
        logger.debug(f'leaf_lipspec {leaf_lipspec}')
        leaf_ratiospec=ratio_specstring.split('//')
        leaf_conformerspec=conformers_specstring.split('//')
        if leaf_conformerspec==['0']: # no value explicitly specified -- replicate 0's to match pattern
            if len(leaf_lipspec)>1 or any([x.count(':')>0 for x in leaf_lipspec]):
                conspec=[]
                for leaf in leaf_lipspec:
                    nlipnames=len(leaf.split(':'))
                    expspec=':'.join('0'*nlipnames)
                    logger.debug(f'nlipnames {nlipnames} expsepec {expspec}')
                    conspec.append(expspec)
                logger.debug(f'{leaf_conformerspec}')
                leaf_conformerspec=conspec
        else:
            leaf_conformerspec=conformers_specstring.split('//')

        assert len(leaf_lipspec)==len(leaf_ratiospec),f'lipid names and mole fractions are not congruent'
        assert len(leaf_lipspec)==len(leaf_conformerspec),f'lipid names and conformer indices are not congruent {leaf_lipspec} {leaf_conformerspec}'
        if len(leaf_lipspec)==1:  # symmetrical bilayer
            leaf_lipspec.append(leaf_lipspec[0])
            leaf_ratiospec.append(leaf_ratiospec[0])
            leaf_conformerspec.append(leaf_conformerspec[0])
        
        global_lipid_names=[]
        # build leaflets
        LL['lipids']=[]
        UL['lipids']=[]
        for li,(leaflet_name,lipid_name_string,lipid_molfrac_string,lipid_conformer_string) in enumerate(zip([LL,UL],leaf_lipspec,leaf_ratiospec,leaf_conformerspec)):
            lipid_names=lipid_name_string.split(':')
            lipid_conformers=np.array(list(map(int,lipid_conformer_string.split(':'))))
            assert len(lipid_names)==len(lipid_conformers),f'lipid names and conformer indices are not congruent within a leaflet specification'
            lipid_molfracs=np.array(list(map(float,lipid_molfrac_string.split(':'))))
            sumlm=np.sum(lipid_molfracs)
            lipid_molfracs/=sumlm
            assert len(lipid_names)==len(lipid_molfracs),f'lipid names and mole fractions are not congruent within a leaflet specification'
            for ll in lipid_names:
                if not ll in global_lipid_names: global_lipid_names.append(ll)
            leaflet_name['lipids']=[dict(name=name,frac=frac,conformer=conformer) for name,frac,conformer in zip(lipid_names,lipid_molfracs,lipid_conformers)]

        logger.debug(f'leaflets')
        my_logger(LL,logger.debug)
        my_logger(UL,logger.debug)

        lipid_data={}
        addl_streamfiles=[]
        for l in global_lipid_names:
            pdbstruct=self.pdb_collection.get_pdb(l)
            if pdbstruct==None:
                logger.error(f'No PDB available for lipid {l}')
            lipid_data[l]=pdbstruct
            for p in lipid_data[l].get_parameters():
                if p.endswith('.str') and not p in addl_streamfiles:
                    addl_streamfiles.append(p)

        tot_lip_mols=mem_area/SAPL
        bilayer_charge=0.0
        for leaflet in [LL,UL]:
            for component in leaflet['lipids']:
                component['n']=int(np.floor(component['frac']*tot_lip_mols))
                bilayer_charge+=lipid_data[component['name']].get_charge()*component['n']
                component['local_name']=lipid_data[component['name']].checkout(index=component['conformer'])
        global_charge=pro_charge+bilayer_charge
        my_logger(self.slices,logger.debug)
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
            for specs in leaflet['lipids']:
                name=specs['name']
                ref_atoms=lipid_data[name].get_ref_atoms()
                hs=' '.join([f"{x['serial']}" for x in ref_atoms['heads']])
                ts=' '.join([f"{x['serial']}" for x in ref_atoms['tails']])
                n=specs['n']
                pm.addline(f'structure {specs["local_name"]}')
                pm.addline(f'number {n}',indents=1)
                length=lipid_data[name].get_ref_length(index=specs['conformer'])
                Dz=np.cos(np.deg2rad(rotation_pm))*length
                fuzz=Dz-length
                fuzz_out,fuzz_in=fuzz*fuzz_factor,fuzz*(1-fuzz_factor)
                leaflet_thickness=leaflet['z-hi']-leaflet['z-lo']
                logger.debug(f'{name}: {specs["local_name"]} length {length} fuzz_in {fuzz_in:.3f} fuzz_out {fuzz_out:.3f}')
                if length>leaflet_thickness:
                    if leaflet is LL:
                        below_plane_z=leaflet['z-hi']-fuzz_out-length
                        above_plane_z=leaflet['z-hi']+fuzz_in
                        below_plane_atoms=hs
                        above_plane_atoms=ts
                    elif leaflet is UL:
                        below_plane_z=leaflet['z-lo']-fuzz_in
                        above_plane_z=leaflet['z-lo']+fuzz_out+length   
                        below_plane_atoms=ts
                        above_plane_atoms=hs
                    pm.addline(f'inside box {box_ll[0]:.3f} {box_ll[1]:.3f} {box_ll[2]:.3f} {box_ur[0]:.3f} {box_ur[1]:.3f} {box_ur[2]:.3f}',indents=1)
                    pm.addline(f'atoms {below_plane_atoms}',indents=1)
                    pm.addline(f'below plane 0. 0. 1. {below_plane_z:.3f}',indents=2)
                    pm.addline( 'end atoms',indents=1)
                    pm.addline(f'atoms {above_plane_atoms}',indents=1)
                    pm.addline(f'above plane 0. 0. 1. {above_plane_z:.3f}',indents=2)
                    pm.addline( 'end atoms',indents=1)
                else:
                    if leaflet is LL:
                        constrain_rotation=180.0
                    elif leaflet is UL:
                        constrain_rotation=0.0
                    inside_z_lo=leaflet['z-lo']
                    inside_z_hi=leaflet['z-hi']
                    pm.addline(f'inside box {box_ll[0]:.3f} {box_ll[1]:.3f} {inside_z_lo:.3f} {box_ur[0]:.3f} {box_ur[1]:.3f} {inside_z_hi:.3f}',indents=1)
                    pm.addline(f'constrain_rotation x {constrain_rotation} {rotation_pm}',indents=1)
                    pm.addline(f'constrain_rotation y {constrain_rotation} {rotation_pm}',indents=1)
                pm.addline(f'nloop {nloop}',indents=1)
                pm.addline(f'end structure')
        pm.writefile()
        self.result=pm.runscript()
        logger.debug(f'first packmol result {self.result}')
        if self.result!=0:
            return super().do()
        
        anion_qtot=cation_qtot=0
        if sg=='+':
            anion_qtot=int(global_charge)
        else:
            cation_qtot=int(np.abs(np.round(global_charge)))

        anion_n=int(anion_qtot//np.abs(int(anion_q)))
        cation_n=int(cation_qtot//np.abs(int(cation_q)))
        ions={anion_name:{'n':anion_n,'claimed':0},cation_name:{'n':cation_n,'claimed':0}}
        my_logger(ions,logger.debug)

        # computes the volumes in the upper and lower chambers that are occupied
        # by lipid and/or protein atoms and adjusts the apparent chamber volumes
        # accordingly
        cdf=coorddf_from_pdb(packmol_output_pdb)
        density_calc_resolution=2. # Angstrom
        nxbins=int(boxdim[0]/density_calc_resolution)
        nybins=int(boxdim[1]/density_calc_resolution)
        for S in [LC,UC]:
            slice_atoms=cdf[(cdf['z']<S['z-hi'])&(cdf['z']>S['z-lo'])]
            my_logger(slice_atoms.head(),logger.debug)
            logger.debug(f'{slice_atoms.shape[0]} atoms')
            S['AVAILABLE-VOLUME']=S['INIT-VOLUME']
            if slice_atoms.shape[0]>0:
                nzbins=int((S['THICKNESS'])/density_calc_resolution)
                h,edges=np.histogramdd(slice_atoms[['x','y','z']].to_numpy(),bins=(nxbins,nybins,nzbins))
                binV=boxdim[0]*boxdim[1]*(S['THICKNESS'])/(nxbins*nybins*nzbins)
                logger.debug(f'bin volume {binV:.3f} {sA3_}')
                hit_bin_count=h[h>0.0].size
                logger.debug(f'{hit_bin_count}/{h.size} bins have some atom density')
                excluded_volume=hit_bin_count*binV
                logger.debug(f'excluded volume {excluded_volume} {sA3_}')
                S['EXCLUDED-VOLUME']=excluded_volume
                S['AVAILABLE-VOLUME']-=S['EXCLUDED-VOLUME']*scale_excluded_volume
                logger.debug(f'available volume {S["AVAILABLE-VOLUME"]} {sA3_} ({nmolec_in_cuA(18.0,1.0,S["AVAILABLE-VOLUME"])} water-equivs)')

        my_logger(self.slices,logger.debug)

        total_available_volume=LC['AVAILABLE-VOLUME']+UC['AVAILABLE-VOLUME']
        logger.debug(f'total volume available for solvent {total_available_volume:.3f} {sA3_} ({total_available_volume/boxV*100:2f} of box)')
        for chamber in [LC,UC]:
            chamber['VOLUME-FRAC']=chamber['AVAILABLE-VOLUME']/total_available_volume
            chamber['comp']={}
            for sn,sspec in solvent.items():
                mass=sspec['mass']
                x=sspec['mol-frac']
                chamber['comp'][sn]=int(x*nmolec_in_cuA(mass,solution_gcc,chamber['AVAILABLE-VOLUME']))
                logger.debug(f'solvent component {sn}: mass {mass} x {x} n {chamber["comp"][sn]}')
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

        my_logger(self.slices,logger.debug)

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
        self.result=self.psfgen(psf=psf,pdb=pdb,addpdb=packmol_output_pdb,additional_topologies=addl_streamfiles)
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
