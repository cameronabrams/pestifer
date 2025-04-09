# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import os 

import numpy as np

from .util.units import _UNITS_, _SYMBOLS_, cuA_of_nmolec

sA_ =_SYMBOLS_['ANGSTROM']
sA2_=_UNITS_['SQUARE-ANGSTROMS']
sA3_=_UNITS_['CUBIC-ANGSTROMS']

logger=logging.getLogger(__name__)

def bilayer_stringsplit(input_string,delimiter0='//',delimiter1=':',return_type=str,normalize=False,symmetrize=None):
    """Splits a string into two lists of strings, separated by the first delimiter.
    The two lists are separated by the second delimiter.
    The first delimiter is used to split the string into two parts, and the second delimiter
    is used to split each part into a list of strings.
    The function returns two lists of strings.
    If the input string is empty, the function returns two empty lists. """
    logger.debug(f'Processing \'{input_string}\'')
    uls,lls=(input_string.split(delimiter0)+[input_string])[:2]
    ul,ll=uls.split(delimiter1),lls.split(delimiter1)
    if symmetrize=='upper':
        ll=ul
    elif symmetrize=='lower':
        ul=ll
    if return_type==str:
        return ul,ll
    ul,ll=list(map(return_type,ul)),list(map(return_type,ll))
    ul,ll=np.array(ul),np.array(ll)
    if normalize:
        ul_xs,ll_xs=np.sum(ul),np.sum(ll)
        if ul_xs>0.0:
            ul/=ul_xs
        if ll_xs>0.0:
            ll/=ll_xs
    return ul,ll

class Bilayer:
    def __init__(self,composition_dict={},lipid_specstring='',lipid_ratio_specstring='',lipid_conformers_specstring='',    
                leaflet_patch_nlipids=100,solvent_specstring='TIP3',solvent_ratio_specstring='1.0',solvent_to_lipid_ratio_specstring='32.0',
                neutralizing_salt=['POT','CLA'],pdb_collection=None,resi_database=None,symmetrize=None):

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
        logger.debug(f'Passed in solvent_specstring \'{solvent_specstring}\' and solvent_ratio_specstring \'{solvent_ratio_specstring}\'')
        if not composition_dict:
            # old-style bilayer composition specification with memgen-style specstrings
            ul_lip,ll_lip=bilayer_stringsplit(lipid_specstring,symmetrize=symmetrize)
            ul_x,ll_x    =bilayer_stringsplit(lipid_ratio_specstring,return_type=float,normalize=True,symmetrize=symmetrize)
            ul_c,ll_c    =bilayer_stringsplit(lipid_conformers_specstring,return_type=int,symmetrize=symmetrize)
            uc_s,lc_s    =bilayer_stringsplit(solvent_specstring,symmetrize=symmetrize)
            uc_x,lc_x    =bilayer_stringsplit(solvent_ratio_specstring,return_type=float,normalize=True,symmetrize=symmetrize)
            uc_slr,lc_slr=bilayer_stringsplit(solvent_to_lipid_ratio_specstring,return_type=float,symmetrize=symmetrize)
            assert len(ul_lip)==len(ul_x),f'Upper leaflet has {len(ul_lip)} lipids but {len(ul_x)} mole fractions specified'
            assert len(ul_lip)==len(ul_c),f'Upper leaflet has {len(ul_lip)} lipids but {len(ul_c)} conformers specified'
            assert len(ul_x)==len(ul_c),f'Upper leaflet has {len(ul_x)} mole fractions but {len(ul_c)} conformers specified'
            assert len(ll_lip)==len(ll_x),f'Upper leaflet has {len(ul_lip)} lipids but {len(ul_x)} mole fractions specified'
            assert len(ll_lip)==len(ll_c),f'Upper leaflet has {len(ul_lip)} lipids but {len(ul_c)} conformers specified'
            assert len(ll_x)==len(ll_c),f'Upper leaflet has {len(ul_x)} mole fractions but {len(ul_c)} conformers specified'
            # build the equivalent new-style composition dictionary; 'patn' is the number of molecules for a minimal patch
            composition_dict={
                'upper_chamber':[{'name':n,'frac':x,'patn':int(x*leaflet_patch_nlipids*s),'MW':resi_database['water_ions'][n].mass(),'charge':resi_database['water_ions'][n].charge} for n,x,s in zip(uc_s,uc_x,uc_slr)],
                'lower_chamber':[{'name':n,'frac':x,'patn':int(x*leaflet_patch_nlipids*s),'MW':resi_database['water_ions'][n].mass(),'charge':resi_database['water_ions'][n].charge} for n,x,s in zip(lc_s,lc_x,lc_slr)],
                'upper_leaflet':[{'name':n,'frac':x,'conf':c,'patn':int(x*leaflet_patch_nlipids),'charge':resi_database['lipid'][n].charge} for n,x,c in zip(ul_lip,ul_x,ul_c)],
                'lower_leaflet':[{'name':n,'frac':x,'conf':c,'patn':int(x*leaflet_patch_nlipids),'charge':resi_database['lipid'][n].charge} for n,x,c in zip(ll_lip,ll_x,ll_c)],
            }
        else: # composition dictionary is provided
            logger.debug(f'Composition dictionary for bilayer is provided')
            ul_lip,ll_lip=[x['name'] for x in composition_dict['upper_leaflet']],[x['name'] for x in composition_dict['lower_leaflet']]
            ul_x,ll_x=[x['frac'] for x in composition_dict['upper_leaflet']],[x['frac'] for x in composition_dict['lower_leaflet']]
            # new-style bilayer composition specification with composition dictionary
            for l in ['upper_leaflet','lower_leaflet']:
                L=composition_dict[l]
                for d in L:
                    if not 'patn' in d:
                        d['patn']=int(d['frac']*leaflet_patch_nlipids)
                    if not 'charge' in d:
                        d['charge']=resi_database['lipid'][d['name']].charge
                    if not 'MW' in d:
                        d['MW']=resi_database['lipid'][d['name']].mass()
            # user need not have specified the solvent composition in the upper and lower chambers
            if 'upper_chamber' not in composition_dict or 'lower_chamber' not in composition_dict:
                logger.debug(f'Provided composition dictionary does not include solvent chamber specifications')
                logger.debug(f'Using specstrings \'{solvent_specstring}\' and \'{solvent_ratio_specstring}\' to build solvent composition')
                uc_s,lc_s=bilayer_stringsplit(solvent_specstring,symmetrize=symmetrize)
                uc_x,lc_x=bilayer_stringsplit(solvent_ratio_specstring,return_type=float,normalize=True,symmetrize=symmetrize)
                if 'upper_chamber' not in composition_dict:
                    composition_dict['upper_chamber']=[{'name':n,'frac':x,'MW':resi_database['water_ions'][n].mass(),'charge':resi_database['water_ions'][n].charge} for n,x in zip(uc_s,uc_x)]
                if 'lower_chamber' not in composition_dict:
                    composition_dict['lower_chamber']=[{'name':n,'frac':x,'MW':resi_database['water_ions'][n].mass(),'charge':resi_database['water_ions'][n].charge} for n,x in zip(lc_s,lc_x)]

            uc_slr,lc_slr=bilayer_stringsplit(solvent_to_lipid_ratio_specstring,return_type=float,symmetrize=symmetrize)
            for c,slr in zip(['upper_chamber','lower_chamber'],[uc_slr,lc_slr]):
                L=composition_dict[c]
                logger.debug(f'{c} {L}')
                for d,dlr in zip(L,slr):
                    if not 'patn' in d:
                        d['patn']=int(d['frac']*leaflet_patch_nlipids*dlr)
                    if not 'charge' in d:
                        d['charge']=resi_database['water_ions'][d['name']].charge
                    if not 'MW' in d:
                        d['MW']=resi_database['water_ions'][d['name']].mass()


        if symmetrize=='upper':
            logger.debug(f'Symmetrizing bilayer to upper leaflet')
            composition_dict['lower_leaflet']=composition_dict['upper_leaflet']
            composition_dict['lower_chamber']=composition_dict['upper_chamber']
        elif symmetrize=='lower':
            logger.debug(f'Symmetrizing bilayer to lower leaflet')
            composition_dict['upper_leaflet']=composition_dict['lower_leaflet']
            composition_dict['upper_chamber']=composition_dict['lower_chamber']

        # if the bilayer is asymmetric (each leaflet has a unique composition), we cannot assume a priori that
        # each leaflet in a patch has the same number of lipids.  We set a flag to indicate that the patch is 
        # asymmetric and that the number of lipids in each leaflet may need to be adjusted after equilibration
        # and measurment of the pressure profile.
        ul_lx,ll_lx=[(x['name'],x['frac']) for x in composition_dict['upper_leaflet']],[(x['name'],x['frac']) for x in composition_dict['lower_leaflet']]
        self.asymmetric=set(ul_lx)!=set(ll_lx)
        if symmetrize is not None:
            assert self.asymmetric==False,f'You have specified a symmetric bilayer but the upper and lower leaflets are not the same: {ul_lx} {ll_lx}'

        lipid_names=[x['name'] for x in composition_dict['upper_leaflet']]
        lipid_names+= [x['name'] for x in composition_dict['lower_leaflet']]
        self.lipid_names=list(set(lipid_names))
        
        solvent_names=[x['name'] for x in composition_dict['upper_chamber']]
        solvent_names+= [x['name'] for x in composition_dict['lower_chamber']]
        self.solvent_names=list(set(solvent_names))
        self.species_names=self.lipid_names+self.solvent_names

        if pdb_collection is not None:
            self.species_data={}
            self.addl_streamfiles=[]
            for l in self.species_names:
                logger.debug(f'Getting pdb for {l}')
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
                logger.debug(f'Checked out {species_name} as {species["local_name"]}')

    def build_patch(self,SAPL=60.0,xy_aspect_ratio=1.0,midplane_z=0.0,half_mid_zgap=1.0,solution_gcc=1.0,rotation_pm=10.0):
        patch_area=SAPL*self.leaflet_patch_nlipids
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
        self.origin=np.array([0.5*Lx,0.5*Ly,0.5*(zmin+zmax)])

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

    def delete_lipid(self,vm,count,leaflet='upper'):
        selstr="and z>0.0" if leaflet=='upper' else "and z<0.0"
        pdb=self.statevars['pdb']
        psf=self.statevars['psf']
        basenamepdb,dum=os.path.splitext(pdb)
        basenamepsf,dum=os.path.splitext(psf)
        vm.newscript(f'delete-lipid')
        vm.addline(f'package require psfgen')
        vm.addline(f'readpsf {psf} pdb {pdb}')
        vm.addline(f'mol new {psf}')
        vm.addline(f'mol addfile {pdb}')
        vm.addline(f'set sel [atomselect top "lipid {selstr}"]')
        vm.addline(f'set lres [lsort -unique [$sel get residue]]')
        vm.addline(f'set nres [llength $lres]')
        vm.addline(f'set indices [list]')
        vm.addline(r'for {set i 0} {$i < $len} {incr i} {')
        vm.addline(f'    lappend indices $i')
        vm.addline(r'}')
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