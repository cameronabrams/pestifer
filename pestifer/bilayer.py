# Author: Cameron F. Abrams, <cfa2@drexel.edu>
import logging
import numpy as np

from psfutil.psfcontents import PSFContents
from util.units import _UNITS_, _SYMBOLS_
from util.units import cuA_of_nmolec, nmolec_in_cuA

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
                leaflet_patch_nlipids=100,solvent_specstring='TIP3',solvent_ratio_specstring='1.0',solvent_to_lipid_ratio_specstring='8.0',
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
                'upper_chamber':[{'name':n,'frac':x,'patn':int(x*leaflet_patch_nlipids*s),'MW':resi_database['water_ions'][n].mass(),'charge':resi_database['water_ions'][n].charge} for n,x,s in zip(uc_s,uc_x,uc_slr)],
                'lower_chamber':[{'name':n,'frac':x,'patn':int(x*leaflet_patch_nlipids*s),'MW':resi_database['water_ions'][n].mass(),'charge':resi_database['water_ions'][n].charge} for n,x,s in zip(lc_s,lc_x,lc_slr)],
                'upper_leaflet':[{'name':n,'frac':x,'conf':c,'patn':int(x*leaflet_patch_nlipids),'charge':resi_database['lipid'][n].charge} for n,x,c in zip(ul_lip,ul_x,ul_c)],
                'lower_leaflet':[{'name':n,'frac':x,'conf':c,'patn':int(x*leaflet_patch_nlipids),'charge':resi_database['lipid'][n].charge} for n,x,c in zip(ll_lip,ll_x,ll_c)],
            }
        else:
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
                uc_s,lc_s=bilayer_stringsplit(solvent_specstring)
                uc_x,lc_x=bilayer_stringsplit(solvent_ratio_specstring,return_type=float)
                uc_xs,lc_xs=np.sum(uc_x),np.sum(lc_x)
                uc_x/=uc_xs
                lc_x/=lc_xs
                if 'upper_chamber' not in composition_dict:
                    composition_dict['upper_chamber']=[{'name':n,'frac':x,'MW':resi_database['water_ions'][n].mass(),'charge':resi_database['water_ions'][n].charge} for n,x in zip(uc_s,uc_x)]
                if 'lower_chamber' not in composition_dict:
                    composition_dict['lower_chamber']=[{'name':n,'frac':x,'MW':resi_database['water_ions'][n].mass(),'charge':resi_database['water_ions'][n].charge} for n,x in zip(lc_s,lc_x)]

            uc_slr,lc_slr=bilayer_stringsplit(solvent_to_lipid_ratio_specstring,return_type=float)
            logger.debug(f'{uc_slr} {type(uc_slr)} {lc_slr} {type(lc_slr)}')
            for c,slr in zip(['upper_chamber','lower_chamber'],[uc_slr,lc_slr]):
                L=composition_dict[c]
                logger.debug(f'{c} {L}')
                for d,dlr in zip(L,slr):
                    logger.debug(f'{d} {type(d)} {d["frac"]} {type(d["frac"])} {d["name"]} {type(d["name"])} {slr} {type(slr)}')
                    if not 'patn' in d:
                        d['patn']=int(d['frac']*leaflet_patch_nlipids*dlr)
                    if not 'charge' in d:
                        d['charge']=resi_database['water_ions'][d['name']].charge
                    if not 'MW' in d:
                        d['MW']=resi_database['water_ions'][d['name']].mass()

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
                self.species_data['local_name']=self.species_data[species_name].checkout()

    def build_patch(self,SAPL=60.0,xy_aspect_ratio=1.0,midplane_z=0.0,half_mid_zgap=1.0):
        patch_area=SAPL*self.leaflet_patch_nlipids
        Lx=np.sqrt(patch_area/xy_aspect_ratio)
        Ly=xy_aspect_ratio*Lx
        self.patch_area=patch_area
        lcvol=cuA_of_nmolec(self.LC['avgMW'],1.0,self.LC['patn'])
        ucvol=cuA_of_nmolec(self.UC['avgMW'],1.0,self.UC['patn'])
        lcdepth=lcvol/self.patch_area
        ucdepth=ucvol/self.patch_area
        zmin=midplane_z-self.LL['maxlength']-lcdepth-half_mid_zgap
        zmax=midplane_z+self.UL['maxlength']+ucdepth+half_mid_zgap
        self.LC['z-lo']=zmin
        self.LC['z-hi']=midplane_z-lcdepth-half_mid_zgap
        self.LL['z-lo']=midplane_z-lcdepth-half_mid_zgap
        self.LL['z-hi']=midplane_z-half_mid_zgap
        self.UL['z-lo']=midplane_z+half_mid_zgap
        self.UL['z-hi']=midplane_z+ucdepth+half_mid_zgap
        self.UC['z-lo']=midplane_z+ucdepth+half_mid_zgap
        self.UC['z-hi']=zmax
        self.patch_ll_corner=np.array([0,0,zmin])
        self.patch_ur_corner=np.array([Lx,Ly,zmax])
