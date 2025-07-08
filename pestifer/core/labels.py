# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Labels and segtypes for residues
This module defines the segment types and residue names used in
CHARMM and PDB files, along with their mappings. It also provides
a class for managing these labels and mappings.
"""

segtypes= {
    'protein': {
        'macro': False,
        'resnames': [
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
            'HIS', 'HSP', 'HSD', 'HSE', 'ILE', 'LEU', 'LYS', 'MET',
            'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'],
        'rescodes': {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
            'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
            'HIS': 'H', 'HSD': 'H', 'HSE': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
            'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
            'TYR': 'Y', 'VAL': 'V', 'HSP': 'H'},
        'invrescodes': {
            'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
            'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
            'H': 'HSD', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
            'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
            'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}},
    'ion': {
        'macro': True,
        'resnames': [
            'LIT', 'SOD', 'MG' , 'POT', 'CAL', 'RUB', 'CES',
            'BAR', 'ZN' , 'CAD', 'CL',  'SO4', 'PO4', 'H2PO',
            'ZN2', 'CLA', 'FE']},
    'ligand': {
        'macro': True,
        'resnames' : [
            'EIC', 'VCG', '83G', 'HEM']},
    'nucleicacid': {
        'macro': True,
        'resnames': [
            'ADE', 'THY', 'CYT', 'GUA', 'URA', 'DA', 'DC', 'DG', 'DT', 'DU',]},
    'glycan': {
        'macro':True,
        'resnames': [
            'BGC'   , 'MAL'   , 'GLC'   , 'NDG'   , 'SGN'   , '4YS'   , 'NAG'   ,
            'GCS'   , 'GCU'   , 'QUI'   , 'OLI'   , 'BMA'   , 'MAN'   , 'BEM'   ,
            'MAV'   , 'RAM'   , 'TYV'   , 'AHR'   , 'ARA'   , 'GLA'   , 'GAL'   ,
            'NGA'   , 'ADA'   , 'GL0'   , 'GUP'   , 'GUL'   , 'LGU'   , 'ALT'   ,
            'WOO'   , 'ALL'   , 'TAL'   , 'IDO'   , 'IDS'   , 'FUL'   , 'FUC'   ,
            'LYX'   , 'ABE'   , 'XYP'   , 'LXC'   , 'XYS'   , 'XYL'   , 'PAR'   ,
            'RIB'   , 'DIG'   , 'COL'   , 'BAC'   , 'API'   , 'FRU'   , 'TAG'   ,
            'SOR'   , 'PSI'   , 'DHA'   , 'KDN'   , 'KDO'   , 'NEU'   , 'SIA'   ,
            'MUR'   , 'GMH'   , 'AGLC'  , 'BGLC'  , 'AALT'  , 'BALT'  , 'AALL'  ,
            'BALL'  , 'AGAL'  , 'BGAL'  , 'AGUL'  , 'BGUL'  , 'AIDO'  , 'BIDO'  ,
            'AMAN'  , 'BMAN'  , 'ATAL'  , 'BTAL'  , 'AXYL'  , 'BXYL'  , 'AFUC'  ,
            'BFUC'  , 'ARHM'  , 'BRHM'  , 'AGLCA' , 'BGLCA' , 'BGLCA0', 'AIDOA' ,
            'BIDOA' , 'AGLCNA', 'BGLCNA', 'BGLCN0', 'AGALNA', 'BGALNA', 'ANE5AC',
            'BNE5AC', 'ABEQ'  , 'ARHMOA', 'MGLYOL', 'MERYOL', 'DTHROL', 'LTHROL',
            'MRIBOL', 'DARAOL', 'LARAOL', 'MXYLOL', 'MALLOL', 'DALTOL', 'LALTOL',
            'DGLUOL', 'LGLUOL', 'DMANOL', 'LMANOL', 'DGULOL', 'LGULOL', 'DIDIOL',
            'LIDIOL', 'MGALOL', 'ALLOSE', 'PSICOS', 'INI1'  , 'INI2'  , 'INI3'  ,
            'INI4'  , 'INI5'  , 'ADEO'  , 'BDEO'  , 'ARIB'  , 'BRIB'  , 'AARB'  ,
            'BARB'  , 'ALYF'  , 'BLYF'  , 'AXYF'  , 'BXYF'  , 'AFRU'  , 'BFRU']},
        'lipid': {
            'macro': True,
            'resnames': [
            '23SM'  , '2HEX'  , '2PTE'  , '3HEX'  , 'ABLIPA', 'ABLIPB', 'ACHO'  , 
            'ACTF'  , 'ADR'   , 'ADRP'  , 'ALIN'  , 'ALINP' , 'APPC'  , 'ARA'   ,
            'ARAN'  , 'ARANP' , 'ARAP'  , 'ASM'   , 'BCLIPA', 'BCLIPB', 'BCLIPC',
            'BEH'   , 'BEHP'  , 'BSM'   , 'BTE1'  , 'BTE2'  , 'BUTA'  , 'C6DHPC',
            'C7DHPC', 'C6DH'  , 'C7DH'  , 'CER160', 'CER180', 'CER181', 'CER2'  ,
            'CER200', 'CER220', 'CER240', 'CER241', 'CER3E' , 'CHAPS' , 'CHAPSO',
            'CHL1'  , 'CHM1'  , 'CHNS'  , 'CHOL'  , 'CHSD'  , 'CHSP'  , 'CJLIPA',
            'CTLIPA', 'CYFOS3', 'CYFOS4', 'CYFOS5', 'CYFOS6', 'CYFOS7', 'CYSF'  ,
            'CYSG'  , 'CYSL'  , 'CYSP'  , 'DAPA'  , 'DAPC'  , 'DAPE'  , 'DAPG'  ,
            'DAPS'  , 'DCPC'  , 'DDA'   , 'DDAO'  , 'DDAOP' , 'DDAP'  , 'DDMG'  ,
            'DDOPC' , 'DDOPE' , 'DDOPS' , 'DDPC'  , 'DEPA'  , 'DEPC'  , 'DEPE'  ,
            'DEPG'  , 'DEPS'  , 'DGLA'  , 'DGLAP' , 'DGPA'  , 'DGPC'  , 'DGPE'  , 
            'DGPG'  , 'DGPS'  , 'DHA'   , 'DHAP'  , 'DHPCE' , 'DIHE'  , 'DIPA'  ,
            'DIPE'  , 'DIPG'  , 'DIPS'  , 'DLA'   , 'DLIPA' , 'DLIPC' , 'DLIPE' ,
            'DLPA'  , 'DLPC'  , 'DLPE'  , 'DLPG'  , 'DLPS'  , 'DLiPC' , 'DLiPE' ,
            'DLiPI' , 'DLIPI' , 'DMP'   , 'DMPA'  , 'DMPC'  , 'DMPCE' , 'DMPE'  ,
            'DMPEE' , 'DMPG'  , 'DMPI'  , 'DMPI13', 'DMPI14', 'DMPI15', 'DMPI24',
            'DMPI25', 'DMPI2A', 'DMPI2B', 'DMPI2C', 'DMPI2D', 'DMPI33', 'DMPI34',
            'DMPI35', 'DMPS'  , 'DNPA'  , 'DNPC'  , 'DNPE'  , 'DNPG'  , 'DNPS'  ,
            'DOMG'  , 'DOPA'  , 'DOPC'  , 'DOPCE' , 'DOPE'  , 'DOPEE' , 'DOPG'  ,
            'DOPP1' , 'DOPP2' , 'DOPP3' , 'DOPS'  , 'DPA'   , 'DPAP'  , 'DPPA'  ,
            'DPPC'  , 'DPPE'  , 'DPPEE' , 'DPPG'  , 'DPPGK' , 'DPPS'  , 'DPT'   ,
            'DPTP'  , 'DSPA'  , 'DSPC'  , 'DSPE'  , 'DSPG'  , 'DSPS'  , 'DTPA'  ,
            'DUPC'  , 'DXPC'  , 'DXPE'  , 'DYPA'  , 'DYPC'  , 'DYPE'  , 'DYPG'  ,
            'DYPS'  , 'ECLIPA', 'ECLIPB', 'ECLIPC', 'EDA'   , 'EDAP'  , 'EICO'  ,
            'EICOP' , 'EPA'   , 'EPAP'  , 'EPEN'  , 'ERG'   , 'ERU'   , 'ERUP'  ,
            'ETA'   , 'ETAC'  , 'ETAM'  , 'ETAP'  , 'ETE'   , 'ETEP'  , 'ETHE'  ,
            'FOIS11', 'FOIS9' , 'FOS10' , 'FOS12' , 'FOS13' , 'FOS14' , 'FOS15' ,
            'FOS16' , 'GH2F'  , 'GLA'   , 'GLAP'  , 'GLPH'  , 'GLYC'  , 'GLYM'  ,
            'GPC'   , 'HEPT'  , 'HEXA'  , 'HPA'   , 'HPAP'  , 'HPLIPA', 'HPLIPB',
            'HTA'   , 'HTAP'  , 'IPAC'  , 'IPPC'  , 'KPLIPA', 'KPLIPB', 'KPLIPC',
            'LAPAO' , 'LAPAOP', 'LAU'   , 'LAUP'  , 'LDAO'  , 'LDAOP' , 'LIGN'  ,
            'LIGNP' , 'LILIPA', 'LIN'   , 'LINP'  , 'LLPA'  , 'LLPC'  , 'LLPE'  ,
            'LLPS'  , 'LMPG'  , 'LNACL1', 'LNACL2', 'LNBCL1', 'LNBCL2', 'LNCCL1',
            'LNCCL2', 'LNDCL1', 'LNDCL2', 'LOACL1', 'LOACL2', 'LOCCL1', 'LOCCL2',
            'LPC12' , 'LPC14' , 'LPPC'  , 'LPPG'  , 'LSM'   , 'LYSM'  , 'MAS'   ,
            'MBUT'  , 'MCLIPA', 'MEA'   , 'MEAP'  , 'MEEF'  , 'MPRO'  , 'MP_1'  ,
            'MP_2'  , 'MSO4'  , 'MYR'   , 'MYRO'  , 'MYROP' , 'MYRP'  , 'NC4'   ,
            'NC5'   , 'NDEC'  , 'NER'   , 'NERP'  , 'NGLIPA', 'NGLIPB', 'NGLIPC',
            'NHEX'  , 'NSM'   , 'OLE'   , 'OLEP'  , 'OSM'   , 'OYPE'  , 'PAL'   ,
            'PALIPA', 'PALIPB', 'PALIPC', 'PALIPD', 'PALIPE', 'PALO'  , 'PALOP' ,
            'PALP'  , 'PC'    , 'PDAG'  , 'PDOPC' , 'PDOPE' , 'PENT'  , 'PGHG'  ,
            'PLPA'  , 'PLPC'  , 'PLPE'  , 'PLPG'  , 'PLPI'  , 'PLPI13', 'PLPI14',
            'PLPI15', 'PLPI24', 'PLPI25', 'PLPI2A', 'PLPI2B', 'PLPI2C', 'PLPI2D',
            'PLPI33', 'PLPI34', 'PLPI35', 'PLPS'  , 'PMCL1' , 'PMCL2' , 'PMPE'  ,
            'PMPG'  , 'PMSPE' , 'PMSPG' , 'PNPI'  , 'PNPI13', 'PNPI14', 'PNPI15',
            'PNPI24', 'PNPI25', 'PNPI2A', 'PNPI2B', 'PNPI2C', 'PNPI2D', 'PNPI33',
            'PNPI34', 'PNPI35', 'POPA'  , 'POPC'  , 'POPCE' , 'POPE'  , 'POPEE' ,
            'POPG'  , 'POPI'  , 'POPI13', 'POPI14', 'POPI15', 'POPI24', 'POPI25',
            'POPI2A', 'POPI2B', 'POPI2C', 'POPI2D', 'POPI33', 'POPI34', 'POPI35',
            'POPP1' , 'POPP2' , 'POPP3' , 'POPS'  , 'PPPE'  , 'PRPE'  , 'PSM'   ,
            'PVCL2' , 'PVPE'  , 'PVPG'  , 'PYPE'  , 'PYPG'  , 'PYPI'  , 'PhPC'  ,
            'QMPE'  , 'SAPA'  , 'SAPC'  , 'SAPE'  , 'SAPG'  , 'SAPI'  , 'SAPI13',
            'SAPI14', 'SAPI15', 'SAPI24', 'SAPI25', 'SAPI2A', 'SAPI2B', 'SAPI2C',
            'SAPI2D', 'SAPI33', 'SAPI34', 'SAPI35', 'SAPS'  , 'SB3_10', 'SB3_12',
            'SB3_14', 'SDA'   , 'SDAP'  , 'SDPA'  , 'SDPC'  , 'SDPE'  , 'SDPG'  ,
            'SDPS'  , 'SDS'   , 'SELIPA', 'SELIPB', 'SELIPC', 'SFLIPA', 'SITO'  ,
            'SLPA'  , 'SLPC'  , 'SLPE'  , 'SLPG'  , 'SLPS'  , 'SOPA'  , 'SOPC'  ,
            'SOPE'  , 'SOPG'  , 'SOPS'  , 'SSM'   , 'STE'   , 'STEP'  , 'STIG'  ,
            'TAG'   , 'TEA'   , 'TETD'  , 'THA'   , 'THAP'  , 'THCHL' , 'THDPPC',
            'TIPA'  , 'TLCL1' , 'TLCL2' , 'TMCL1' , 'TMCL2' , 'TOCL1' , 'TOCL2' ,
            'TPA'   , 'TPAP'  , 'TPG'   , 'TPT'   , 'TPTP'  , 'TRI'   , 'TRIP'  ,
            'TRIPAO', 'TRPAOP', 'TSPC'  , 'TTA'   , 'TTAP'  , 'TYCL1' , 'TYCL2' ,
            'UDAO'  , 'UDAOP' , 'UFOS10', 'VCLIPA', 'VCLIPB', 'VCLIPC', 'VCLIPD',
            'VCLIPE', 'YOPA'  , 'YOPC'  , 'YOPE'  , 'YOPG'  , 'YOPS'  , 'YPLIPA',
            'YPLIPB']},
        'water': {
            'macro': False,
            'resnames': ['HOH', 'TIP3', 'WAT']},
        'other': {
            'macro': False,
            'resnames': []}
    }

_atom_aliases = [
    "ILE CD1 CD",
    "BGLCNA C7 C",
    "BGLCNA O7 O",
    "BGLCNA C8 CT",
    "BGLCNA N2 N",
    "ANE5 C10 C",
    "ANE5 C11 CT",
    "ANE5 N5 N",
    "ANE5 O1A O11",
    "ANE5 O1B O12",
    "ANE5 O10 O",
    "VCG C01 C1",
    "VCG C01 C1",
    "VCG C02 C2",
    "VCG C03 C3",
    "VCG C04 C4",
    "VCG C05 C5",
    "VCG C06 C6",
    "VCG C07 C7",
    "VCG C08 C8",
    "VCG C09 C9",
    "TIP3 O OH2"
]
_residue_aliases = [
    "HIS HSD",
    "PO4 H2PO4",
    "H2PO H2PO4",
    "MAN AMAN",
    "BMA BMAN",
    "NAG BGLCNA",
    "NDG AGLCNA",
    "BGLC BGLCNA",
    "FUC AFUC",
    "FUL BFUC",
    "GAL BGAL",
    "ANE5 ANE5AC",
    "SIA ANE5AC",
    "EIC LIN",
    "HOH TIP3",
    "ZN ZN2",
    "CL CLA",
    "C6DH C6DHPC",
    "C7DH C7DHPC",
    "DT THY",
    "DA ADE",
    "DC CYT",
    "DG GUA",
    "DU URA",
    "HEM HEME"
]

class LabelMappers:
    """
    Class to hold label mappers for residue names and segment types.
    """
    def __init__(self):
        # process _data to initialize the label mappers
        self.aliases={}
        self.aliases['atom'] = _atom_aliases
        self.aliases['residue'] = _residue_aliases
        self.segtypes=segtypes
        self.segtype_of_resname = {}
        self.charmm_resname_of_pdb_resname = {}
        self.pdb_resname_of_charmm_resname = {}
        self.res_321 = segtypes['protein']['rescodes']
        self.res_123 = segtypes['protein']['invrescodes']
        for segtype in segtypes:
            for resname in segtypes[segtype]['resnames']:
                self.segtype_of_resname[resname] = segtype
        for alias in _residue_aliases:
            parts = alias.split()
            resname, alias1 = parts
            self.charmm_resname_of_pdb_resname[resname] = alias1
            self.pdb_resname_of_charmm_resname[alias1] = resname
    
    def update_atomselect_macros(self,fp):
        """
        Update the atomselect macros in the file ``fp`` based on the ``segtypes`` dict.
        This is a developer-only feature.  Access to this method is provided by the ``pestifer modify-package`` command (see :func:`pestifer.core.pestifer.modify_package`).
        """
        for segtype, data in self.segtypes.items():
            if data['macro']:
                resnames = data['resnames']
                macro_content = ' '.join(resnames)
                fp.write(f"update_atomselect_macro {segtype} \"resname {macro_content}\" 0\n")

Labels=LabelMappers()
""" 
Global instance of :class:`LabelMappers` class to access segment types and residue names.
This instance provides access to the segment types and residue names used in CHARMM and PDB files.
It allows for easy mapping between residue names and their corresponding segment types,
as well as providing access to the CHARMM residue names for PDB residue names.  Any module file
that imports ``Labels`` from ``pestifer.core.labels`` will have access to this instance.
This is useful for tasks that require residue name and segment type management, such as
preparing input files for molecular simulations or analyzing protein structures.
"""