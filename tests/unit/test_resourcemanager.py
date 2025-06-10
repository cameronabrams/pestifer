import unittest
import os
from pestifer.resourcemanager import ResourceManager
from pestifer.charmmffcontent import CHARMMFFResiDatabase, CHARMMFFContent
from pestifer.pdbrepository import PDBInput
from pestifer import resources

class TestResourceManager(unittest.TestCase):
    def test_resource_manager(self):
        RM=ResourceManager()
        self.assertEqual(RM.resources_path,os.path.dirname(resources.__file__))
        for r in RM.base_resources:
            self.assertTrue(os.path.exists(RM.resource_path[r]))
    
    def test_resource_get_example_yaml(self):
        RM=ResourceManager()
        c=RM.get_example_yaml_by_index(1)
        self.assertEqual(os.path.basename(c),'01-bpti.yaml')

    def test_resource_ycleptic(self):
        RM=ResourceManager()
        self.assertEqual(os.path.basename(RM.ycleptic_config),"base.yaml")
        self.assertEqual(os.path.basename(RM.get_ycleptic_config()),"base.yaml")

    def test_resource_tcl(self):
        RM=ResourceManager()
        self.assertTrue(os.path.exists(RM.get_tcldir()))
        self.assertTrue(os.path.exists(RM.get_tcl_pkgdir()))
        self.assertTrue(os.path.exists(RM.get_tcl_scriptsdir()))

    def test_charmmff_content(self):
        RM=ResourceManager()
        RM.charmmff_content.clean_local_charmmff_files()
        self.assertTrue(RM.charmmff_content!=None)
        self.assertTrue(RM.charmmff_content.filenamemap!={})
        basenames=[k for k in RM.charmmff_content.filenamemap.keys()]
        self.assertTrue(len(basenames)>0)
        self.assertEqual(len(set(basenames)),len(basenames))  # check for duplicates
        self.assertTrue(len(RM.charmmff_content.streams)>0)
        self.assertEqual(RM.charmmff_content.streams,['prot', 'cphmd', 'carb', 'na', 'lipid', 'misc'])
        self.assertTrue(RM.charmmff_content.custom_files!=None)

        RM.charmmff_content.copy_charmmfile_local('par_all36m_prot.prm')
        self.assertTrue(os.path.exists('par_all36m_prot.prm'))
        RM.charmmff_content.copy_charmmfile_local('top_all36_prot.rtf')
        self.assertTrue(os.path.exists('top_all36_prot.rtf'))
        RM.charmmff_content.copy_charmmfile_local('toppar_water_ions.str')
        self.assertTrue(os.path.exists('toppar_water_ions.str'))
        with open('toppar_water_ions.str','r') as f:
            contents=f.read()
            self.assertTrue('! commented out by pestifer' in contents)
        RM.charmmff_content.copy_charmmfile_local('toppar_all36_carb_glycopeptide.str')
        self.assertTrue(os.path.exists('toppar_all36_carb_glycopeptide.str'))
        with open('toppar_all36_carb_glycopeptide.str','r') as f:
            contents=f.read()
            self.assertFalse('! commented out by pestifer' in contents)
        RM.charmmff_content.copy_charmmfile_local('toppar_all36_moreions.str')
        self.assertTrue(os.path.exists('toppar_all36_moreions.str'))
        RM.charmmff_content.clean_local_charmmff_files()
        self.assertFalse(os.path.exists('par_all36m_prot.prm'))
        self.assertFalse(os.path.exists('top_all36_prot.rtf'))
        self.assertFalse(os.path.exists('toppar_water_ions.str'))
        self.assertFalse(os.path.exists('toppar_all36_carb_glycopeptide.str'))
        self.assertFalse(os.path.exists('toppar_all36_moreions.str'))

    def test_pdb_repository(self):
        RM=ResourceManager()
        PD=RM.charmmff_content.pdb_repository
        self.assertTrue(PD!=None)
        c=PD.checkout('PSM')
        self.assertTrue(c!=None)
        ch=c.get_charge()
        self.assertTrue(ch==0.0)
        p=c.get_parameters()
        self.assertTrue('toppar_all36_lipid_sphingo.str' in p)
        self.assertTrue(len(c.info['conformers'])==10)
        self.assertEqual(c.info['conformers'][0]['head-tail-length'],26.779)
        self.assertEqual(c.info['conformers'][0]['max-internal-length'],30.398)
        c.get_pdb(0)
        self.assertTrue(os.path.exists('PSM-00.pdb'))
        os.remove('PSM-00.pdb')
        c.get_pdb(0,noh=True)
        self.assertTrue(os.path.exists('PSM-00-noh.pdb'))
        os.remove('PSM-00-noh.pdb')
        c=PD.checkout('TIP3')
        self.assertTrue(c!=None)
        self.assertTrue(c.info=={})
        c.get_pdb(0)
        self.assertTrue(os.path.exists('TIP3.pdb'))
        c=PD.checkout('FAKE')
        self.assertTrue(c==None)
        PD.registercollection('pdb_depot','user')
        c=PD.checkout('FAKE')
        self.assertTrue(c!=None)
        c.get_pdb(0)
        self.assertTrue(os.path.exists('FAKE-00.pdb'))
        os.remove('FAKE-00.pdb')

    def test_charmmff_residatabase(self):
        RM=ResourceManager()
        C=RM.charmmff_content
        CDB=CHARMMFFResiDatabase(C)
        self.assertTrue(CDB!=None)
        self.assertEqual(len(CDB.residues),1101)
        self.assertTrue('ALA' in CDB.residues)
        self.assertTrue('TIP3' in CDB.residues)
        self.assertTrue('FAKE' not in CDB.residues)
        ala=CDB.residues['ALA']
        self.assertTrue(ala!=None)
        self.assertTrue(ala.metadat['stream']=='prot')
        self.assertTrue(ala.metadat['substream']=='')
        self.assertEqual(ala.mass(),89.093)
