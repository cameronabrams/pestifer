# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import pytest
import os
import shutil
import unittest

from unittest.mock import patch

from pestifer.core.artifacts import ArtifactDict
from pestifer.core.config import Config
# from pestifer.charmmff.charmmffresidatabase import CHARMMFFResiDatabase
from pestifer.molecule.bilayer import Bilayer, specstrings_builddict

class TestBilayer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Ensure the resource manager is initialized before any tests run
        cls.C = Config().configure_new()
        cls.RM = cls.C.RM
        cls.charmmff_content = cls.RM.charmmff_content
        cls.pdbrepository = cls.charmmff_content.pdbrepository

    @classmethod
    def tearDownClass(cls):
        del cls.charmmff_content

    def test_bilayer_init_empty(self):
        self.charmmff_content.deprovision()
        with patch("pestifer.molecule.bilayer.logger.debug") as mock_logger:
            test_bilayer=Bilayer(composition_dict={})
            mock_logger.assert_called_once_with('Empty bilayer')
            self.assertEqual(type(test_bilayer), Bilayer)
            self.assertEqual(test_bilayer.artifacts, ArtifactDict())

    def test_bilayer_init_nonempty(self):
        self.charmmff_content.deprovision()
        if os.path.exists('__test_bilayer_init_nonempty'):
            shutil.rmtree('__test_bilayer_init_nonempty')
        os.mkdir('__test_bilayer_init_nonempty')
        os.chdir('__test_bilayer_init_nonempty')
        self.assertFalse(self.charmmff_content.provisioned)
        test_bilayer =Bilayer(composition_dict={
            'upper_leaflet': [{'name':'POPC','frac':1.0,'conf':0}],
            'lower_leaflet': [{'name':'POPE','frac':1.0,'conf':0}]},
            charmmffcontent=self.charmmff_content)
        assert all([x in test_bilayer.lipid_names for x in ['POPE', 'POPC']])
        assert all([x in ['POPE', 'POPC'] for x in test_bilayer.lipid_names])
        os.chdir('..')

    def test_bilayer_memgen_to_composition_simple_symm(self):
        self.charmmff_content.deprovision()
        cdict=specstrings_builddict(lipid_specstring='POPC',lipid_ratio_specstring='1',lipid_conformers_specstring='1')
        assert len(cdict['upper_leaflet']) == 1
        assert len(cdict['lower_leaflet']) == 1
        assert len(cdict['upper_chamber']) == 1
        assert len(cdict['lower_chamber']) == 1
        assert cdict['upper_leaflet'][0]['name'] == 'POPC'
        assert cdict['upper_leaflet'][0]['frac'] == 1.0
        assert cdict['upper_leaflet'][0]['conf'] == 1
        assert cdict['lower_leaflet'][0]['name'] == 'POPC'
        assert cdict['lower_leaflet'][0]['frac'] == 1.0
        assert cdict['lower_leaflet'][0]['conf'] == 1
        assert 'patn' not in cdict['upper_leaflet'][0]
        assert 'patn' not in cdict['lower_leaflet'][0]
        assert 'charge' not in cdict['upper_leaflet'][0]
        assert 'charge' not in cdict['lower_leaflet'][0]
        assert 'MW' not in cdict['upper_leaflet'][0]
        assert 'MW' not in cdict['lower_leaflet'][0]
        assert len(cdict['upper_chamber']) == 1
        assert len(cdict['lower_chamber']) == 1
        assert cdict['upper_chamber'][0]['name'] == 'TIP3'
        assert cdict['upper_chamber'][0]['frac'] == 1.0
        assert cdict['lower_chamber'][0]['name'] == 'TIP3'
        assert cdict['upper_chamber'][0]['frac'] == 1.0
        assert 'conf' not in cdict['upper_chamber'][0]
        assert 'conf' not in cdict['lower_chamber'][0]
        assert 'patn' not in cdict['upper_chamber'][0]
        assert 'patn' not in cdict['lower_chamber'][0]
        assert 'charge' not in cdict['upper_chamber'][0]
        assert 'charge' not in cdict['lower_chamber'][0]
        assert 'MW' not in cdict['upper_chamber'][0]
        assert 'MW' not in cdict['lower_chamber'][0]

    def test_bilayer_memgen_to_composition_simple_asymm(self):
        self.charmmff_content.deprovision()
        cdict=specstrings_builddict(lipid_specstring='POPC//POPE',lipid_ratio_specstring='1',lipid_conformers_specstring='1')
        assert len(cdict['upper_leaflet']) == 1
        assert len(cdict['lower_leaflet']) == 1
        assert len(cdict['upper_chamber']) == 1
        assert len(cdict['lower_chamber']) == 1
        assert cdict['upper_leaflet'][0]['name'] == 'POPC'
        assert cdict['upper_leaflet'][0]['frac'] == 1.0
        assert cdict['upper_leaflet'][0]['conf'] == 1
        assert cdict['lower_leaflet'][0]['name'] == 'POPE'
        assert cdict['lower_leaflet'][0]['frac'] == 1.0
        assert cdict['lower_leaflet'][0]['conf'] == 1
        assert 'patn' not in cdict['upper_leaflet'][0]
        assert 'patn' not in cdict['lower_leaflet'][0]
        assert 'charge' not in cdict['upper_leaflet'][0]
        assert 'charge' not in cdict['lower_leaflet'][0]
        assert 'MW' not in cdict['upper_leaflet'][0]
        assert 'MW' not in cdict['lower_leaflet'][0]
        assert len(cdict['upper_chamber']) == 1
        assert len(cdict['lower_chamber']) == 1
        assert cdict['upper_chamber'][0]['name'] == 'TIP3'
        assert cdict['upper_chamber'][0]['frac'] == 1.0
        assert cdict['lower_chamber'][0]['name'] == 'TIP3'
        assert cdict['upper_chamber'][0]['frac'] == 1.0
        assert 'conf' not in cdict['upper_chamber'][0]
        assert 'conf' not in cdict['lower_chamber'][0]
        assert 'patn' not in cdict['upper_chamber'][0]
        assert 'patn' not in cdict['lower_chamber'][0]
        assert 'charge' not in cdict['upper_chamber'][0]
        assert 'charge' not in cdict['lower_chamber'][0]
        assert 'MW ' not in cdict['upper_chamber'][0]
        assert 'MW' not in cdict['lower_chamber'][0]

    def test_bilayer_memgen_to_composition_complex_symm(self):
        self.charmmff_content.deprovision()
        cdict=specstrings_builddict(lipid_specstring='POPC:POPE',lipid_ratio_specstring='0.25:0.75',lipid_conformers_specstring='1')
        assert len(cdict['upper_leaflet']) == 2
        assert len(cdict['lower_leaflet']) == 2
        assert cdict['upper_leaflet'][0]['name'] == 'POPC'
        assert cdict['upper_leaflet'][0]['frac'] == 0.25
        assert cdict['upper_leaflet'][0]['conf'] == 1
        assert cdict['upper_leaflet'][1]['name'] == 'POPE'
        assert cdict['upper_leaflet'][1]['frac'] == 0.75
        assert cdict['upper_leaflet'][1]['conf'] == 1
        assert cdict['lower_leaflet'][0]['name'] == 'POPC'
        assert cdict['lower_leaflet'][0]['frac'] == 0.25
        assert cdict['lower_leaflet'][0]['conf'] == 1
        assert cdict['lower_leaflet'][1]['name'] == 'POPE'
        assert cdict['lower_leaflet'][1]['frac'] == 0.75
        assert cdict['lower_leaflet'][1]['conf'] == 1
        assert 'patn' not in cdict['upper_leaflet'][0]
        assert 'patn' not in cdict['upper_leaflet'][1]
        assert 'patn' not in cdict['lower_leaflet'][0]
        assert 'patn' not in cdict['lower_leaflet'][1]
        assert 'charge' not in cdict['upper_leaflet'][0]
        assert 'charge' not in cdict['upper_leaflet'][1]
        assert 'charge' not in cdict['lower_leaflet'][0]
        assert 'charge' not in cdict['lower_leaflet'][1]
        assert 'MW' not in cdict['upper_leaflet'][0]
        assert 'MW' not in cdict['upper_leaflet'][1]
        assert 'MW' not in cdict['lower_leaflet'][0]
        assert 'MW' not in cdict['lower_leaflet'][1]
        assert len(cdict['upper_chamber']) == 1
        assert len(cdict['lower_chamber']) == 1
        assert cdict['upper_chamber'][0]['name'] == 'TIP3'
        assert cdict['upper_chamber'][0]['frac'] == 1.0
        assert cdict['lower_chamber'][0]['name'] == 'TIP3'
        assert cdict['upper_chamber'][0]['frac'] == 1.0
        assert 'conf' not in cdict['upper_chamber'][0]
        assert 'conf' not in cdict['lower_chamber'][0]
        assert 'patn' not in cdict['upper_chamber'][0]
        assert 'patn' not in cdict['lower_chamber'][0]
        assert 'charge' not in cdict['upper_chamber'][0]
        assert 'charge' not in cdict['lower_chamber'][0]
        assert 'MW' not in cdict['upper_chamber'][0]
        assert 'MW' not in cdict['lower_chamber'][0]

    def test_bilayer_memgen_to_composition_complex_asymm(self):
        cdict=specstrings_builddict(lipid_specstring='POPC:POPE//POPC:CHL1',lipid_ratio_specstring='0.25:0.75//0.50:0.50',lipid_conformers_specstring='1')
        self.charmmff_content.deprovision()
        assert len(cdict['upper_leaflet']) == 2
        assert len(cdict['lower_leaflet']) == 2
        assert cdict['upper_leaflet'][0]['name'] == 'POPC'
        assert cdict['upper_leaflet'][0]['frac'] == 0.25
        assert cdict['upper_leaflet'][0]['conf'] == 1
        assert cdict['upper_leaflet'][1]['name'] == 'POPE'
        assert cdict['upper_leaflet'][1]['frac'] == 0.75
        assert cdict['upper_leaflet'][1]['conf'] == 1
        assert cdict['lower_leaflet'][0]['name'] == 'POPC'
        assert cdict['lower_leaflet'][0]['frac'] == 0.50
        assert cdict['lower_leaflet'][0]['conf'] == 1
        assert cdict['lower_leaflet'][1]['name'] == 'CHL1'
        assert cdict['lower_leaflet'][1]['frac'] == 0.50
        assert cdict['lower_leaflet'][1]['conf'] == 1
        assert 'patn' not in cdict['upper_leaflet'][0]
        assert 'patn' not in cdict['upper_leaflet'][1]
        assert 'patn' not in cdict['lower_leaflet'][0]
        assert 'patn' not in cdict['lower_leaflet'][1]
        assert 'charge' not in cdict['upper_leaflet'][0]
        assert 'charge' not in cdict['upper_leaflet'][1]
        assert 'charge' not in cdict['lower_leaflet'][0]
        assert 'charge' not in cdict['lower_leaflet'][1]
        assert 'MW' not in cdict['upper_leaflet'][0]
        assert 'MW' not in cdict['upper_leaflet'][1]
        assert 'MW' not in cdict['lower_leaflet'][0]
        assert 'MW' not in cdict['lower_leaflet'][1]
        assert len(cdict['upper_chamber']) == 1
        assert len(cdict['lower_chamber']) == 1
        assert cdict['upper_chamber'][0]['name'] == 'TIP3'
        assert cdict['upper_chamber'][0]['frac'] == 1.0
        assert cdict['lower_chamber'][0]['name'] == 'TIP3'
        assert cdict['upper_chamber'][0]['frac'] == 1.0
        assert 'conf' not in cdict['upper_chamber'][0]
        assert 'conf' not in cdict['lower_chamber'][0]
        assert 'patn' not in cdict['upper_chamber'][0]
        assert 'patn' not in cdict['lower_chamber'][0]
        assert 'charge' not in cdict['upper_chamber'][0]
        assert 'charge' not in cdict['lower_chamber'][0]
        assert 'MW' not in cdict['upper_chamber'][0]
        assert 'MW' not in cdict['lower_chamber'][0]

    def test_bilayer_init_memgen_style(self):
        self.charmmff_content.deprovision()
        if os.path.exists('__test_bilayer_init_memgen_style'):
            shutil.rmtree('__test_bilayer_init_memgen_style')
        os.mkdir('__test_bilayer_init_memgen_style')
        os.chdir('__test_bilayer_init_memgen_style')
        cdict=specstrings_builddict(lipid_specstring='POPC',lipid_ratio_specstring='1.0',lipid_conformers_specstring='1')
        test_bilayer=Bilayer(composition_dict=cdict,charmmffcontent=self.charmmff_content)
        assert test_bilayer.lipid_names == ['POPC']
        assert all([x in test_bilayer.species_names for x in ['POPC', 'POT', 'TIP3', 'CLA']])
        cdict=specstrings_builddict(lipid_specstring='POPC//POPE',lipid_ratio_specstring='1',lipid_conformers_specstring='1')
        self.charmmff_content.deprovision()
        test_bilayer=Bilayer(composition_dict=cdict,charmmffcontent=self.charmmff_content)

        assert all([x in test_bilayer.lipid_names for x in ['POPE', 'POPC']])
        assert all([x in ['POPE', 'POPC'] for x in test_bilayer.lipid_names])
        assert all([x in test_bilayer.species_names for x in ['POPE', 'POPC', 'TIP3', 'POT', 'CLA']])
        assert all([x in ['POPE', 'POPC', 'TIP3', 'POT', 'CLA'] for x in test_bilayer.species_names])
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['name'] == 'POPC'
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['name'] == 'POPE'
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['frac'] == 1.0
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['frac'] == 1.0
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['conf'] == 1
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['conf'] == 1
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['patn'] == 100
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['patn'] == 100
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['charge'] == 0.0
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['charge'] == 0.0
        assert test_bilayer.slices['upper_chamber']['composition'][0]['name'] == 'TIP3'
        assert test_bilayer.slices['upper_chamber']['composition'][0]['frac'] == 1.0
        assert test_bilayer.slices['upper_chamber']['composition'][0]['charge'] == 0.0
        assert test_bilayer.slices['upper_chamber']['composition'][0]['patn'] == 3200
        assert test_bilayer.slices['upper_chamber']['composition'][0]['MW'] == 18.0154
        assert test_bilayer.slices['lower_chamber']['composition'][0]['name'] == 'TIP3'
        assert test_bilayer.slices['lower_chamber']['composition'][0]['frac'] == 1.0
        assert test_bilayer.slices['lower_chamber']['composition'][0]['charge'] == 0.0
        assert test_bilayer.slices['lower_chamber']['composition'][0]['patn'] == 3200
        assert test_bilayer.slices['lower_chamber']['composition'][0]['MW'] == 18.0154
        assert test_bilayer.asymmetric == True

        cdict=specstrings_builddict(lipid_specstring='POPC:CHL1//POPE:CHL1',lipid_ratio_specstring='0.5:0.5',lipid_conformers_specstring='1:1')
        self.charmmff_content.deprovision()
        test_bilayer=Bilayer(composition_dict=cdict,charmmffcontent=self.charmmff_content)
        assert all([x in test_bilayer.lipid_names for x in ['POPE', 'POPC','CHL1']])
        assert all([x in ['POPE', 'POPC','CHL1', 'POT', 'CLA'] for x in test_bilayer.lipid_names])
        assert all([x in test_bilayer.species_names for x in ['POPE', 'POPC','CHL1','TIP3', 'POT', 'CLA']])
        assert all([x in ['POPE', 'POPC','CHL1','TIP3', 'POT', 'CLA'] for x in test_bilayer.species_names])
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['name'] == 'POPC'
        assert test_bilayer.slices['upper_leaflet']['composition'][1]['name'] == 'CHL1'
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['name'] == 'POPE'
        assert test_bilayer.slices['lower_leaflet']['composition'][1]['name'] == 'CHL1'
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['frac'] == 0.5
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['frac'] == 0.5
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['conf'] == 1
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['conf'] == 1
        assert test_bilayer.slices['upper_leaflet']['composition'][1]['frac'] == 0.5
        assert test_bilayer.slices['lower_leaflet']['composition'][1]['frac'] == 0.5
        assert test_bilayer.slices['upper_leaflet']['composition'][1]['conf'] == 1
        assert test_bilayer.slices['lower_leaflet']['composition'][1]['conf'] == 1
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['patn'] == 50
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['patn'] == 50
        assert test_bilayer.slices['upper_leaflet']['composition'][1]['patn'] == 50
        assert test_bilayer.slices['lower_leaflet']['composition'][1]['patn'] == 50
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['charge'] == 0.0
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['charge'] == 0.0
        assert test_bilayer.slices['upper_leaflet']['composition'][1]['charge'] == 0.0
        assert test_bilayer.slices['lower_leaflet']['composition'][1]['charge'] == 0.0
        assert test_bilayer.slices['upper_chamber']['composition'][0]['name'] == 'TIP3'
        assert test_bilayer.slices['upper_chamber']['composition'][0]['frac'] == 1.0
        assert test_bilayer.slices['upper_chamber']['composition'][0]['charge'] == 0.0
        assert test_bilayer.slices['upper_chamber']['composition'][0]['patn'] == 3200
        assert test_bilayer.slices['upper_chamber']['composition'][0]['MW'] == 18.0154
        assert test_bilayer.slices['lower_chamber']['composition'][0]['name'] == 'TIP3'
        assert test_bilayer.slices['lower_chamber']['composition'][0]['frac'] == 1.0
        assert test_bilayer.slices['lower_chamber']['composition'][0]['charge'] == 0.0
        assert test_bilayer.slices['lower_chamber']['composition'][0]['patn'] == 3200
        assert test_bilayer.slices['lower_chamber']['composition'][0]['MW'] == 18.0154
        assert test_bilayer.asymmetric == True

        cdict=specstrings_builddict(lipid_specstring='POPS:CHL1//POPE:CHL1',lipid_ratio_specstring='0.75:0.25',lipid_conformers_specstring='3:4//7:2')
        self.charmmff_content.deprovision()
        test_bilayer=Bilayer(composition_dict=cdict,charmmffcontent=self.charmmff_content)
        assert all([x in test_bilayer.lipid_names for x in ['POPS', 'POPE','CHL1']])
        assert all([x in ['POPS', 'POPE','CHL1'] for x in test_bilayer.lipid_names])
        assert all([x in test_bilayer.species_names for x in ['POPS', 'POPE','CHL1','TIP3', 'POT', 'CLA']])
        assert all([x in ['POPS', 'POPE','CHL1','TIP3','POT', 'CLA'] for x in test_bilayer.species_names])
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['name'] == 'POPS'
        assert test_bilayer.slices['upper_leaflet']['composition'][1]['name'] == 'CHL1'
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['name'] == 'POPE'
        assert test_bilayer.slices['lower_leaflet']['composition'][1]['name'] == 'CHL1'
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['frac'] == 0.75
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['frac'] == 0.75
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['conf'] == 3
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['conf'] == 7
        assert test_bilayer.slices['upper_leaflet']['composition'][1]['frac'] == 0.25
        assert test_bilayer.slices['lower_leaflet']['composition'][1]['frac'] == 0.25
        assert test_bilayer.slices['upper_leaflet']['composition'][1]['conf'] == 4
        assert test_bilayer.slices['lower_leaflet']['composition'][1]['conf'] == 2
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['patn'] == 75
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['patn'] == 75
        assert test_bilayer.slices['upper_leaflet']['composition'][1]['patn'] == 25
        assert test_bilayer.slices['lower_leaflet']['composition'][1]['patn'] == 25
        assert test_bilayer.slices['upper_leaflet']['composition'][0]['charge'] == -1.0
        assert test_bilayer.slices['lower_leaflet']['composition'][0]['charge'] == 0.0
        assert test_bilayer.slices['upper_leaflet']['composition'][1]['charge'] == 0.0
        assert test_bilayer.slices['lower_leaflet']['composition'][1]['charge'] == 0.0
        assert test_bilayer.slices['upper_chamber']['composition'][0]['name'] == 'TIP3'
        assert test_bilayer.slices['upper_chamber']['composition'][0]['frac'] == 1.0
        assert test_bilayer.slices['upper_chamber']['composition'][0]['charge'] == 0.0
        assert test_bilayer.slices['upper_chamber']['composition'][0]['patn'] == 3200
        assert test_bilayer.slices['upper_chamber']['composition'][0]['MW'] == 18.0154
        assert test_bilayer.slices['lower_chamber']['composition'][0]['name'] == 'TIP3'
        assert test_bilayer.slices['lower_chamber']['composition'][0]['frac'] == 1.0
        assert test_bilayer.slices['lower_chamber']['composition'][0]['charge'] == 0.0
        assert test_bilayer.slices['lower_chamber']['composition'][0]['patn'] == 3200
        assert test_bilayer.slices['lower_chamber']['composition'][0]['MW'] == 18.0154
        assert test_bilayer.slices['upper_chamber']['composition'][1]['name'] == 'POT'
        assert test_bilayer.slices['upper_chamber']['composition'][1]['charge'] == 1.0
        assert test_bilayer.slices['upper_chamber']['composition'][1]['patn'] == 38
        assert test_bilayer.slices['upper_chamber']['composition'][1]['MW'] == 39.0983
        assert test_bilayer.slices['lower_chamber']['composition'][1]['name'] == 'POT'
        assert test_bilayer.slices['lower_chamber']['composition'][1]['charge'] == 1.0
        assert test_bilayer.slices['lower_chamber']['composition'][1]['patn'] == 37
        assert test_bilayer.slices['lower_chamber']['composition'][1]['MW'] == 39.0983
        assert test_bilayer.asymmetric == True
        self.RM.charmmff_content.clean_local_charmmff_files()
        os.chdir('..')

    def test_bilayer_spec_out(self):
        self.charmmff_content.deprovision()
        if os.path.exists('__test_bilayer_spec_out'):
            shutil.rmtree('__test_bilayer_spec_out')
        os.mkdir('__test_bilayer_spec_out')
        os.chdir('__test_bilayer_spec_out')
        cdict=specstrings_builddict(lipid_specstring='POPC:CHL1//POPE:CHL1', lipid_ratio_specstring='0.75:0.25//0.33:0.67', lipid_conformers_specstring='3:4//7:2')
        test_bilayer = Bilayer(composition_dict=cdict, charmmffcontent=self.charmmff_content)
        test_bilayer.spec_out()
        assert test_bilayer.patch_ll_corner[0]==pytest.approx(0.0, rel=1e-2)
        assert test_bilayer.patch_ll_corner[1]==pytest.approx(0.0, rel=1e-2)
        assert test_bilayer.patch_ll_corner[2]==pytest.approx(0.0, rel=1e-2)
        assert test_bilayer.patch_ur_corner[0]==pytest.approx(86.60, rel=1e-2)
        assert test_bilayer.patch_ur_corner[1]==pytest.approx(86.60, rel=1e-2)
        assert test_bilayer.patch_ur_corner[2]==pytest.approx(83.08, rel=1e-2)
        self.RM.charmmff_content.clean_local_charmmff_files()

        os.chdir('..')
        
