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

    def test_bilayer_init_with_salt_gives_every_species_a_conf_key(self):
        # regression: with salt_concentration>0 the constructor appends ion species to the chambers;
        # the phased-conformer routing keys every grid-placement checkout on 'conf_key', so a late
        # ion append that omitted it used to KeyError (surfaced building ex17, which is salted).
        self.charmmff_content.deprovision()
        d = '__test_bilayer_salt'
        if os.path.exists(d):
            shutil.rmtree(d)
        os.mkdir(d)
        os.chdir(d)
        try:
            b = Bilayer(composition_dict={
                'upper_leaflet': [{'name': 'POPC', 'frac': 1.0, 'conf': 0}],
                'lower_leaflet': [{'name': 'POPE', 'frac': 1.0, 'conf': 0}]},
                salt_concentration=0.154,
                charmmffcontent=self.charmmff_content)
            for layer in b.slices.values():
                for sp in layer['composition']:
                    assert 'conf_key' in sp, f'species {sp.get("name")} missing conf_key'
        finally:
            os.chdir('..')

    def test_constructor_does_not_mutate_composition_dict(self):
        # Regression: the constructor fills per-species 'patn' from leaflet_nlipids, but must
        # do so on its own copy. If it stamps 'patn' into the caller's dict, a small patch's
        # counts leak into a larger membrane built from the same dict -- under-populating and
        # under-hydrating it (the dehydrated bilayer then collapses laterally).
        self.charmmff_content.deprovision()
        cdict = specstrings_builddict(lipid_specstring='POPC', lipid_ratio_specstring='1',
                                      lipid_conformers_specstring='1')
        Bilayer(composition_dict=cdict, leaflet_nlipids={'upper': 100, 'lower': 100},
                charmmffcontent=self.charmmff_content)
        for key in ('upper_leaflet', 'lower_leaflet', 'upper_chamber', 'lower_chamber'):
            for entry in cdict[key]:
                self.assertNotIn('patn', entry, f'{key} entry was mutated with patn')
                self.assertNotIn('charge', entry, f'{key} entry was mutated with charge')

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
        # Orthohexagonal lattice: the box is no longer square -- its dimensions follow from the
        # integer nx x ny cell grid (here 10 x 10, zero vacancies for 100 lipids) with row pitch
        # sqrt(3)/2 of the column pitch, so Lx != Ly while patch_area (= n*SAPL = 7500) is preserved.
        assert test_bilayer.patch_ur_corner[0]==pytest.approx(93.06, rel=1e-2)
        assert test_bilayer.patch_ur_corner[1]==pytest.approx(80.59, rel=1e-2)
        assert test_bilayer.patch_area==pytest.approx(7500.0, rel=1e-2)
        # box z reserves each leaflet's true coordinate z-extent (choline + caps), not just the
        # head-to-tail-tip length -- so lipids don't overflow into and thin the water chambers.
        # (78.14 with the athermal-MC conformer set, whose melted tails give a smaller z-extent than
        # the old extended-rod conformers that gave 87.86.)
        assert test_bilayer.patch_ur_corner[2]==pytest.approx(78.14, rel=1e-2)
        self.RM.charmmff_content.clean_local_charmmff_files()

        os.chdir('..')

    def test_bilayer_write_grid_pdb(self):
        self.charmmff_content.deprovision()
        d = '__test_bilayer_write_grid_pdb'
        if os.path.exists(d):
            shutil.rmtree(d)
        os.mkdir(d)
        os.chdir(d)
        cdict = specstrings_builddict(lipid_specstring='POPC:CHL1//POPC:CHL1',
                                      lipid_ratio_specstring='0.7:0.3//0.7:0.3',
                                      lipid_conformers_specstring='0:0')
        b = Bilayer(composition_dict=cdict, leaflet_nlipids=dict(upper=64, lower=64),
                    charmmffcontent=self.charmmff_content)
        b.spec_out(SAPL=55.0)
        out = b.write_grid_pdb('grid_patch.pdb', seed=1, clash_cutoff=1.2)
        assert os.path.exists(out)
        natoms = nwat = nlip_lower = nlip_upper = 0
        lip_xyz, sol_xyz = [], []
        with open(out) as fh:
            for ln in fh:
                if ln.startswith(('ATOM', 'HETATM')):
                    natoms += 1
                    xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
                    z = xyz[2]
                    rn = ln[17:21].strip()
                    if rn in ('TIP3', 'POT', 'CLA', 'SOD'):
                        nwat += 1
                        sol_xyz.append(xyz)
                    else:
                        lip_xyz.append(xyz)
                        if z < b.midplane_z:
                            nlip_lower += 1
                        else:
                            nlip_upper += 1
        # a complete patch: lipids in both leaflets and solvent in the chambers
        assert natoms > 0
        assert nwat > 0
        assert nlip_lower > 0 and nlip_upper > 0
        # no chamber-solvent atom may sit within the clash cutoff of a lipid atom
        # (the near-coincidence that broke VMD bond perception / psfgen coordinate read)
        import numpy as _np
        from scipy.spatial import cKDTree
        dmin = cKDTree(_np.array(lip_xyz)).query(_np.array(sol_xyz))[0].min()
        # tolerate a hair below the 1.2 A cutoff: the packer's own cKDTree removal and this independent
        # recompute take different float paths, so a pair right on the boundary can read 1.1995 here
        # while the removal saw >= 1.2. A real surviving clash is far smaller than this.
        assert dmin >= 1.2 - 1e-2, f"solvent-lipid clash survived: min distance {dmin:.3f} A"
        self.RM.charmmff_content.clean_local_charmmff_files()
        os.chdir('..')

    def test_bilayer_write_grid_pdb_samples_conformers(self):
        """Lever 1: with no pinned conformer the packer draws per-lipid across the whole
        conformer ensemble instead of stamping one frozen shape onto every lipid (the
        over-packing root cause).  Each placed lipid is fingerprinted by its maximum
        internal pairwise distance, which is invariant under the rigid placement transforms
        (mirror + in-plane spin + translation), so distinct fingerprints prove distinct
        source conformers were used."""
        import numpy as _np
        self.charmmff_content.deprovision()
        d = '__test_bilayer_grid_conformer_sampling'
        if os.path.exists(d):
            shutil.rmtree(d)
        os.mkdir(d)
        os.chdir(d)
        # single species, no conformer pin (empty conformers specstring) -> ensemble draw
        cdict = specstrings_builddict(lipid_specstring='POPC//POPC',
                                      lipid_ratio_specstring='1.0//1.0')
        b = Bilayer(composition_dict=cdict, leaflet_nlipids=dict(upper=48, lower=48),
                    charmmffcontent=self.charmmff_content)
        # the ensemble must actually be available for the draw to be meaningful
        assert len(b.species_data['POPC'].pdbcontents) > 1
        # and no species should carry a pinned conformer
        for layer in ('upper_leaflet', 'lower_leaflet'):
            for species in b.slices[layer]['composition']:
                assert 'conf' not in species
        b.spec_out(SAPL=55.0)
        out = b.write_grid_pdb('grid_patch.pdb', seed=1, clash_cutoff=1.2)
        residues = {}
        with open(out) as fh:
            for ln in fh:
                if ln.startswith(('ATOM', 'HETATM')):
                    if ln[17:21].strip() in ('TIP3', 'POT', 'CLA', 'SOD'):
                        continue
                    residues.setdefault(ln[22:26], []).append(
                        (float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
        assert len(residues) > 1
        fingerprints = set()
        for xyz in residues.values():
            a = _np.array(xyz)
            dmax = _np.sqrt(((a[:, None, :] - a[None, :, :]) ** 2).sum(-1)).max()
            fingerprints.add(round(float(dmax), 1))
        assert len(fingerprints) > 1, (
            f'grid placement used a single conformer shape ({fingerprints}); '
            f'expected a per-lipid draw across the {len(b.species_data["POPC"].pdbcontents)}'
            f'-conformer ensemble')
        self.RM.charmmff_content.clean_local_charmmff_files()
        os.chdir('..')
        
