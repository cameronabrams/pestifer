import os
import unittest
from pestifer.core.artifacts import *
from pathlib import Path


@dataclass
class Stamper:
    name: str



class TestArtifacts(unittest.TestCase):

    
    def test_artifactfile_is_abstract(self):
        with self.assertRaises(TypeError):
            ArtifactFile()

    def test_TXTFile(self):
        Path('test.txt').write_text('This is a test file.')
        txt_file = TXTFile('test')
        self.assertTrue(txt_file.exists())
        self.assertEqual(txt_file.path, Path('test.txt'))
        os.remove(txt_file.path)  # Clean up the test file
        self.assertFalse(txt_file.exists(), "TXTFile should be removed after test.")
        txt_file=TXTFile()
        self.assertEqual(txt_file.mime_type, 'text/plain')

    def test_artifactdata(self):
        data = ArtifactData('This is some test data.')
        self.assertEqual(data.key, None)
        self.assertEqual(data.data, 'This is some test data.')
        self.assertIsInstance(data, ArtifactData)
        
    def test_artifact_stamp(self):
        stamper = Stamper(name='TestTask')
        artifact = ArtifactData('Test data')
        self.assertFalse(artifact.has_stamp())
        artifact.stamp(stamper)
        self.assertTrue(artifact.has_stamp())
        self.assertEqual(artifact.produced_by, stamper)

class TestArtifactList(unittest.TestCase):
    
    def test_artifact_list(self):
        al = ArtifactList()
        self.assertEqual(len(al), 0)
        data1 = ArtifactData('Data 1')
        data2 = ArtifactData('Data 2')
        al.append(data1)
        al.append(data2)
        self.assertEqual(len(al), 2)
        self.assertIn(data1, al)
        self.assertIn(data2, al)

    def test_charmmff_files(self):
        Path('p1.prm').write_text('CHARMM parameter file content.')
        Path('p1.rtf').write_text('CHARMM topology file content.')
        Path('p1.str').write_text('CHARMM stream file content.')
        Path('p2.prm').write_text('CHARMM parameter file content.')
        Path('p2.rtf').write_text('CHARMM topology file content.')
        Path('p2.str').write_text('CHARMM stream file content.')
        Path('p3.prm').write_text('CHARMM parameter file content.')
        Path('p3.rtf').write_text('CHARMM topology file content.')
        Path('p3.str').write_text('CHARMM stream file content.')
        p=CharmmffParFiles([CharmmffParFile(x) for x in ['p1','p2']])
        self.assertEqual(len(p), 2)
        self.assertTrue(p.all_exist(), "All CHARMM parameter files should exist.")
        p.append(CharmmffParFile('p3'))
        self.assertEqual(len(p), 3)
        t=CharmmffTopFiles()
        t.append(CharmmffTopFile('p1'))
        self.assertEqual(len(t), 1)
        t.append(CharmmffTopFile('p2'))
        self.assertEqual(len(t), 2)
        t.append(CharmmffTopFile('p3'))
        self.assertEqual(len(t), 3)
        s=CharmmffStreamFiles()
        s.append(CharmmffStreamFile('p1'))
        self.assertEqual(len(s), 1)
        ss=[CharmmffStreamFile(x) for x in ['p2','p3']]
        s.extend(ss)
        self.assertEqual(len(s), 3)
        for L in [p,t,s]:
            for f in L:
                self.assertTrue(f.exists())
                os.remove(f.path)  # Clean up the test files
        self.assertFalse(any(f.exists() for f in p), "All CHARMM files should be removed after test.")
        self.assertFalse(any(f.exists() for f in t), "All CHARMM files should be removed after test.")
        self.assertFalse(any(f.exists() for f in s), "All CHARMM files should be removed after test.")

    def test_artifactlist_filter_by_produced_by(self):
        stamper1 = Stamper(name='Task1')
        stamper2 = Stamper(name='Task2')
        al = ArtifactList()
        al.append(ArtifactData('Data 1', produced_by=stamper1))
        al.append(ArtifactData('Data 2', produced_by=stamper2))
        al.append(ArtifactData('Data 3', produced_by=stamper1))
        
        filtered = al.filter_by_produced_by(stamper1)
        self.assertEqual(len(filtered), 2)
        self.assertTrue(all(a.produced_by == stamper1 for a in filtered))

class TestArtifactDict(unittest.TestCase):
    
    def test_artifact_dict(self):
        ad = ArtifactDict()
        self.assertEqual(len(ad), 0)
        data1 = ArtifactData('Data 1')
        data2 = ArtifactData('Data 2')
        ad['key1'] = data1
        ad['key2'] = data2
        self.assertEqual(len(ad), 2)
        self.assertIn('key1', ad)
        self.assertIn('key2', ad)
        self.assertEqual(ad['key1'], data1)
        self.assertEqual(ad['key2'], data2)

    def test_artifact_dict_filter_by_produced_by(self):
        stamper1 = Stamper(name='Task1')
        stamper2 = Stamper(name='Task2')
        ad = ArtifactDict()
        ad['key1'] = ArtifactData('Data 1', produced_by=stamper1)
        ad['key2'] = ArtifactData('Data 2', produced_by=stamper2)
        ad['key3'] = ArtifactData('Data 3', produced_by=stamper1)
        
        filtered = ad.filter_by_produced_by(stamper1)
        self.assertEqual(len(filtered), 2)
        self.assertTrue(all(a.produced_by == stamper1 for a in filtered.values()))
        self.assertIn('key1', filtered)
        self.assertIn('key3', filtered)
        self.assertNotIn('key2', filtered)

    def test_artifact_dict_update_item(self):
        ad = ArtifactDict()
        data1 = ArtifactData('Data 1')
        ad['key1'] = data1
        self.assertEqual(ad['key1'], data1)
        data2 = ArtifactData('Data 2')
        ad.update_item(data2, key='key2')
        self.assertEqual(ad['key2'], data2)

