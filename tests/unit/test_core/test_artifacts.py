import os
import unittest
from pestifer.core.artifacts import *
from pathlib import Path

class TestArtifacts(unittest.TestCase):

    def test_artifact_is_abstract(self):
        with self.assertRaises(TypeError):
            Artifact()
    
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
        self.assertEqual(data.value, 'This is some test data.')
        self.assertIsInstance(data, ArtifactData)

    def test_artifactfilelist_is_abstract(self):
        with self.assertRaises(TypeError):
            ArtifactFileList()

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