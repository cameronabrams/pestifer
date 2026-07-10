import unittest
import shutil
import os
from pathlib import Path
from pestifer.core.examplemanager import ExampleManager

class TestExampleManager(unittest.TestCase):

    def setUp(self):
        assert os.getcwd().endswith('tests/unit/test_core/test_examplemanager')
        assert os.path.isdir('../fixtures')
        assert os.path.isdir('../fixtures/test_examplemanager_inputs')
        assert os.path.isdir('../fixtures/test_examplemanager_inputs/project')
        assert os.path.isdir('../fixtures/test_examplemanager_inputs/userspace')
        if os.path.isdir('./project'):
            shutil.rmtree('./project', ignore_errors=True)
        if os.path.isdir('./userspace'):
            shutil.rmtree('./userspace', ignore_errors=True)
        shutil.copytree('../fixtures/test_examplemanager_inputs/project', './project/', dirs_exist_ok=True)
        shutil.copytree('../fixtures/test_examplemanager_inputs/userspace', './userspace/', dirs_exist_ok=True)
        self.manager = ExampleManager(
            examples_path='project/package/resources/examples',
            docs_source_path='project/docs/source'
        )

    def tearDown(self):
        shutil.rmtree('project')
        shutil.rmtree('userspace')

    def test_example_manager_init(self):
        self.assertIsNotNone(self.manager)
        self.assertEqual(len(self.manager.examples), 0)
        self.assertEqual(self.manager.path, Path(os.path.join(os.getcwd(), 'project', 'package', 'resources', 'examples')))
        self.assertTrue(os.path.isdir(self.manager.path))

    def _build_example_set(self):
        os.chdir('userspace')
        self.example1 = self.manager.append_example(1, 'exA.yaml', auxiliary_inputs=['exA_companion1.pdb', 'exA_companion1.psf'])
        self.example2 = self.manager.append_example(2, 'exB.yaml')
        os.chdir('..')

    def test_example_manager_build_example_set(self):
        self._build_example_set()
        self.assertEqual(len(self.manager.examples), 2)
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, self.manager.examples[0].scriptpath)))
        # auxiliary companions live in <inputs>/aux/ (so the inputs dir holds only the main yaml)
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, self.manager.examples[0].auxpath, 'exA_companion1.pdb')))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, self.manager.examples[0].auxpath, 'exA_companion1.psf')))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, self.manager.examples[1].inputspath, 'exB.yaml')))
        self.assertTrue(os.path.isfile(self.manager.sphinx_example_manager.examples_rst))
        self.assertTrue(os.path.isdir(self.manager.sphinx_example_manager.examples_folder_path))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.sphinx_example_manager.examples_folder_path, '01', 'exA.rst')))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.sphinx_example_manager.examples_folder_path, '02', 'exB.rst')))
    
    def test_yaml_aux_example_is_registered_and_checked_out(self):
        # regression: an example with additional YAML helper scripts as auxiliary_inputs must
        # still be discovered (the scanner keys off the single main yaml in inputs/, with aux
        # scripts tucked in inputs/aux/) and fully retrieved by checkout.
        import tempfile
        os.chdir('userspace')
        self.manager.append_example(1, 'exA.yaml', auxiliary_inputs=['exB.yaml'])
        os.chdir('..')
        # the aux yaml lives under inputs/aux, not inputs, so inputs holds exactly one yaml
        ex = self.manager.examples[0]
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, ex.auxpath, 'exB.yaml')))
        # a fresh manager re-scanning the resources must still find the example
        from pestifer.core.examplemanager import ExampleManager
        m2 = ExampleManager(examples_path='project/package/resources/examples',
                            docs_source_path='project/docs/source')
        self.assertEqual(len(m2.examples), 1)
        self.assertEqual(m2.examples[0].shortname, 'exA')
        # checkout retrieves both the main script and the aux helper into the CWD
        out = tempfile.mkdtemp(); cwd = os.getcwd(); os.chdir(out)
        try:
            m2.checkout_example(1)
            self.assertTrue(os.path.isfile('exA.yaml'))
            self.assertTrue(os.path.isfile('exB.yaml'))
        finally:
            os.chdir(cwd); shutil.rmtree(out)

    def test_append_example_from_path_yields_bare_shortname(self):
        # append_example may be given a path outside the CWD (e.g. `example add /some/dir/x.yaml`);
        # the shortname must be the bare basename (no dir, no ext) or the toctree entry and the
        # resource folder key are corrupted (regression: a full path leaked into the toctree).
        script = os.path.abspath('userspace/exB.yaml')  # a path with a directory component
        ex = self.manager.append_example(3, script)
        self.assertEqual(ex.shortname, 'exB')
        # the yaml was copied into the resource folder, and the toctree entry is 'NN/exB'
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, ex.inputspath, 'exB.yaml')))
        toctree = open(self.manager.sphinx_example_manager.examples_rst).read()
        self.assertIn('03/exB', toctree)
        self.assertNotIn(os.path.abspath('userspace'), toctree)

    def test_example_manager_delete_example(self):
        self._build_example_set()
        self.manager.delete_example(2)
        self.assertEqual(len(self.manager.examples), 1)
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, self.manager.examples[0].scriptpath)))
        self.assertFalse(os.path.isfile(os.path.join(self.manager.path, self.example2.scriptpath)))
        self.assertFalse(os.path.isfile(os.path.join(self.manager.sphinx_example_manager.examples_folder_path, '02', 'exB.rst')))

    def test_example_manager_update_example_inplace(self):
        self._build_example_set()
        ex_rst = os.path.join(self.manager.sphinx_example_manager.examples_folder_path, '01', 'exA.rst')
        self.assertTrue(os.path.isfile(ex_rst))
        with open(ex_rst, 'r') as f:
            rst_current = f.read()
        self.assertNotIn('PDB ID 9ggs', rst_current)
        os.chdir('userspace')
        os.mkdir('exa_updated')
        os.chdir('exa_updated')
        Path('exA_companion1.pdb').touch()
        Path('exA_companion1.psf').touch()
        self.manager.update_example(1, title='Updated title', db_id='9ggs', author_name='Mel Brooks')
        os.chdir('..')
        shutil.rmtree('exa_updated')
        os.chdir('..')
        self.assertTrue(os.path.exists('project'))
        example = self.manager.examples.get_example_by_example_id(1)
        self.assertEqual(example.title, 'Updated title')
        self.assertEqual(example.author_name, 'Mel Brooks')
        self.assertEqual(example.shortname, 'exA')
        ex_rst = os.path.join(self.manager.sphinx_example_manager.examples_folder_path, '01', 'exA.rst')
        self.assertTrue(os.path.isfile(ex_rst))
        with open(ex_rst, 'r') as f:
            rst_current = f.read()
        self.assertIn('PDB ID 9ggs', rst_current)

    def test_example_manager_update_example_rename(self):
        cwd = os.getcwd()
        self._build_example_set()
        os.chdir('userspace')
        os.mkdir('exa_updated')
        os.chdir('exa_updated')
        with open('../exA.yaml', 'r') as f:
            lines= f.readlines()
        with open('exAA.yaml', 'w') as f:
            for line in lines:
                if 'title:' in line:
                    f.write('title: exAA_updated\n')
                else:
                    f.write(line)
        self.manager.update_example(1,'exAA')
        os.chdir('..')
        shutil.rmtree('exa_updated')
        os.chdir(cwd)
        example = self.manager.examples.get_example_by_example_id(1)
        self.assertEqual(example.title, 'exAA_updated')
        self.assertEqual(example.shortname, 'exAA')
        ex_rst = os.path.join(self.manager.sphinx_example_manager.examples_folder_path, '01', 'exAA.rst')
        self.assertTrue(os.path.isfile(ex_rst))
        with open(ex_rst, 'r') as f:
            rst_current = f.read()
        self.assertIn('examples/01/inputs/exAA.yaml', rst_current)
        