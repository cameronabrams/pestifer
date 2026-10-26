from argparse import Namespace
from pestifer.subcommands.build import RunSubcommand


class TestBuildDefaultLogFile:
    """The build/run subcommand derives its diagnostics-log name from the input config
    so successive builds in one directory get unique logs instead of clobbering a shared
    ``pestifer_diagnostics.log``."""

    def test_derives_log_from_config_basename(self):
        sc = RunSubcommand()
        assert sc.default_log_file(Namespace(config='foo.yaml')) == 'foo-diagnostics.log'

    def test_strips_directory_component(self):
        sc = RunSubcommand()
        assert sc.default_log_file(Namespace(config='some/dir/helper-01-base.yaml')) == 'helper-01-base-diagnostics.log'

    def test_extensionless_config(self):
        sc = RunSubcommand()
        assert sc.default_log_file(Namespace(config='bar')) == 'bar-diagnostics.log'

    def test_falls_back_to_static_default_without_config(self):
        sc = RunSubcommand()
        assert sc.default_log_file(Namespace(config=None)) == sc.log_file
        assert sc.default_log_file(Namespace()) == sc.log_file


class TestAFailedBuildKeepsItsSeeds:
    """A build that stops writes ``run-record-failed.json``.  The 3.29.8 sweep lost example 14
    replica 1's seeds to a terminate refusal, because the record was written only on success.

    Driven through ``RunSubcommand.func`` itself, with the Controller replaced, because the
    helper that builds the record is not what failed: the call site never asked for it.
    """

    def _run(self, tmp_path, monkeypatch, *, raise_in_do_tasks=False, exit_code=0):
        import json
        from types import SimpleNamespace
        from unittest import mock
        from pestifer.core.errors import PestiferBuildError
        import pestifer.subcommands.build as B

        monkeypatch.chdir(tmp_path)
        (tmp_path / 'x.yaml').write_text('title: x\n')
        config = mock.MagicMock()
        config.__getitem__ = lambda _s, k: {'namd': {'seed': 27021972}} if k == 'user' else {}
        config.userfile = 'x.yaml'
        config.namd_type = 'cpu'
        config.configure_new.return_value = config
        ligate = SimpleNamespace(index=3, taskname='ligate', outcome={'ccd_seed': 1234},
                                 duration=1.0)
        terminate = SimpleNamespace(index=12, taskname='terminate', outcome={}, duration=1.0)

        class FakeController:
            def configure(self, *_a, **_k):
                self.tasks = [ligate, terminate]
                self.exit_code = 0
                self.current_task = None
                return self

            def do_tasks(self):
                self.current_task = terminate
                if raise_in_do_tasks:
                    raise PestiferBuildError('terminate: RING_PIERCING_FOUND refusing ...')
                self.exit_code = exit_code
                return {}

        monkeypatch.setattr(B, 'Config', lambda **_k: config)
        monkeypatch.setattr(B, 'Controller', FakeController)
        monkeypatch.setattr(B, 'log_environment', lambda _c: {})
        monkeypatch.setattr(B, 'log_citations', lambda _c: [])
        args = Namespace(config='x.yaml', output_dir='./', ncpus=1, gpu=False, seed=None,
                         check=False, restart=False, fresh=False, from_task=None,
                         complete_config=False, controller_index=0)
        if raise_in_do_tasks:
            import pytest
            with pytest.raises(PestiferBuildError):
                RunSubcommand.func(args)
        else:
            RunSubcommand.func(args)
        failed = tmp_path / 'run-record-failed.json'
        done = tmp_path / 'run-record.json'
        return (json.loads(failed.read_text()) if failed.exists() else None,
                json.loads(done.read_text()) if done.exists() else None)

    def test_a_raised_refusal_writes_the_failed_record(self, tmp_path, monkeypatch):
        failed, done = self._run(tmp_path, monkeypatch, raise_in_do_tasks=True)
        assert done is None
        assert failed['status'] == 'failed'
        assert failed['exit_code'] != 0
        assert failed['failed_task'] == {'index': 12, 'task': 'terminate'}
        assert 'RING_PIERCING_FOUND' in failed['error']
        assert failed['seed'] == 27021972
        assert {'index': 3, 'task': 'ligate', 'ccd_seed': 1234} in failed['protocol']
        assert failed['system'] == {}

    def test_a_nonzero_result_writes_the_failed_record(self, tmp_path, monkeypatch):
        failed, done = self._run(tmp_path, monkeypatch, exit_code=1)
        assert done is None
        assert failed['status'] == 'failed'
        assert failed['failed_task']['task'] == 'terminate'

    def test_a_finished_build_writes_only_the_ordinary_record(self, tmp_path, monkeypatch):
        (tmp_path / 'run-record-failed.json').write_text('{}')   # from an earlier attempt
        failed, done = self._run(tmp_path, monkeypatch)
        assert failed is None
        assert done['status'] == 'completed'
