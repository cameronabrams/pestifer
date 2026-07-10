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
