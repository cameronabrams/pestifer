# These tests exercise pestifer's Sphinx extension (pestifer.sphinxext), which imports
# `sphinx`/`docutils`. Sphinx is a docs-build dependency, not a runtime or test one, so on an
# environment without it (e.g. a tool-free CI runner) skip this directory instead of failing
# to import. Where sphinx is installed (a dev checkout), the tests collect and run normally.
import importlib.util

if importlib.util.find_spec("sphinx") is None:
    collect_ignore_glob = ["test_*.py"]
