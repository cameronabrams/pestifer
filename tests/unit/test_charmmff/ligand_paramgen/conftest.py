# These tests import pestifer.charmmff.ligand_paramgen, which requires `rdkit` (from the optional
# `ligand-paramgen` extra). Skip the directory when rdkit is not installed rather than erroring
# on import; it collects and runs where the extra is present.
import importlib.util

if importlib.util.find_spec("rdkit") is None:
    collect_ignore_glob = ["test_*.py"]
