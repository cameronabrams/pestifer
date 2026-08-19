# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
pytest fixtures for the integration tests.

The parsing and invariant code lives in :mod:`tests.integration.helpers` so tests can import it
directly; this module only exposes it as a fixture for tests that prefer injection.
"""

import pytest

from .helpers import assert_psf_sane, parse_psf, parse_psf_pdb, parse_pdb_coords  # noqa: F401


@pytest.fixture
def psf_sanity():
    """The structural invariant battery, injected rather than imported."""
    return assert_psf_sane
