# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Custom exception hierarchy for pestifer.  Catching PestiferError at the CLI
level allows clean user-facing messages without tracebacks for expected build
failures.
"""

class PestiferError(Exception):
    """Base class for all expected pestifer failures."""

class PestiferBuildError(PestiferError):
    """Raised when a build task reaches an unresolvable condition."""
