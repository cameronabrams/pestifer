# Author: ChatGPT with modifications by Cameron F. Abrams, <cfa22@drexel.edu>
"""
Pipeline context for managing passing of information from one task to another.
"""

from dataclasses import dataclass
from pathlib import Path

@dataclass
class PDBFile:
    path: Path

@dataclass
class CIFFile:
    path: Path

@dataclass
class PSFFile:
    path: Path

@dataclass
class COORFile:
    path: Path

@dataclass
class XSCFile:
    path: Path

@dataclass
class VELFile:
    path: Path

@dataclass
class NAMDConfigFile:
    path: Path

@dataclass
class CSVFile:
    path: Path

@dataclass
class Artifact:
    key: str
    value: object
    value_type: type
    type: str     # e.g., 'intermediate', 'final', 'log', etc.
    produced_by: object
    propagate: bool  # flag for whether to pass to next task

class PipelineContext:
    def __init__(self):
        self.artifacts: dict[str, Artifact] = {}
        self.history: list[dict] = []

    def register(self, key, value, value_type, produced_by, *, type='intermediate', propagate=True):
        self.artifacts[key] = Artifact(
            key=key,
            value=value,
            value_type=value_type,
            type=type,
            produced_by=produced_by,
            propagate=propagate
        )

    def get(self, key: str, expected_type=None):
        artifact = self.artifacts.get(key,None)
        if artifact is None:
            return None
        if expected_type and not isinstance(artifact.value, expected_type):
            raise TypeError(f"{key} is not a {expected_type}")
        return artifact.value

    def get_propagated_inputs(self):
        return {k: a.value for k, a in self.artifacts.items() if a.propagate}
