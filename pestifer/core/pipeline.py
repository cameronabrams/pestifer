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
    file: object  # e.g., PDBFile, CIFFile, PSFFile, COORFile, XSCFile, VELFile, NAMDConfigFile, CSVFile
    type: str     # e.g., 'intermediate', 'final', 'log', etc.
    produced_by: str
    propagate: bool  # flag for whether to pass to next task

class PipelineContext:
    def __init__(self):
        self.artifacts: dict[str, Artifact] = {}
        self.history: list[dict] = []

    def register(self, key, file_obj, produced_by, *, type='intermediate', propagate=True):
        self.artifacts[key] = Artifact(
            key=key,
            file=file_obj,
            type=type,
            produced_by=produced_by,
            propagate=propagate
        )

    def get(self, key):
        return self.artifacts[key].file

    def get_propagated_inputs(self):
        return {k: a.file for k, a in self.artifacts.items() if a.propagate}
