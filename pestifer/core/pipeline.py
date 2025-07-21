# Author: ChatGPT with modifications by Cameron F. Abrams, <cfa22@drexel.edu>
"""
Pipeline context for managing passing of information from one task to another.
"""

from dataclasses import dataclass
from pathlib import Path
import logging
logger=logging.getLogger(__name__)

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
class DCDFile:
    path: Path

@dataclass
class XSTFile:
    path: Path

@dataclass
class XSTFile:
    path: Path

@dataclass
class NAMDConfigFile:
    path: Path

@dataclass
class CVFile:
    path: Path

@dataclass
class CSVFile:
    path: Path

@dataclass
class LogFile:
    path: Path

@dataclass
class TclFile:
    path: Path

@dataclass
class DATAFile:
    path: Path

@dataclass
class PackMolInputFile:
    path: Path

@dataclass
class TopoFileName:
    name: str

@dataclass
class ParamFileName:
    name: str

@dataclass
class PQRFile:
    path: Path

@dataclass
class ImageFile:
    path: Path

@dataclass
class Artifact:
    key: str
    value: object
    value_type: type
    produced_by: object

    def append(self,new_value,unique=False):
        if isinstance(self.value, list):
            if unique and new_value in self.value:
                return
            self.value.append(new_value)
        else:
            if unique and self.value == new_value:
                return
            self.value = [self.value, new_value]

class PipelineContext:
    def __init__(self):
        self.artifacts: dict[str, Artifact] = {}
        self.history: list[dict] = []

    def register(self, key, value, value_type, produced_by):
        """ 
        Create and register an artifact in the pipeline context.

        Parameters
        ----------
        key : str
            Unique identifier for the artifact.
        value : object
            The value of the artifact, can be any object.
        value_type : type
            The expected type of the value.
        produced_by : object
            The task or object that produced this artifact.
        """
        if key in self.artifacts:
            """ move to history if already exists """
            self.history.append(self.artifacts[key])
        logger.debug(f'Registering artifact {key} of type {value_type} produced by {produced_by.index}')
        self.artifacts[key] = Artifact(
            key=key,
            value=value,
            value_type=value_type,
            produced_by=produced_by
        )

    def get(self, key: str, expected_type=None):
        artifact = self.artifacts.get(key,None)
        if artifact is None:
            return None
        if expected_type and not isinstance(artifact.value, expected_type):
            raise TypeError(f"{key} is not a {expected_type}")
        return artifact.value

    def get_artifact_series(self, key: str = ''):
        """
        Retrieve a series of artifacts with the same key.
        
        Parameters
        ----------
        key : str
            The key for the artifacts to retrieve.
        produced_by : object, optional
            The task or object that produced the artifacts to retrieve.

        Returns
        -------
        list
            A list of artifact values associated with the given key, in reverse registration order.
        """
        history = [h.value for h in self.history if h.key == key]
        current = self.artifacts.get(key, None)
        if current:
            return history + [current.value]
        return history
    
    def get_artifact_values_collection_as_dict(self, produced_by=None):
        """
        Retrieve a collection of artifacts produced by a specific task or object.
        
        Parameters
        ----------
        produced_by : object, optional
            The task or object that produced the artifacts to retrieve.

        Returns
        -------
        dict
            A dictionary of artifact values associated with the given produced_by, keyed by artifact key.
        """
        history = {h.key: h.value for h in self.history if (not produced_by or h.produced_by == produced_by)}
        current = {k: a.value for k, a in self.artifacts.items() if (not produced_by or a.produced_by == produced_by)}
        return {**history, **current}

    def context_to_string(self):
        """ 
        Report the current artifacts in the pipeline context.
        
        Returns
        -------
        str
            A string representation of the current artifacts in the pipeline context.
        """
        return "\n".join([f"{k}: {a.value}" for k, a in self.artifacts.items()])

    def history_to_string(self):
        """ 
        Report the history of artifacts in the pipeline context.
        
        Returns
        -------
        str
            A string representation of the history of artifacts in the pipeline context.
        """
        return "\n".join([f"{h.key}: {h.value} (produced by {h.produced_by.index})" for h in self.history])
