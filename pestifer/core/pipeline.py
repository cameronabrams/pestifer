# Author: ChatGPT with modifications by Cameron F. Abrams, <cfa22@drexel.edu>
"""
Pipeline context for managing passing of information from one task to another via artifacts.
"""
from __future__ import annotations
import logging

from .artifacts import Artifact, FileArtifact, ArtifactDict, ArtifactList, FileArtifactList, StateArtifacts

logger = logging.getLogger(__name__)

class PipelineContext:
    def __init__(self, controller_index: int = 0):
        self.head: ArtifactDict = ArtifactDict()
        self.history: ArtifactList = ArtifactList()
        self.controller_index = controller_index

    def __repr__(self):
        return f"PipelineContext(controller_index={self.controller_index})"

    def register(self, artifact: Artifact | ArtifactList | ArtifactDict, key: str = None, requesting_object: object | None = None):
        """ 
        Create and register an artifact in the pipeline context.

        Parameters
        ----------
        artifact : Artifact
            The artifact to register.
        key : str
            Unique identifier for the artifact; if not provided, it will be derived from the artifact's key attribute.
        requesting_object : object, optional
            The object that is requesting the registration of the artifact. This can be used for logging or tracking purposes.
        """
        if key: # If a key is provide, override the artifact's key
            artifact.key = key
        else:
            if hasattr(artifact, 'key'):
                key = artifact.key
        if not key:
            raise ValueError("Key must be provided for the artifact.")
        if key in self.head:
            """ move to history if already exists """
            self.history.append(self.head[key])
        logger.debug(f'Registering artifact {type(artifact)}')
        # stamp the arifact if it isn't already stamped
        if not artifact.has_stamp():
            # logger.debug(f'Stamping artifact {key} with requesting object {requesting_object}')
            if requesting_object:
                artifact.stamp(requesting_object)
            else:
                artifact.stamp(self)
        self.head[key] = artifact
    
    def bury(self, artifact: Artifact):
        """
        Bury an artifact in the history without registering it as a current artifact
        
        Parameters
        ----------
        artifact : Artifact
            The artifact to bury.
        """
        self.history.append(artifact)

    def get_current_artifact(self, key: str):
        return self.head.get(key, None)

    def get_current_artifact_data(self, key: str):
        """
        Retrieve the data of an artifact by its key.

        Parameters
        ----------
        key : str
            The key for the artifact to retrieve.

        Returns
        -------
        Any
            The value of the artifact, or None if not found.
        """
        artifact = self.get_current_artifact(key)
        if artifact:
            return artifact.data
        return None

    def get_artifact_series_by_key(self, key: str = '') -> ArtifactList | FileArtifactList:
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
        ArtifactList
            A list of artifacts associated with the given key, in reverse registration order.
        """
        series = ArtifactList([h for h in self.history if h.key == key])
        current = self.get_current_artifact(key)
        if current:
            series.append(current)  # most recent at end
        if all([isinstance(x, FileArtifact) for x in series]):
            return FileArtifactList(series)
        return series

    def get_artifact_collection_as_lists(self, produced_by: object | None = None) -> dict[str, ArtifactList | ArtifactList | FileArtifactList]:
        """
        Retrieve a collection of artifacts produced by a specific task or object.
        
        Parameters
        ----------
        produced_by : object | None
            The task or object that produced the artifacts to retrieve.

        Returns
        -------
        dict[str, ArtifactList | ArtifactList | FileArtifactList]
            A dict containing:
            - A list of all data artifacts.
            - A list of all ArtifactFileList artifacts.
            - A list of file artifacts not in ArtifactFileList artifacts.
        """
        data_artifacts = ArtifactList()
        filelist_artifacts = ArtifactList()
        file_artifacts = FileArtifactList()
        # logger.debug(f'Getting artifact collection for produced_by={repr(produced_by)} from history of {len(self.history)} artifacts and current head of {len(self.head)} artifacts.')
        # for h in self.history:
        #     logger.debug(f'History artifact: {repr(h)}')
        # for a in self.head.values():
        #     logger.debug(f'Current artifact: {repr(a)}')
        if produced_by is None:
            # If no produced_by is specified, use all artifacts in history and head
            my_history = self.history
            my_current = self.head
        else:
            my_history = self.history.filter_by_produced_by(produced_by=produced_by)
            my_current = self.head.filter_by_produced_by(produced_by=produced_by)
        # logger.debug(f'Filtered history artifacts: {len(my_history)}')
        # logger.debug(f'Filtered current artifacts: {len(my_current)}')
        # history = [h for h in self.history if (not produced_by or h.produced_by == produced_by)]
        # current = [a for a in self.head.values() if (not produced_by or a.produced_by == produced_by)]
        all_artifacts = my_history + my_current.to_list()
        # logger.debug(f'Found {len(all_artifacts)} artifacts in the pipeline context for produced_by={repr(produced_by)}.')
        for artifact in all_artifacts:
            if isinstance(artifact, FileArtifactList):
                filelist_artifacts.append(artifact)
            elif isinstance(artifact, FileArtifact):
                file_artifacts.append(artifact)
            elif isinstance(artifact, StateArtifacts):
                file_artifacts.extend(artifact.to_list())
            else:
                data_artifacts.append(artifact)

        # logger.debug(f'Found {len(data_artifacts)} data artifacts, {len(filelist_artifacts)} file list artifacts, and {len(file_artifacts)} file artifacts.')

        return {'data': data_artifacts, 'filelists': filelist_artifacts, 'files': file_artifacts}

    def context_to_string(self) -> str:
        """ 
        Report the current artifacts in the pipeline context.
        
        Returns
        -------
        str
            A string representation of the current artifacts in the pipeline context.
        """
        return "\n".join([f"{k}: {a.value}" for k, a in self.head.items()])

    def history_to_string(self) -> str:
        """ 
        Report the history of artifacts in the pipeline context.
        
        Returns
        -------
        str
            A string representation of the history of artifacts in the pipeline context.
        """
        return "\n".join([f"{h.key}: {h.value} (produced by {h.produced_by})" for h in self.history])

    def stash(self, key: str) -> str:
        """
        Stash the current artifact under a new key.
        
        Parameters
        ----------
        key : str
            The key under which to stash the current artifact.
        """
        artifact = self.head.pop(key, None)
        if artifact:
            self.history.append(artifact)

    def import_artifacts(self, other: PipelineContext):
        """
        Import artifacts from another pipeline context.
        """
        logger.debug(f'Importing artifacts from {repr(other)} into {repr(self)}')
        for other_artifact_key, other_artifact in other.head.items():
            if not other_artifact_key in self.head:
                self.head[other_artifact_key] = other_artifact
            else:
                self.history.append(other_artifact)
        self.history.extend(other.history)