# Author: ChatGPT with modifications by Cameron F. Abrams, <cfa22@drexel.edu>
"""
Pipeline context for managing passing of information from one task to another via artifacts.
"""
from __future__ import annotations
import logging

from .artifacts import Artifact, FileArtifact, ArtifactDict, ArtifactList, FileArtifactDict, FileArtifactList, StateArtifacts

logger = logging.getLogger(__name__)

class PipelineContext:
    def __init__(self, controller_index: int = 0):
        self.head: ArtifactDict = ArtifactDict(key='Head')
        self.history: ArtifactList = ArtifactList(key='History')
        self.controller_index = controller_index

    def __repr__(self):
        return f"PipelineContext(controller_index={self.controller_index})"

    def register(self, data: object, key: str, requestor: object, artifact_type: type = Artifact | ArtifactList | ArtifactDict, **kwargs) -> Artifact:
        """
        Artifact registrar.  If an artifact with the requested key already exists, and the data is the same, the existing artifact
        is stamped with the requestor and returned.  If the data is different, the existing artifact is moved to history and a new artifact
        is created, stamped with the requestor, registered as the current artifact, and returned.
        If no artifact with the requested key exists, a new artifact is created, stamped with the requestor, registered as the current artifact, and returned."""

        if key in self.head:
            existing_artifact = self.head[key]
            if existing_artifact.data == data:
                # same data, just stamp and return
                existing_artifact.stamp(requestor)
                logger.debug(f'Artifact with key {key} already exists with same data; stamped with {requestor}.')
                return existing_artifact
            else:
                # different data, move existing to history and create new artifact
                self.history.append(existing_artifact)
                logger.debug(f'Artifact with key {key} already exists with different data; moved to history.')
        # create new artifact
        new_aritifact = artifact_type(data=data, key=key, **kwargs).stamp(requestor)
        self.head[key] = new_aritifact
        logger.debug(f'Registered new artifact with key {key} produced by {requestor}; optional parameters: {kwargs}.')
        return new_aritifact
    
    def register_if_exists(self, data: object, requestor: object, artifact_type: type = FileArtifact, **kwargs) -> FileArtifact | None:
        """
        Artifact registrar that only registers the artifact if it exists (for file artifacts).
        If an artifact with the requested key already exists, and the data is the same, the existing artifact
        is stamped with the requestor and returned.  If the data is different, the existing artifact is moved to history and a new artifact
        is created, stamped with the requestor, registered as the current artifact, and returned.
        If no artifact with the requested key exists, a new artifact is created, stamped with the requestor, registered as the current artifact, and returned.

        Parameters
        ----------
        data : object
            The data for the artifact to register.
        key : str
            Unique identifier for the artifact.
        requestor : object
            The object that is requesting the registration of the artifact. This can be used for logging or tracking purposes.
        artifact_type : type
            The type of artifact to create (default is FileArtifact).
        **kwargs : dict
            Additional keyword arguments to pass to the artifact constructor.

        Returns
        -------
        FileArtifact | None
            The registered artifact if it exists, otherwise None.
        """
        fa = artifact_type(data, **kwargs)
        if not isinstance(fa, FileArtifact):
            raise ValueError(f"Artifact type {artifact_type} is not a FileArtifact; cannot use register_if_exists.")
        if fa.exists():
            return self.register(data=data, key=fa.key, requestor=requestor, artifact_type=artifact_type, **kwargs)
        else:
            logger.debug(f'File artifact {fa.name} does not exist; not registering.')
            return None

    def rekey(self, old_key: str, new_key: str):
        """
        Change the key of an existing artifact in the head.

        Parameters
        ----------
        old_key : str
            The current key of the artifact to rekey.
        new_key : str
            The new key to assign to the artifact.

        Raises
        ------
        ValueError
            If the old_key does not exist in the head or if the new_key already exists in the head.
        """
        if old_key not in self.head:
            raise ValueError(f"Cannot rekey artifact: old key '{old_key}' does not exist.")
        if new_key in self.head:
            raise ValueError(f"Cannot rekey artifact: new key '{new_key}' already exists.")
        artifact = self.head.pop(old_key)
        artifact.key = new_key
        self.head[new_key] = artifact
        logger.debug(f'Rekeyed artifact from {old_key} to {new_key}.')

    # def register(self, artifact: Artifact | ArtifactList | ArtifactDict, key: str = None, requesting_object: object = None):
    #     """ 
    #     Register an artifact in the pipeline context.  If a key is provided, overwrite the artifact's key.
    #     If 'living' is True, and an artifact with the same key already exists, the exising artifact is
    #     updated with the contents of the new artifact.  Otherise, the existing artifact is moved to history
    #     and the new artifact is registered as the current artifact.

    #     Parameters
    #     ----------
    #     artifact : Artifact
    #         The artifact to register.
    #     key : str
    #         Unique identifier for the artifact; if not provided, it will be derived from the artifact's key attribute.
    #     requesting_object : object, optional
    #         The object that is requesting the registration of the artifact. This can be used for logging or tracking purposes.
    #     """
    #     if key and hasattr(artifact, 'key'): # If a key is provide, override the artifact's key
    #         artifact.key = key
    #     else:
    #         if hasattr(artifact, 'key'):
    #             key = artifact.key
    #     if not key:
    #         raise ValueError(f"{type(artifact).__name__} has no key attribute; key must be provided for the artifact as a keyword argument to register().")
        
    #     if key in self.head:
    #         """ move to history if already exists """
    #         logger.debug(f'Stashed artifact {key} in history')
    #         self.history.append(self.head[key])
    #     # stamp the artifact if it isn't already stamped
    #     # logger.debug(f'Stamping artifact {key} with requesting object {requesting_object}')
    #     self.head[key] = artifact
    #     logger.debug(f'Registered artifact {repr(artifact)} at key {key}')

    def show_artifacts(self, header: str = 'Current Artifacts'):
        """
        Debugging utility to show current and historical artifacts.
        """
        logger.debug('*'*72)
        logger.debug(header)
        logger.debug('*'*72)
        logger.debug(f'Head:')
        self.show_artifact(self.head)
        logger.debug(f'History:')
        self.show_artifact(self.history)
        logger.debug('*'*72)

    def show_artifact(self, artifact: Artifact, depth = 1):
        if isinstance(artifact, ArtifactList):
            logger.debug(f'{"    "*depth}- "{artifact.key}" ({id(artifact)}) List with {len(artifact)} items:')
            for item in artifact:
                self.show_artifact(item, depth + 1)
        elif isinstance(artifact, ArtifactDict):
            logger.debug(f'{"    "*depth}- "{artifact.key}" ({id(artifact)}) Dict with {len(artifact)} items:')
            for key, item in artifact.items():
                self.show_artifact(item, depth + 1)
        elif isinstance(artifact, Artifact):
            pytestable_str = ''
            if hasattr(artifact, 'pytestable') and artifact.pytestable:
                pytestable_str = ' (pytestable)'
            logger.debug(f'{"    "*depth}- {artifact.key} ({id(artifact)}): data={artifact.data}{pytestable_str}')
        else:
            logger.debug(f'{"    "*depth}- ***Unknown artifact type: {type(artifact)}***')

    def bury(self, artifact: Artifact):
        """
        Bury an artifact in the history without registering it as a current artifact
        
        Parameters
        ----------
        artifact : Artifact
            The artifact to bury.
        """
        self.history.append(artifact)

    def get_current_artifact(self, key: str, **kwargs) -> Artifact | None:
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

    def get_all_file_artifacts(self, produced_by: object | None = None) -> FileArtifactList:
        """
        Retrieve a collection of artifacts produced by a specific task or object.
        
        Parameters
        ----------
        produced_by : object | None
            The task or object that produced the artifacts to retrieve.

        Returns
        -------
        FileArtifactList
            A list of all file artifacts produced by the specified task.
        """
        file_artifacts = FileArtifactList()
        if produced_by is None:
            # If no produced_by is specified, use all artifacts in history and head
            my_history = self.history
            my_current = self.head
        else:
            my_history = self.history.filter_by_produced_by(produced_by=produced_by)
            my_current = self.head.filter_by_produced_by(produced_by=produced_by)
        all_artifacts = my_history + my_current.to_list()
        for artifact in all_artifacts:
            if isinstance(artifact, FileArtifact):
                file_artifacts.append(artifact)
            elif isinstance(artifact, FileArtifactList):
                file_artifacts.extend(artifact.data)
            elif isinstance(artifact, FileArtifactDict):
                file_artifacts.extend(artifact.data.values())
            else:
                logger.debug(f'Ignoring non-file artifact "{artifact.key}" of type {type(artifact)}')
        file_artifacts = file_artifacts.unique_paths()

        return FileArtifactList([x for x in file_artifacts if x.exists()])

    def context_to_string(self) -> str:
        """ 
        Report the current artifacts in the pipeline context.
        
        Returns
        -------
        str
            A string representation of the current artifacts in the pipeline context.
        """
        return "\n".join([f"{k}: {a.data}" for k, a in self.head.items()])

    def history_to_string(self) -> str:
        """ 
        Report the history of artifacts in the pipeline context.
        
        Returns
        -------
        str
            A string representation of the history of artifacts in the pipeline context.
        """
        return "\n".join([f"{h.key}: {h.data} (produced by {h.produced_by})" for h in self.history])

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