# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
A class for handling artifacts in the Pestifer core. Artifacts are files or 
data generated during the execution of tasks and are managed by the 
:class:`~pestifer.core.pipeline.PipelineContext`.

Tasks are the primary creators of Artifacts, and the pipeline context is responsible 
for managing their lifecycles. Any task may create and register an Artifact.  All 
tasks have a ``register`` method that interfaces with the pipeline to register an 
Artifact, and a ``get_current_artifact`` method to retrieve an Artifact by its key. 

All Artifacts must have a "key" that the pipeline uses to track them.  A key value is normally created when an Artifact is created:

.. code-block::python

  artifact = Artifact(data=my_data, key="my_artifact_key")
  self.register(artifact)

Any artifact can be retrieved from the pipeline context, in a later task, say, using its key:

.. code-block::python

  retrieved_artifact = self.get_current_artifact("my_artifact_key")

Containers of Artifacts can also be registered and retrieved.  However, importantly, Artifact containers do not have a single ``key`` attribute.  A key is associated with an Artifact container when a task calls its register method:

.. code-block::python

  artifact_container = ArtifactDict(a1=Artifact(data=my_data, key="my_artifact_key"), a2=Artifact(data=my_data2, key="my_artifact_key2"))
  self.register(artifact_container, key='my_artifact_container_key')

Artifact containers allow for groups of artifacts to be associated with a single key, making it easier to manage related artifacts, and to separate them from other artifacts in the pipeline.

The data in an artifact can be any object, but most useful are Paths.  Artifacts with data that are Paths are called FileArtifacts.  When a task creates a file, it is a good idea to create a FileArtifact and register it.

.. code-block::python

  file_artifact = FileArtifact(data=file_path, key="my_file_artifact_key")
  self.register(file_artifact)

This is how the pipeline is used to "pass" files from one task to any other.  Because each Artifact's registration associates it with a specific task instance, later tasks can query the pipeline for artifacts created by previous tasks using their keys and task ids.

"""
from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
import difflib
import logging
import os
import tarfile

from unidiff import PatchSet

from pestifer.core.labels import Labels

from ..molecule.atom import AtomList
from ..psfutil.psfcontents import PSFContents
from ..util.stringthings import my_logger

from pidibble.pdbparse import PDBParser

logger = logging.getLogger(__name__)

@dataclass
class Artifact():
    """
    Base Artifact class.
    """
    data: object | None = None
    """ The data contained in the artifact. """
    key: str | None = None
    """ A unique identifier for the artifact. """
    produced_by: object | None = None
    """ The task that generated the artifact. """    
    provenance: list[object] = field(default_factory=list)
    """ A list of tasks that registered the artifact, ordered chronologically. """

    def __repr__(self) -> str:
        retstr = f'{self.__class__.__name__}(key=\'{self.key}\', produced_by=\'{self.produced_by}\', type(data)={type(self.data)})'
        return retstr
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Artifact):
            return False
        if not self.__class__.__name__ == other.__class__.__name__:
            return False
        return self.key == other.key and self.data == other.data

    def copy(self, **kwargs) -> Artifact:
        """
        Create a copy of the artifact, optionally overriding attributes with keyword arguments.

        Parameters
        ----------
        **kwargs : dict
            Attributes to override in the copied artifact.

        Returns
        -------
        Artifact
            A new Artifact instance with the same attributes as the original, except for any overridden by kwargs.
        """
        attrs = {
            'data': self.data,
            'key': self.key,
            'produced_by': self.produced_by,
            'provenance': self.provenance.copy()
        }
        attrs.update(kwargs)
        return self.__class__(**attrs)

    def has_stamp(self) -> bool:        
        """
        Check if the artifact has a stamp (owner information).
        
        Returns
        -------
        bool
            True if the artifact has a stamp, False otherwise.
        """
        return self.produced_by is not None

    def stamp(self, owner: object) -> Artifact:
        """
        Stamp the artifact with the owner information.
        Any attributes that are instances of Artifact are also stamped.

        Parameters
        ----------
        owner : object
            The owner of the artifact, which can be any object that produced this artifact.
        """
        if owner is None:
            raise ValueError(f"Owner of {repr(self)} (type {type(self)}) must not be None.")
        if self.produced_by is owner:
            return self
        if self.produced_by is not None:
            self.provenance.append(self.produced_by)
        self.produced_by = owner
        return self

@dataclass
class DataArtifact(Artifact):
    """
    Represents a data artifact in the Pestifer core.
    """
    description: str | None = None

@dataclass
class ArtifactList(Artifact):
    """
    A list of Artifacts.
    """
    data: list[Artifact] = field(default_factory=list)

    def __getitem__(self, index: int) -> Artifact:
        return self.data[index]
    
    def __setitem__(self, index: int, value: Artifact) -> None:
        self.data[index] = value

    def __delitem__(self, index: int) -> None:
        del self.data[index]
    
    def __iter__(self):
        return iter(self.data)
    
    def __len__(self) -> int:
        return len(self.data)

    def __add__(self, other: ArtifactList) -> ArtifactList:
        if not isinstance(other, ArtifactList):
            raise TypeError(f"Can only add ArtifactList to ArtifactList, not {type(other)}")
        new_list = ArtifactList(self.data + other.data, produced_by=self.produced_by, key=self.key)
        return new_list

    def append(self, item: Artifact) -> None:
        """
        Append an Artifact to the list, if it is not already present.
        """
        if not isinstance(item, Artifact):
            raise TypeError(f"Expected Artifact, got {type(item)}")
        if item not in self.data:
            self.data.append(item)
            self.produced_by = item.produced_by
        else:
            logger.debug(f"ArtifactList: not appending duplicate {item}")

    def extend(self, items: list[Artifact]) -> None:
        """
        Extend the list by appending elements from the iterable.
        """
        for item in items:
            self.append(item)

    def remove(self, item: Artifact) -> None:
        """
        Remove an Artifact from the list.
        """
        if not isinstance(item, Artifact):
            raise TypeError(f"Expected Artifact, got {type(item)}")
        if item in self.data:
            self.data.remove(item)
            self.produced_by = None
        else:
            logger.debug(f"ArtifactList: not removing non-existent {item}")
    
    def sort(self, *, key=None, reverse: bool = False) -> None:
        """
        Sort the list in place.
        """
        self.data.sort(key=key, reverse=reverse)

    @classmethod
    def from_dict(cls, artifact_dict: ArtifactDict) -> ArtifactList:
        """
        Create an ArtifactList from an ArtifactDict.
        """
        artifact_list = cls()
        for artifact in artifact_dict.values():
            artifact_list.append(artifact)
        return artifact_list
    
    def to_dict(self) -> ArtifactDict:
        """
        Create an ArtifactDict from the ArtifactList.
        """
        return ArtifactDict.from_list(self)

    def stamp(self, owner: object) -> ArtifactList:
        """
        Stamp all artifacts in the list with the owner information.
        
        Parameters
        ----------
        owner : object
            The owner of the artifacts, which can be any object that produced these artifacts. If not provided, the artifacts will not be stamped.
        
        Returns
        -------
        ArtifactList
            The artifact list with all artifacts stamped with the owner information.
        """
        self.produced_by = owner
        for artifact in self.data:
           artifact.stamp(owner)
        return self

    def filter_by_produced_by(self, produced_by: object) -> ArtifactList:
        """
        Filter the artifact list by the task that produced them.

        Parameters
        ----------
        produced_by : object
            The task or object that produced the artifacts to filter by.

        Returns
        -------
        ArtifactList
            A new ArtifactList containing only the artifacts produced by the specified task.

        """
        return ArtifactList([artifact for artifact in self.data if artifact.produced_by == produced_by])

    def filter_by_artifact_type(self, artifact_type: type[Artifact]) -> ArtifactList:
        """
        Filter the artifact list by the type of artifacts.

        Parameters
        ----------
        artifact_type : type[Artifact]
            The type of artifacts to filter by.
        
        Returns
        -------
        ArtifactList
            A new ArtifactList containing only the artifacts of the specified type.

        """
        return ArtifactList([artifact for artifact in self.data if isinstance(artifact, artifact_type)])
    
    def filter_by_key(self, key: str) -> ArtifactList:
        """
        Filter the artifact list by the key of artifacts.

        Parameters
        ----------
        key : str
            The key of artifacts to filter by.
        
        Returns
        -------
        ArtifactList
            A new ArtifactList containing only the artifacts with the specified key.

        """
        return ArtifactList([artifact for artifact in self.data if artifact.key == key])

@dataclass
class ArtifactDict(Artifact):
    """
    Dictionary of Artifacts.
    """
    data: dict[str, Artifact] = field(default_factory=dict)

    def __getitem__(self, key: str) -> Artifact:
        return self.data[key]
    
    def __setitem__(self, key: str, value: Artifact) -> None:
        self.data[key] = value

    def __delitem__(self, key: str) -> None:
        del self.data[key]
    
    def __iter__(self):
        return iter(self.data)
    
    def __len__(self) -> int:
        return len(self.data)

    def clear(self) -> None:
        """
        Clear all artifacts from the dictionary.
        """
        self.data.clear()
    
    def get(self, key: str, default: Artifact | None = None) -> Artifact | None:
        """
        Get an artifact by key, returning a default value if the key is not found.
        """
        return self.data.get(key, default)

    def values(self) -> list[Artifact]:
        """
        Get a list of all artifacts in the dictionary.
        """
        return list(self.data.values())
    
    def keys(self) -> list[str]:
        """
        Get a list of all keys in the dictionary.
        """
        return list(self.data.keys())
    
    def items(self) -> list[tuple[str, Artifact]]:
        """
        Get a list of all key-artifact pairs in the dictionary.
        """
        return list(self.data.items())
    
    def pop(self, key, default = None) -> Artifact:
        """
        Remove and return an arbitrary artifact from the dictionary.
        """
        return self.data.pop(key, default)

    @classmethod
    def from_list(cls, artifact_list: ArtifactList) -> ArtifactDict:
        """
        Create an ArtifactDict from an ArtifactList.
        """
        artifact_dict = cls()
        for artifact in artifact_list:
            artifact_dict[artifact.key] = artifact
        return artifact_dict
    
    def to_list(self) -> ArtifactList:
        """
        Create an ArtifactList from the ArtifactDict.
        """
        return ArtifactList.from_dict(self)

    def stamp(self, owner: object) -> ArtifactDict:
        """
        Stamp all artifacts in the dictionary with the owner information.
        
        Parameters
        ----------
        owner : object
            The owner of the artifacts, which can be any object that produced these artifacts.

        Returns
        -------
        ArtifactDict
            The artifact dictionary with all artifacts stamped with the owner information.
        """
        # if owner is None:
        #     raise ValueError("Owner must not be None.")
        self.produced_by = owner
        for artifact in self.data.values():
            artifact.stamp(owner)
        return self
    
    def update_item(self, artifact: Artifact, key: str | None = None) -> None:
        """
        Updates a new or existing key in the ArtifactDict with an Artifact.
        """
        if key:
            self.data[key] = artifact
        else:
            self.data[artifact.key] = artifact

    def filter_by_produced_by(self, produced_by: object) -> ArtifactDict:
        """
        Filter the artifact dictionary by the producer of artifacts.

        Parameters
        ----------
        produced_by : object
            The producer of artifacts to filter by.

        Returns
        -------
        ArtifactDict
            A new ArtifactDict containing only the artifacts produced by the specified producer.
        """
        return ArtifactDict({key: artifact for key, artifact in self.data.items() if artifact.produced_by == produced_by})

@dataclass
class FileArtifact(Artifact, ABC):
    """
    A class for artifacts that are files.
    """
    data: str | None = None  # name with or without extension
    """ The name of the file artifact, with or without extension. """
    description: str | None = None
    """ A brief description of the file artifact; optional. """
    mime_type: str | None = 'application/octet-stream'
    """ The MIME type of the file artifact; 'application/octet-stream' by default. """
    pytestable: bool = False
    """ Flag indicating whether or not this file artifact is pytestable; i.e., if generated by pytest, there is a gold-standard version available for comparison. """
    nonstate_results: bool = False
    """ Flag indicating whether or not this file artifact is a non-state result file, such as a PNG image of a plot. """

    @property
    @abstractmethod
    def ext(self) -> str:
        """ The file extension of the file artifact, which is used as the artifact key by default. """
        pass

    def __post_init__(self):
        if self.key is None:
            self.key = self.ext
        if self.data and '.' in self.data:
            # correctly assigns data in case this is initialized
            # with a filename that has an extension
            self.data, dot_apparent_ext = os.path.splitext(self.data)
            apparent_ext = dot_apparent_ext.replace('.', '')
            if self.ext != apparent_ext:
                logger.warning(f"Extension mismatch: {self.ext} != {apparent_ext}")

    @property
    def name(self) -> str:
        """ The name of the artifact's file """
        return self.data + '.' + self.ext

    @property
    def path(self) -> Path:
        """ The file path of the file artifact. """
        return Path(self.name)

    def __repr__(self) -> str:
        restr = f"FileArtifact(name={self.name}, produced_by={self.produced_by}"
        if self.pytestable:
            restr += f" *pytestable*"
        restr += ")"
        return restr

    def exists(self) -> bool:
        """
        Check if the file artifact exists.
        """
        return self.path.exists()

    def remove(self) -> None:
        """
        Remove the file artifact.
        """
        if self.exists():
            self.path.unlink()

    def validate(self):
        """
        Validate the file artifact.

        Raises
        ------
        FileNotFoundError
            If the file artifact does not exist.
        """
        if not self.exists():
            raise FileNotFoundError(f"{self.path} not found")

    def diff(self, other: FileArtifact | Path) -> str:
        """
        Compare this file artifact with another for equality.

        Parameters
        ----------
        other : FileArtifact
            The other file artifact to compare against.

        Returns
        -------
        bool
            True if the artifacts are equal, False otherwise.
        """
        if self.mime_type == 'application/octet-stream':
            return f'Cannot compare binary files' # only compare text files
        if not self.pytestable:
            return f'self (type {type(self)}) is not pytestable'
        if not isinstance(other, FileArtifact | Path):
            return f'bad type for other {type(other)}'
        if isinstance(other, FileArtifact):
            if not other.pytestable:
                return f'other (type {type(other)}) is not pytestable'
            if other.mime_type == 'application/octet-stream':
                return f'Cannot compare binary files' # only compare text files
            other_path = other.path
        else:
            other_path = other

        a = self.path.read_text().splitlines(keepends=True)
        b = other_path.read_text().splitlines(keepends=True)

        diffresult = ''.join(difflib.unified_diff(a, b, fromfile=self.path.name, tofile=other_path.name))

        return diffresult

    def compare(self, other: FileArtifact | Path) -> bool:
        return self.diff(other) == ''

@dataclass
class FileArtifactDict(ArtifactDict):
    
    data: dict[str, FileArtifact] = field(default_factory=dict)


@dataclass
class FileArtifactList(ArtifactList):

    data: list[FileArtifact] = field(default_factory=list)

    def all_exist(self) -> bool:
        """
        Check if all artifact files in the caller exist.
        """
        return all(artifact.exists() for artifact in self.data)

    def paths_to_list(self) -> list[Path]:
        """
        Convert the list of artifact files in the caller to a list of their paths.
        """
        return [artifact.path for artifact in self.data]

    def unique_paths(self) -> FileArtifactList:
        """
        Reduce the list so that all paths are unique.
        """
        unique = FileArtifactList()
        for artifact in self.data:
            if artifact.path not in [a.path for a in unique.data]:
                unique.append(artifact)
            else:
                if artifact.pytestable:
                    # identify the matching artifact already in the unique list
                    matching = [a for a in unique.data if a.path == artifact.path]
                    if matching:
                        # if a matching artifact is found, keep the one with pytestable set
                        if artifact.pytestable:
                            unique.remove(matching[0])
                            unique.append(artifact)
                            logger.debug(f"Uniquify: {matching[0]} -> {artifact}")
        return unique

    def make_tarball(self, basename: str, remove: bool = False, arcname_prefix: str = None, unique: bool = False):
        """
        Create a tarball from the list of artifact files.
        
        Parameters
        ----------
        basename : str
            The base name for the tarball file.
        remove : bool
            Whether to remove the original files after creating the tarball; default False.
        """
        if remove:
            remove_files = []
            remove_artifacts = []
        process_list = self.data
        if unique:
            new_artifactlist: FileArtifactList = self.unique_paths()
            process_list = new_artifactlist.data
        with tarfile.open(f"{basename}.tar.gz", "w:gz") as tar:
            for artifact in process_list:
                if artifact.exists():
                    if arcname_prefix is None:
                        arcname = artifact.path.name
                    else:
                        arcname = os.path.join(arcname_prefix, artifact.path.name)
                    tar.add(artifact.path, arcname=arcname)
                    if remove:
                        remove_files.append(artifact.path)
                        remove_artifacts.append(artifact)
        if remove:
            logger.debug(f"Removing files:")
            my_logger([str(f) for f in remove_files], logger.debug, depth=1)
            for f in remove_files:
                try:
                    os.remove(f)
                except Exception as e:
                    logger.warning(f"Failed to remove {f}: {e}")
            for a in remove_artifacts:
                self.remove(a)

@dataclass
class StateArtifacts(FileArtifactDict):
    """
    Compound artifact holding congruent PSF/COOR/PDB/VEL files and XSC files 
    corresponding to the same state.
    """
    key: str = 'state'
    """ The key identifying the state artifacts; 'state' by default. """
    description: str = "Set of congruent topology/coordinate/system files"
    """ A brief description of the state artifacts. """
    pdb: str | Path | PDBFileArtifact | None = None
    """ The PDB file name, Path, or PDBFileArtifact. """
    psf: str | Path | PSFFileArtifact | None = None
    """ The PSF file name, Path, or PSFFileArtifact. """
    coor: str | Path | NAMDCoorFileArtifact | None = None
    """ The coordinate file name, Path, or NAMDCoorFileArtifact. """
    vel: str | Path | NAMDVelFileArtifact | None = None
    """ The velocity file name, Path, or NAMDVelFileArtifact. """
    xsc: str | Path | NAMDXscFileArtifact | None = None
    """ The XSC file name, Path, or NAMDXscFileArtifact. """

    def to_list(self):
        return [getattr(self, attr) for attr in ['pdb', 'psf', 'coor', 'vel', 'xsc'] if getattr(self, attr) is not None]

    def stamp(self, owner: object) -> StateArtifacts:
        self.produced_by = owner
        for attr in ['pdb', 'psf', 'coor', 'vel', 'xsc']:
            if artifact := getattr(self, attr):
                artifact.stamp(owner)
        return self

    def __post_init__(self):
        for attr, artifact_type in [('pdb', PDBFileArtifact), ('psf', PSFFileArtifact), ('coor', NAMDCoorFileArtifact), ('vel', NAMDVelFileArtifact), ('xsc', NAMDXscFileArtifact)]:
            my_attr = getattr(self, attr)
            if my_attr is None:
                my_attr = self.data.get(attr)
                if not my_attr:  # only set to None if not already in data
                    setattr(self, attr, None)
                    continue
                if isinstance(my_attr, artifact_type):
                    setattr(self, attr, my_attr)
                    self.data[attr] = my_attr
                elif isinstance(my_attr, str | Path):
                    setattr(self, attr, artifact_type(data=str(my_attr), key=attr))
                    self.data[attr] = getattr(self, attr)
                else:
                    raise TypeError(f"Invalid type for {attr}: {type(my_attr)}; expected str, Path, or {artifact_type}")
            else:
                if isinstance(my_attr, str | Path):
                    setattr(self, attr, artifact_type(data=str(my_attr), key=attr))
                    self.data[attr] = getattr(self, attr)
                elif isinstance(my_attr, artifact_type):
                    self.data[attr] = my_attr
                else:
                    raise TypeError(f"Invalid type for {attr}: {type(my_attr)}; expected str, Path, or {artifact_type}")

@dataclass
class TXTFileArtifact(FileArtifact):
    """
    A text file artifact.
    """
    description: str = "Text file"
    ext: str = 'txt'
    mime_type: str = 'text/plain'

@dataclass
class CharmmffFileArtifact(FileArtifact):
    """
    A CHARMM force field file artifact.
    """
    description: str = "CHARMM force field file"

@dataclass
class CharmmffFileArtifactList(FileArtifactList):
    """
    A collection of CHARMM file artifacts.
    """
    description: str = "CHARMM file artifacts"
    key: str = 'charmm_files'

    def append(self, item: CharmmffFileArtifact) -> None:
        """
        Append a CHARMM file artifact to the collection.
        """
        if not item.exists() or not isinstance(item, CharmmffFileArtifact):
            raise TypeError(f"Expected CharmmffFileArtifact, got {type(item)}")
        super().append(item)

    def extend(self, items: list[CharmmffFileArtifact]) -> None:
        """
        Extend the collection by appending elements from the iterable.
        """
        for item in items:
            self.append(item)

@dataclass
class CharmmffTopFileArtifact(CharmmffFileArtifact):
    """
    A toplevel CHARMM force field topology file artifact.
    """
    description: str = "CHARMM force field topology file"
    ext: str = 'rtf'
    
@dataclass
class CharmmffParFileArtifact(CharmmffFileArtifact):
    """
    A toplevel CHARMM force field parameter file artifact.
    """
    description: str = "CHARMM force field parameter file"
    ext: str = 'prm'

@dataclass
class CharmmffStreamFileArtifact(CharmmffFileArtifact):
    """
    A CHARMM force field stream file artifact.
    """
    description: str = "CHARMM force field stream file"
    ext: str = 'str'
    # stream: str = ''

@dataclass
class CharmmffTopFileArtifacts(CharmmffFileArtifactList):
    """
    A collection of CHARMM force field topology file artifacts.
    """
    description: str = "CHARMM force field topology files"
    key: str = 'charmmff_topfiles'

    def append(self, item: CharmmffTopFileArtifact | str | Path) -> None:
        """
        Append a CHARMM force field topology file artifact to the collection.
        """
        if isinstance(item, (str, Path)):
            item = CharmmffTopFileArtifact(data=str(item))
        if not isinstance(item, CharmmffTopFileArtifact):
            raise TypeError(f"Expected CharmmffTopFileArtifact, got {type(item)}")
        super().append(item)

@dataclass
class CharmmffParFileArtifacts(CharmmffFileArtifactList):
    """
    A collection of CHARMM force field parameter file artifacts.
    """
    description: str = "CHARMM force field parameter files"

    def append(self, item: CharmmffParFileArtifact | str | Path) -> None:
        """
        Append a CHARMM force field parameter file artifact to the collection.
        """
        if isinstance(item, (str, Path)):
            item = CharmmffParFileArtifact(data=str(item))
        if not isinstance(item, CharmmffParFileArtifact):
            raise TypeError(f"Expected CharmmffParFileArtifact, got {type(item)}")
        super().append(item)

@dataclass
class CharmmffStreamFileArtifacts(CharmmffFileArtifactList):
    """
    A collection of CHARMM force field stream file artifacts.
    """
    description: str = "CHARMM force field stream files"

    def append(self, item: CharmmffStreamFileArtifact | str | Path) -> None:
        """
        Append a CHARMM force field stream file artifact to the collection.
        """
        if isinstance(item, (str, Path)):
            item = CharmmffStreamFileArtifact(data=str(item))
        if not isinstance(item, CharmmffStreamFileArtifact):
            raise TypeError(f"Expected CharmmffStreamFileArtifact, got {type(item)}")
        super().append(item)

@dataclass
class YAMLFileArtifact(TXTFileArtifact):
    """
    A YAML file artifact.
    """
    description: str = "YAML file"
    ext: str = 'yaml'
    mime_type: str = 'application/x-yaml'

@dataclass
class JSONFileArtifact(TXTFileArtifact):
    """
    A JSON file artifact.
    """
    description: str = "JSON file"
    ext: str = 'json'
    mime_type: str = 'application/json'

@dataclass
class TclScriptArtifact(TXTFileArtifact):
    """
    A Tcl script file artifact.
    """
    description: str = "Tcl script file"
    ext: str = 'tcl'
    mime_type: str = 'application/x-tcl'
    
@dataclass
class LogFileArtifact(TXTFileArtifact):
    """
    A generic log file artifact.
    """
    description: str = "Log file generated during task execution"
    ext: str = 'log'

@dataclass
class LogFileArtifactList(FileArtifactList):
    """
    A list of log file artifacts.
    """
    description: str = "List of log files"

@dataclass
class PackmolLogFileArtifact(LogFileArtifact):
    """
    A Packmol log file artifact.
    """
    description: str = "Log file for Packmol execution"

@dataclass
class NAMDLogFileArtifact(LogFileArtifact):
    """
    A NAMD log file artifact.
    """
    description: str = "Log file for NAMD execution"

@dataclass
class NAMDLogFileArtifactList(LogFileArtifactList):
    """
    A list of NAMD log file artifacts.
    """
    description: str = "List of NAMD log files"

@dataclass
class VMDLogFileArtifact(LogFileArtifact):
    """
    A VMD log file artifact.
    """
    description: str = "Log file for VMD execution"

@dataclass
class PsfgenLogFileArtifact(VMDLogFileArtifact):
    """
    A psfgen log file artifact.
    """
    description: str = "Log file for psfgen execution"

@dataclass
class NAMDOutputFileArtifact(FileArtifact):
    """
    A generic binary NAMD output file artifact.
    """
    description: str = "Output file for NAMD execution"
    mime_type: str = 'application/octet-stream'

@dataclass
class NAMDOutputFileArtifactList(FileArtifactList):
    """
    A list of generic binary NAMD output file artifacts.
    """
    description: str = "List of output files for NAMD execution"

@dataclass
class NAMDCoorFileArtifact(NAMDOutputFileArtifact):
    """
    A NAMD binary coordinate file artifact.
    """
    description: str = "Binary coordinate file"
    ext: str = 'coor'

@dataclass
class NAMDCoorFileArtifactList(NAMDOutputFileArtifactList):
    """
    A list of NAMD binary coordinate file artifacts.
    """
    description: str = "List of NAMD binary coordinate files"

@dataclass
class NAMDVelFileArtifact(NAMDOutputFileArtifact):
    """
    A NAMD binary velocity file artifact.
    """
    description: str = "Binary velocity file"
    ext: str = 'vel'

@dataclass
class NAMDXscFileArtifact(LogFileArtifact):
    """
    A NAMD XSC file artifact.
    """
    description: str = "XSC file"
    ext: str = 'xsc'

@dataclass
class NAMDXstFileArtifact(LogFileArtifact):
    """
    A NAMD XST file artifact.
    """
    description: str = "XST file"
    ext: str = 'xst'

@dataclass
class NAMDDcdFileArtifact(NAMDOutputFileArtifact):
    """
    A NAMD binary DCD file artifact.
    """
    description: str = "Binary DCD file"
    ext: str = 'dcd'

@dataclass
class PDBFileArtifact(TXTFileArtifact):
    """
    A generic PDB file artifact.
    """
    description: str = "PDB file"
    ext: str = 'pdb'
    
    def compare(self, other: PDBFileArtifact | Path):
        """ A specific comparison method for pairs of PDB files.  Name-by-name
        congruency of atom list is checked, excluding all waters.
        """
        if not isinstance(other, PDBFileArtifact | Path) or not self.pytestable:
            raise TypeError(f"Expected PDBFileArtifact or Path, got {type(other)}")
        
        # Compare the two PDB file artifacts
        my_struct = PDBParser(filepath=self.name).parse().parsed
        other_struct = PDBParser(filepath=other.name).parse().parsed
        my_atoms = AtomList.from_pdb(my_struct)
        other_atoms = AtomList.from_pdb(other_struct)
        # check congrency by atoms, allow for variance in number of waters
        my_non_waters = filter(lambda x: Labels.segtype_of_resname[x.resname] != 'water', my_atoms)
        other_non_waters = filter(lambda x: Labels.segtype_of_resname[x.resname] != 'water', other_atoms)
        by_names = all(x.name==y.name for x,y in zip(my_non_waters, other_non_waters))
        return by_names

@dataclass
class PDBFileArtifactList(FileArtifactList):
    """
    A list of PDB file artifacts.
    """
    description: str = "List of PDB files"

@dataclass
class PackMolPDBForcedFileArtifact(PDBFileArtifact):
    """
    A Packmol PDB file artifact with forced formatting.
    """
    description: str = "Packmol-generated PDB file with forced formatting"
    ext: str = 'pdb_FORCED'

@dataclass
class CIFFileArtifact(TXTFileArtifact):
    """
    A CIF file artifact.
    """
    description: str = "CIF file"
    ext: str = 'cif'

@dataclass
class PQRFileArtifact(TXTFileArtifact):
    """
    A PQR file artifact (output generated by pdb2pqr).
    """
    description: str = "PQR file"
    ext: str = 'pqr'

@dataclass
class PSFFileArtifact(TXTFileArtifact):
    """
    A PSF file artifact.
    """
    description: str = "PSF file"
    ext: str = 'psf'

    def compare(self, other: PSFFileArtifact | Path) -> bool:
        """
        A specific comparison method for pairs of PSF files.  Waters are ignored.
        """
        if not isinstance(other, PSFFileArtifact | Path):
            raise TypeError(f"Expected PSFFileArtifact, got {type(other)}")
        # Compare the two PSF file artifacts
        my_struct = PSFContents(self.name)
        other_struct = PSFContents(other.name)
        my_non_waters = filter(lambda x: x.segtype != 'water', my_struct.atoms)
        other_non_waters = filter(lambda x: x.segtype != 'water', other_struct.atoms)
        return all(x.atomname == y.atomname for x, y in zip(my_non_waters, other_non_waters))

@dataclass
class PSFFileArtifactList(FileArtifactList):
    """
    A list of PSF file artifacts.
    """
    description: str = "List of PSF files"

@dataclass
class VMDScriptArtifact(TclScriptArtifact):
    """
    A VMD script file artifact.
    """
    description: str = "VMD script file"

@dataclass
class PsfgenInputScriptArtifact(VMDScriptArtifact):
    """
    A PSFgen input script artifact.
    """
    description: str = "PSFgen input script"

@dataclass
class NAMDConfigFileArtifact(TclScriptArtifact):
    """
    A NAMD configuration file artifact.
    """
    description: str = "NAMD configuration file"
    ext: str = 'namd'

    def compare(self, other: NAMDConfigFileArtifact | Path) -> bool:
        """
        A specific comparison method for NAMD Config files created by Pestifer.
        Creation time-stamps are ignored.
        """
        patchlist = PatchSet(self.diff(other).splitlines())
        if len(patchlist) == 0:
            return True
        if len(patchlist) < 2:
            pfile = patchlist[0]
            if len(pfile) == 0:
                return True
            if len(pfile) == 1:
                hunk = pfile[0]
                if len(hunk) == 0:
                    return True
                adds=[]
                dels=[]
                for line in hunk:
                    line=str(line)
                    if line.startswith('+') and not line.startswith('+++'):
                        adds.append(line[1:])
                    elif line.startswith('-') and not line.startswith('---'):
                        dels.append(line[1:])
                if len(adds) == 0 and len(dels) == 0:
                    return True
                if len(adds) > 0 and len(dels) > 0:
                    if len(adds) == 1 and 'Created' in adds[0] and len(dels) == 1 and 'Created' in dels[0]:
                        return True
        return False

@dataclass
class NAMDColvarsConfigArtifact(JSONFileArtifact):
    """
    A NAMD Colvars configuration file artifact.
    """
    description: str = "NAMD Colvars configuration file"
    ext: str = 'in'

@dataclass
class PackmolInputScriptArtifact(TXTFileArtifact):
    """
    A Packmol input script artifact.
    """
    description: str = "Packmol input script"
    ext: str = 'inp'
    

@dataclass
class NAMDColvarsTrajectoryArtifact(LogFileArtifact):
    """
    A NAMD Colvars trajectory output file artifact.
    """
    description: str = "NAMD Colvars trajectory output file"
    ext: str = 'colvars.traj'

@dataclass
class NAMDColvarsStateArtifact(LogFileArtifact):
    """
    A NAMD Colvars state output file artifact.
    """
    description: str = "NAMD Colvars state output file"
    ext: str = 'colvars.state'

@dataclass
class InputFileArtifact(TXTFileArtifact):
    """
    A generic input file artifact.
    """
    description: str = "Generic input file"
    ext: str = 'inp'
    
@dataclass
class DataFileArtifact(TXTFileArtifact):
    """
    A generic data file artifact.
    """
    description: str = "Data file"
    ext: str = 'dat'

@dataclass
class DataFileArtifactList(FileArtifactList):
    """ List of data file artifacts. """
    description: str = "List of data files"

@dataclass
class CSVDataFileArtifact(DataFileArtifact):
    """
    A CSV data file artifact.
    """
    description: str = "CSV data file"
    ext: str = 'csv'

@dataclass
class CSVDataFileArtifactList(DataFileArtifactList):
    """
    A list of CSV data file artifacts.
    """
    description: str = "List of CSV data files"

@dataclass
class PNGImageFileArtifact(FileArtifact):
    """
    A PNG image file artifact.
    """
    description: str = "PNG image file"
    ext: str = 'png'
    mime_type: str = 'image/png'
    nonstate_results: bool = True

@dataclass
class PNGImageFileArtifactList(FileArtifactList):
    """
    A list of PNG image file artifacts.
    """
    description: str = "List of PNG image files"
