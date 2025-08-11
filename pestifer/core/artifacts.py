# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
A class for handling artifacts in the Pestifer core. Artifacts are files or 
data generated during the execution of tasks.  They are managed by the pipeline context.
"""
from __future__ import annotations
from abc import ABC, abstractmethod
from collections import UserList, UserDict
from dataclasses import dataclass, field
from pathlib import Path
import logging
import os
import tarfile

logger = logging.getLogger(__name__)

@dataclass
class Artifact():
    """
    Represents an artifact in the Pestifer core.
    The ``produced_by`` attribute indicates which task generated this artifact.
    The ``package`` attribute indicates whether the most recently generated artifact should be included in the product package, if it is a file.
    The ``data`` attribute can hold any data type, and the ``key`` attribute is used to identify the artifact.
    """
    data: object | None = None
    key: str | None = None
    produced_by: object | None = None
    package: bool = False
    provenance: list[object] = field(default_factory=list)

    def has_stamp(self) -> bool:        
        """
        Check if the artifact has a stamp (owner information).
        
        Returns
        -------
        bool
            True if the artifact has a stamp, False otherwise.
        """
        return self.produced_by is not None

    def stamp(self, owner: object | None = None) -> Artifact:
        """
        Stamp the artifact with the owner information.
        This method is used to set the owner of the artifact, which can be useful for tracking who created or modified it.  If owner is None, do nothing.
        Any attributes that are instances of Artifact are also stamped.

        Parameters
        ----------
        owner : Any, optional
            The owner of the artifact, which can be any object that produced this artifact. If not provided, the artifact will not be stamped.
        """
        if owner:
            if self.has_stamp():
                if self.produced_by is owner:
                    logger.debug(f"Artifact {self.key} is already stamped by {repr(self.produced_by)}. No need to override.")
                    return self
                logger.debug(f"Artifact {self.key} is already stamped by {repr(self.produced_by)}. Overriding with new owner {repr(owner)}.")
            self.provenance.append(self.produced_by)
            self.produced_by = owner
            for attr_name, attr_value in self.__dict__.items():
                if isinstance(attr_value, Artifact):
                    attr_value.stamp(owner)
        else:
            raise ValueError("Owner must be provided for artifact stamping.")
        return self

@dataclass
class DataArtifact(Artifact):
    """
    Represents a data artifact in the Pestifer core.
    """
    description: str | None = None

@dataclass
class ArtifactList(UserList, Artifact):
    """
    Represents a list of artifacts in the Pestifer core.
    This class is used to store a collection of artifacts, allowing for easy management and retrieval.
    """
    data: list[Artifact] = field(default_factory=list)

    @classmethod
    def from_dict(cls, artifact_dict: ArtifactDict) -> ArtifactList:
        artifact_list = cls()
        for artifact in artifact_dict.values():
            artifact_list.append(artifact)
        artifact_list.produced_by = artifact_dict.produced_by
        artifact_list.key = artifact_dict.key
        return artifact_list
    
    def to_dict(self) -> ArtifactDict:
        return ArtifactDict.from_list(self)
    
    def stamp(self, owner: object = None) -> ArtifactList:
        """
        Stamp all artifacts in the list with the owner information.
        
        Parameters
        ----------
        owner : Any, optional
            The owner of the artifacts, which can be any object that produced these artifacts. If not provided, the artifacts will not be stamped.
        
        Returns
        -------
        ArtifactList
            The artifact list with all artifacts stamped with the owner information.
        """
        if owner:
            self.produced_by = owner
            for artifact in self.data:
                artifact.stamp(owner)
            return self
        raise ValueError("Owner must be provided for artifact stamping.")

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
        filtered_list = ArtifactList()
        filtered_list.key = self.key
        filtered_list.produced_by = produced_by
        filtered_list.data = [artifact for artifact in self.data if artifact.produced_by == produced_by]
        return filtered_list
    
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
        filtered_list = ArtifactList()
        filtered_list.key = self.key
        filtered_list.produced_by = self.produced_by
        filtered_list.data = [artifact for artifact in self.data if isinstance(artifact, artifact_type)]
        return filtered_list

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
        filtered_list = ArtifactList()
        filtered_list.key = self.key
        filtered_list.produced_by = self.produced_by
        filtered_list.data = [artifact for artifact in self.data if artifact.key == key]
        return filtered_list

@dataclass
class ArtifactDict(UserDict, Artifact):
    """
    Represents a dictionary of artifacts in the Pestifer core.
    This class is used to store a collection of artifacts with unique keys, allowing for easy management and retrieval.

    """
    data: dict[str, Artifact] = field(default_factory=dict)

    @classmethod
    def from_list(cls, artifact_list: ArtifactList) -> ArtifactDict:
        artifact_dict = cls()
        for artifact in artifact_list:
            artifact_dict[artifact.key] = artifact
        artifact_dict.produced_by = artifact_list.produced_by
        artifact_dict.key = artifact_list.key
        return artifact_dict
    
    def to_list(self) -> ArtifactList:
        return ArtifactList.from_dict(self)
    
    def stamp(self, owner: object = None) -> ArtifactDict:
        """
        Stamp all artifacts in the dictionary with the owner information.
        
        Parameters
        ----------
        owner : Any, optional
            The owner of the artifacts, which can be any object that produced these artifacts. If not provided, the artifacts will not be stamped.
        
        Returns
        -------
        ArtifactDict
            The artifact dictionary with all artifacts stamped with the owner information.
        """
        if owner:
            self.produced_by = owner
            for artifact in self.data.values():
                artifact.stamp(owner)
            return self
        raise ValueError("Owner must be provided for artifact stamping.")

    def update_item(self, artifact: Artifact, key: str | None = None) -> None:
        if key:
            self.data[key] = artifact
        else:
            self.data[artifact.key] = artifact

    def filter_by_produced_by(self, produced_by: object) -> ArtifactDict:
        filtered_dict = ArtifactDict()
        filtered_dict.key = self.key
        filtered_dict.produced_by = produced_by
        for key, artifact in self.data.items():
            if artifact.produced_by == produced_by:
                filtered_dict[key] = artifact
        return filtered_dict

@dataclass
class FileArtifact(Artifact, ABC):
    """
    Represents a file artifact in the Pestifer core.
    """
    data: str | None = None  # name with or without extension
    description: str | None = None
    mime_type: str | None = 'application/octet-stream'

    @property
    @abstractmethod
    def ext(self) -> str:
        """Return a string representing this file's 'type key' fallback."""
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
        return self.data + '.' + self.ext

    @property
    def path(self) -> Path:
        return Path(self.name)

    def exists(self) -> bool:
        return self.path.exists()

    def validate(self):
        if not self.exists():
            raise FileNotFoundError(f"{self.path} not found")
    
@dataclass
class FileArtifactDict(ArtifactDict):
    data: dict[str, FileArtifact] = field(default_factory=dict)

@dataclass
class FileArtifactList(ArtifactList):
    data: list[FileArtifact] = field(default_factory=list)

    def append(self, item: FileArtifact) -> None:
        if not isinstance(item, FileArtifact):
            raise TypeError(f"Expected FileArtifact, got {type(item)}")
        super().append(item)

    def all_exist(self) -> bool:
        """
        Check if all artifact files exist.
        """
        return all(artifact.exists() for artifact in self.data)

    def paths_to_list(self) -> list[Path]:
        """
        Convert the list of artifact files to a list of their paths.
        """
        return [artifact.path for artifact in self.data]

    def make_tarball(self, basename: str, remove: bool = False) -> None:
        """
        Create a tarball from the list of artifact files.
        
        Parameters
        ----------
        basename : str
            The base name for the tarball file.
        """
        if remove:
            remove_files = []
            remove_artifacts = []
        with tarfile.open(f"{basename}.tar.gz", "w:gz") as tar:
            for artifact in self.data:
                if artifact.exists():
                    tar.add(artifact.path, arcname=artifact.path.name)
                    if remove:
                        remove_files.append(artifact.path)
                        remove_artifacts.append(artifact)
        if remove:
            logger.debug(f"Removing files: {remove_files}")
            for f in remove_files:
                try:
                    os.remove(f)
                except Exception as e:
                    logger.warning(f"Failed to remove {f}: {e}")
            for a in remove_artifacts:
                self.remove(a)
        return

@dataclass
class StateArtifacts(Artifact):
    """
    Compound artifact holding congruent PSF/COOR/PDB/VEL files and XSC files 
    corresponding to the same state.
    """
    key: str = 'state'
    description: str = "Set of congruent topology/coordinate/system files"
    psf: PSFFileArtifact | None = None
    pdb: PDBFileArtifact | None = None
    coor: NAMDCoorFileArtifact | None = None
    vel: NAMDVelFileArtifact | None = None
    xsc: NAMDXscFileArtifact | None = None

    def to_list(self) -> list[FileArtifact]:
        """
        Convert the state artifacts to a list of file artifacts, excluding
        any None's.
        """
        return_me = []
        for artifact in [self.psf, self.pdb, self.coor, self.vel, self.xsc]:
            if artifact is not None:
                return_me.append(artifact)
        return return_me

@dataclass
class TXTFileArtifact(FileArtifact):
    """
    Represents a text file artifact.
    """
    description: str = "Text file"
    ext: str = 'txt'
    package: bool = False
    mime_type: str = 'text/plain'

@dataclass
class CharmmffTopFileArtifact(TXTFileArtifact):
    """
    Represents a toplevel CHARMM force field topology file.
    """
    description: str = "CHARMM force field topology file"
    ext: str = 'rtf'
    package: bool = False

@dataclass
class CharmmffParFileArtifact(TXTFileArtifact):
    """
    Represents a toplevel CHARMM force field parameter file.
    """
    description: str = "CHARMM force field parameter file"
    ext: str = 'prm'
    package: bool = False

@dataclass
class CharmmffStreamFileArtifact(TXTFileArtifact):
    """
    Represents a CHARMM force field stream file.
    """
    description: str = "CHARMM force field stream file"
    ext: str = 'str'
    package: bool = False
    stream: str = ''

@dataclass
class CharmmffTopFileArtifacts(FileArtifactList):
    """
    Represents CHARMM force field topology files.
    """
    description: str = "CHARMM force field topology files"
    
    def append(self, item: CharmmffTopFileArtifact) -> None:
        if not isinstance(item, CharmmffTopFileArtifact):
            raise TypeError(f"Expected CharmmffTopFileArtifact, got {type(item)}")
        super().append(item)

@dataclass
class CharmmffParFileArtifacts(FileArtifactList):
    """
    Represents CHARMM force field parameter files.
    """
    description: str = "CHARMM force field parameter files"

    def append(self, item: CharmmffParFileArtifact) -> None:
        if not isinstance(item, CharmmffParFileArtifact):
            raise TypeError(f"Expected CharmmffParFileArtifact, got {type(item)}")
        super().append(item)

@dataclass
class CharmmffStreamFileArtifacts(FileArtifactList):
    """
    Represents CHARMM force field stream files.
    """
    description: str = "CHARMM force field stream files"

    def append(self, item: CharmmffStreamFileArtifact) -> None:
        if not isinstance(item, CharmmffStreamFileArtifact):
            raise TypeError(f"Expected CharmmffStreamFileArtifact, got {type(item)}")
        super().append(item)

@dataclass
class YAMLFileArtifact(TXTFileArtifact):
    """
    Represents a YAML file artifact.
    """
    description: str = "YAML file"
    ext: str = 'yaml'
    package: bool = False

@dataclass
class JSONFileArtifact(TXTFileArtifact):
    """
    Represents a JSON file artifact.
    """
    description: str = "JSON file"
    ext: str = 'json'
    package: bool = False

@dataclass
class TclScriptArtifact(TXTFileArtifact):
    """
    Represents a Tcl script file artifact.
    """
    description: str = "Tcl script file"
    ext: str = 'tcl'
    package: bool = False

@dataclass
class LogFileArtifact(TXTFileArtifact):
    """
    Represents a log file artifact.
    """
    description: str = "Log file generated during task execution"
    ext: str = 'log'
    package: bool = False

@dataclass
class PackmolLogFileArtifact(LogFileArtifact):
    """
    Represents a Packmol log file artifact.
    """
    description: str = "Log file for Packmol execution"

@dataclass
class NAMDLogFileArtifact(LogFileArtifact):
    """
    Represents a NAMD log file artifact.
    """
    description: str = "Log file for NAMD execution"

class VMDLogFileArtifact(LogFileArtifact):
    """
    Represents a VMD log file artifact.
    """
    description: str = "Log file for VMD execution"

@dataclass
class PsfgenLogFileArtifact(VMDLogFileArtifact):
    """
    Represents a psfgen log file artifact.
    """
    description: str = "Log file for psfgen execution"

@dataclass
class NAMDOutputFileArtifact(FileArtifact):
    """
    Represents a NAMD output file artifact.
    """
    description: str = "Output file for NAMD execution"
    mime_type: str = 'application/octet-stream'

@dataclass
class NAMDCoorFileArtifact(NAMDOutputFileArtifact):
    """
    Represents a NAMD coordinate file artifact.
    """
    description: str = "Binary coordinate file"
    ext: str = 'coor'
    package: bool = True

@dataclass
class NAMDVelFileArtifact(NAMDOutputFileArtifact):
    """
    Represents a NAMD velocity file artifact.
    """
    description: str = "Binary velocity file"
    ext: str = 'vel'
    package: bool = True

@dataclass
class NAMDXscFileArtifact(LogFileArtifact):
    """
    Represents a NAMD XSC file artifact.
    """
    description: str = "XSC file"
    ext: str = 'xsc'
    package: bool = True

@dataclass
class NAMDXstFileArtifact(LogFileArtifact):
    """
    Represents a NAMD XST file artifact.
    """
    description: str = "XST file"
    ext: str = 'xst'
    package: bool = False
    mime_type: str = 'text/plain'

@dataclass
class NAMDDcdFileArtifact(NAMDOutputFileArtifact):
    """
    Represents a NAMD DCD file artifact.
    """
    description: str = "Binary DCD file"
    ext: str = 'dcd'
    package: bool = False

@dataclass
class PDBFileArtifact(TXTFileArtifact):
    """
    Represents a base PDB file artifact.
    """
    description: str = "PDB file"
    ext: str = 'pdb'
    package: bool = True

@dataclass
class PDBFileArtifactList(FileArtifactList):
    """
    Represents a list of PDB files.
    """
    description: str = "List of PDB files"

    def append(self, item: PDBFileArtifact) -> None:
        if not isinstance(item, PDBFileArtifact):
            raise TypeError(f"Expected PDBFileArtifact, got {type(item)}")
        super().append(item)

@dataclass
class CIFFileArtifact(TXTFileArtifact):
    """
    Represents a CIF file artifact.
    """
    description: str = "CIF file"
    ext: str = 'cif'
    package: bool = False

@dataclass
class PQRFileArtifact(TXTFileArtifact):
    """
    Represents a PQR file artifact.
    """
    description: str = "PQR file"
    ext: str = 'pqr'
    package: bool = False

@dataclass
class PSFFileArtifact(TXTFileArtifact):
    """
    Represents a PSF file artifact.
    """
    description: str = "PSF file"
    ext: str = 'psf'
    package: bool = True

@dataclass
class VMDScriptArtifact(TclScriptArtifact):
    """
    Represents a VMD script file artifact.
    """
    description: str = "VMD script file"
    package: bool = False

@dataclass
class PsfgenInputScriptArtifact(VMDScriptArtifact):
    """
    Represents a PSFgen input script artifact.
    """
    description: str = "PSFgen input script"
    package: bool = False

@dataclass
class NAMDConfigFileArtifact(TclScriptArtifact):
    """
    Represents a NAMD configuration file artifact.
    """
    description: str = "NAMD configuration file"
    ext: str = 'namd'
    package: bool = True

@dataclass
class NAMDColvarsConfigArtifact(JSONFileArtifact):
    """
    Represents a NAMD Colvars configuration file artifact.
    This class is used to store the configuration file for NAMD Colvars, which contains the parameters and settings for collective variables in the molecular dynamics simulation.
    The NAMD Colvars configuration file is typically a text file with commands and parameters for the NAMD Colvars tool.
    """
    description: str = "NAMD Colvars configuration file"
    ext: str = 'in'
    package: bool = False

@dataclass
class PackmolInputScriptArtifact(TXTFileArtifact):
    """
    Represents a Packmol input script artifact.
    """
    description: str = "Packmol input script"
    ext: str = 'inp'
    package: bool = False

@dataclass
class NAMDColvarsTrajectoryArtifact(LogFileArtifact):
    """
    Represents a NAMD Colvars trajectory output file artifact.
    """
    description: str = "NAMD Colvars trajectory output file"
    ext: str = 'colvars.traj'
    package: bool = False

@dataclass
class NAMDColvarsStateArtifact(LogFileArtifact):
    """
    Represents a NAMD Colvars state output file artifact.
    """
    description: str = "NAMD Colvars state output file"
    ext: str = 'colvars.state'
    package: bool = False

@dataclass
class InputFileArtifact(TXTFileArtifact):
    """
    Represents a generic input file artifact.
    """
    description: str = "Generic input file"
    ext: str = 'inp'
    package: bool = False
    
@dataclass
class DataFileArtifact(TXTFileArtifact):
    """
    Represents a data file artifact.
    """
    description: str = "Data file"
    ext: str = 'dat'
    package: bool = False

@dataclass
class CSVDataFileArtifact(DataFileArtifact):
    """
    Represents a CSV data file artifact.
    """
    description: str = "CSV data file"
    ext: str = 'csv'
    package: bool = False

@dataclass
class PNGImageFileArtifact(FileArtifact):
    """
    Represents a PNG image file artifact.
    """
    description: str = "PNG image file"
    ext: str = 'png'
    package: bool = False
    mime_type: str = 'image/png'