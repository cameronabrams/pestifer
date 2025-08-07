# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
A class for handling artifacts in the Pestifer core. Artifacts are files or data generated during the execution of tasks.
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
    This class is used to store information about an artifact, including its key, value, type, and the task that produced it.
    The ``produced_by`` attribute indicates which task generated this artifact.
    The ``package`` attribute indicates whether the most recently generated artifact should be included in the product package, if it is a file.
    The ``data`` attribute can hold any data type, and the ``key`` attribute is used to identify the artifact.
    """
    data: object | None = None
    key: str | None = None
    produced_by: object | None = None
    package: bool = False

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

        Parameters
        ----------
        owner : Any, optional
            The owner of the artifact, which can be any object that produced this artifact. If not provided, the artifact will not be stamped.
        """
        if owner:
            if self.has_stamp():
                logger.warning(f"Artifact {self.key} is already stamped by {self.produced_by}. Overriding with new owner {owner}.")
            self.produced_by = owner
        else:
            raise ValueError("Owner must be provided for artifact stamping.")
        return self

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
    This class is used to store information about a file artifact, including its path and an optional description.
    """
    data: str | None = None  # base name without extension
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
            self.data = self.data.split('.')[0]

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
class TXTFileArtifact(FileArtifact):
    """
    Represents a text file artifact.
    """
    description: str = "Text file"
    ext: str = 'txt'
    package: bool = False
    mime_type: str = 'text/plain'

@dataclass
class DataArtifact(Artifact):
    """
    Represents a data artifact in the Pestifer core.
    This class is used to store information about a data artifact, including its value and the task that produced it.
    """
    description: str | None = None

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
        return all(artifact.exists() for artifact in self)

    def paths_to_list(self) -> list[Path]:
        """
        Convert the list of artifact files to a list of their paths.
        """
        return [artifact.path for artifact in self]

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
            for artifact in self:
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
class CharmmffTopFileArtifact(TXTFileArtifact):
    """
    Represents a toplevel CHARMM force field topology file.
    This class is used to store the topology file for CHARMM force fields that are used in any particular psfgen invocation.
    """
    description: str = "CHARMM force field topology file"
    ext: str = 'rtf'
    package: bool = False

@dataclass
class CharmmffParFileArtifact(TXTFileArtifact):
    """
    Represents a toplevel CHARMM force field parameter file.
    This class is used to store the parameter file for CHARMM force fields that are used in any particular psfgen invocation.
    """
    description: str = "CHARMM force field parameter file"
    ext: str = 'prm'
    package: bool = False

@dataclass
class CharmmffStreamFileArtifact(TXTFileArtifact):
    """
    Represents a CHARMM force field stream file.
    This class is used to store the stream file for CHARMM force fields that are used in any particular psfgen invocation.
    """
    description: str = "CHARMM force field stream file"
    ext: str = 'str'
    package: bool = False
    stream: str = ''

@dataclass
class CharmmffTopFileArtifacts(FileArtifactList):
    """
    Represents CHARMM force field topology files.
    This class is used to store the topology files for CHARMM force fields that are used in any particular psfgen invocation.
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
    This class is used to store the parameter files for CHARMM force fields that are used in any particular psfgen invocation.
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
    This class is used to store the stream files for CHARMM force fields that are used in any particular psfgen invocation.
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
    This class is used to store a YAML file that may be generated during the execution of tasks.
    The YAML file can contain various configurations or data in a structured format.
    """
    description: str = "YAML file"
    ext: str = 'yaml'
    package: bool = False

@dataclass
class JSONFileArtifact(TXTFileArtifact):
    """
    Represents a JSON file artifact.
    This class is used to store a JSON file that may be generated during the execution of tasks.
    The JSON file can contain structured data in JavaScript Object Notation format.
    """
    description: str = "JSON file"
    ext: str = 'json'
    package: bool = False

@dataclass
class TclScriptArtifact(TXTFileArtifact):
    """
    Represents a Tcl script file artifact.
    This class is used to store a Tcl script file that may be generated during the execution of tasks.
    The Tcl script file can contain various commands and configurations for molecular simulations.
    """
    description: str = "Tcl script file"
    ext: str = 'tcl'
    package: bool = False

@dataclass
class LogFileArtifact(TXTFileArtifact):
    """
    Represents a log file artifact.
    This class is used to store the log file generated during the execution of tasks.
    The log file contains information about the execution of tasks and any errors that may have occurred.
    """
    description: str = "Log file for task execution"
    ext: str = 'log'
    package: bool = False

@dataclass
class PackmolLogFileArtifact(LogFileArtifact):
    """
    Represents a Packmol log file artifact.
    This class is used to store the log file generated by Packmol during the execution of tasks.
    The Packmol log file contains information about the packing of molecules and any errors that may have occurred.
    """
    description: str = "Log file for Packmol execution"

@dataclass
class NAMDLogFileArtifact(LogFileArtifact):
    """
    Represents a NAMD log file artifact.
    This class is used to store the log file generated by NAMD during the execution of tasks.
    The NAMD log file contains information about the molecular dynamics simulation and any errors that may have occurred.
    """
    description: str = "Log file for NAMD execution"

class VMDLogFileArtifact(LogFileArtifact):
    """
    Represents a VMD log file artifact.
    This class is used to store the log file generated by VMD during the execution of tasks.
    The VMD log file contains information about the visualization of molecular structures and any errors that may have occurred.
    """
    description: str = "Log file for VMD execution"

@dataclass
class PsfgenLogFileArtifact(VMDLogFileArtifact):
    """
    Represents a psfgen log file artifact.
    This class is used to store the log file generated by psfgen during the execution of tasks.
    The psfgen log file contains information about the generation of PSF files and any errors that may have occurred.
    """
    description: str = "Log file for psfgen execution"

@dataclass
class NAMDOutputFileArtifact(FileArtifact):
    """
    Represents a NAMD output file artifact.
    This class is used to store the output file generated by NAMD during the execution of tasks.
    The NAMD output file contains information about the results of the molecular dynamics simulation.
    """
    description: str = "Output file for NAMD execution"
    mime_type: str = 'application/octet-stream'

@dataclass
class NAMDCoorFileArtifact(NAMDOutputFileArtifact):
    """
    Represents a NAMD coordinate file artifact.
    This class is used to store the coordinate file generated by NAMD during the execution of tasks.
    The NAMD coordinate file contains the atomic coordinates after the simulation.
    """
    description: str = "Binary coordinate file"
    ext: str = 'coor'
    package: bool = True

@dataclass
class NAMDVelFileArtifact(NAMDOutputFileArtifact):
    """
    Represents a NAMD velocity file artifact.
    This class is used to store the velocity file generated by NAMD during the execution of tasks.
    The NAMD velocity file contains the atomic velocities after the simulation.
    """
    description: str = "Binary velocity file"
    ext: str = 'vel'
    package: bool = True

@dataclass
class NAMDXscFileArtifact(LogFileArtifact):
    """
    Represents a NAMD XSC file artifact.
    This class is used to store the XSC file generated by NAMD during the execution of tasks.
    The NAMD XSC file contains the simulation box dimensions and periodic boundary conditions.
    """
    description: str = "XSC file"
    ext: str = 'xsc'
    package: bool = True

@dataclass
class NAMDXstFileArtifact(LogFileArtifact):
    """
    Represents a NAMD XST file artifact.
    This class is used to store the XST file generated by NAMD during the execution of tasks.
    The NAMD XST file contains the extended system information, including coordinates and velocities.
    """
    description: str = "XST file"
    ext: str = 'xst'
    package: bool = False
    mime_type: str = 'text/plain'

@dataclass
class NAMDDcdFileArtifact(NAMDOutputFileArtifact):
    """
    Represents a NAMD DCD file artifact.
    This class is used to store the DCD file generated by NAMD during the execution of tasks.
    The NAMD DCD file contains the trajectory information of the molecular dynamics simulation.
    """
    description: str = "Binary DCD file"
    ext: str = 'dcd'
    package: bool = False

@dataclass
class PDBFileArtifact(TXTFileArtifact):
    """
    Represents a base PDB file artifact.
    This class is used to store any PDB file generated during the execution of tasks.
    The base PDB file contains the initial atomic coordinates before any simulation.
    """
    description: str = "PDB file"
    ext: str = 'pdb'
    package: bool = True

@dataclass
class PDBFileArtifactList(FileArtifactList):
    """
    Represents a list of PDB files.
    This class is used to store a collection of PDB files, allowing for easy management and retrieval.
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
    This class is used to store a CIF file that may be fetched by a ``pestifer.tasks.fetch.FetchTask``.
    """
    description: str = "CIF file"
    ext: str = 'cif'
    package: bool = False

@dataclass
class PQRFileArtifact(TXTFileArtifact):
    """
    Represents a PQR file artifact.
    This class is used to store a PQR file that may be generated during the execution of ``pestifer.tasks.pdb2pqr.PDB2PQRTask``.
    """
    description: str = "PQR file"
    ext: str = 'pqr'
    package: bool = False

@dataclass
class PSFFileArtifact(TXTFileArtifact):
    """
    Represents a PSF file artifact.
    This class is used to store the PSF file generated during the execution of tasks.
    The PSF file contains the topology information of the molecular system.
    """
    description: str = "PSF file"
    ext: str = 'psf'
    package: bool = True

@dataclass
class VMDScriptArtifact(TclScriptArtifact):
    """
    Represents a VMD script file artifact.
    This class is used to store the VMD script file that may be generated during the execution of tasks.
    The VMD script file contains commands for visualizing molecular structures in VMD.
    """
    description: str = "VMD script file"
    package: bool = False

@dataclass
class PsfgenInputScriptArtifact(VMDScriptArtifact):
    """
    Represents a PSFgen input script artifact.
    This class is used to store the input script for PSFgen, which contains the specifications for generating PSF files.
    The PSFgen input script is typically a text file with commands and parameters for the PSFgen tool.
    """
    description: str = "PSFgen input script"
    package: bool = False

@dataclass
class NAMDConfigFileArtifact(TclScriptArtifact):
    """
    Represents a NAMD configuration file artifact.
    This class is used to store the configuration file for NAMD, which contains the parameters and settings for the molecular dynamics simulation.
    The NAMD configuration file is typically a text file with commands and parameters for the NAMD tool.
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
    This class is used to store the input script for Packmol, which contains the specifications for packing molecules into a defined space.
    The Packmol input script is typically a text file with commands and parameters for the Packmol tool.
    """
    description: str = "Packmol input script"
    ext: str = 'inp'
    package: bool = False

@dataclass
class NAMDColvarsOutputArtifact(LogFileArtifact):
    """
    Represents a NAMD Colvars output file artifact.
    This class is used to store the output file generated by NAMD Colvars, which contains the results of the collective variables during the molecular dynamics simulation.
    The NAMD Colvars output file is typically a text file with data on the collective variables.
    """
    description: str = "NAMD Colvars output file"
    ext: str = 'colvars.traj'
    package: bool = False

@dataclass
class InputFileArtifact(TXTFileArtifact):
    """
    Represents a generic input file artifact.
    This class is used to store a data file that may be generated during the execution of tasks.
    The data file can contain various types of information, such as simulation parameters or results.
    """
    description: str = "Generic input file"
    ext: str = 'inp'
    package: bool = False
    
@dataclass
class DataFileArtifact(TXTFileArtifact):
    """
    Represents a data file artifact.
    This class is used to store a data file that may be generated during the execution of tasks.
    The data file can contain various types of information, such as simulation parameters or results.
    """
    description: str = "Data file"
    ext: str = 'dat'
    package: bool = False

@dataclass
class CSVDataFileArtifact(DataFileArtifact):
    """
    Represents a CSV data file artifact.
    This class is used to store a CSV file that may be generated during the execution of tasks.
    The CSV data file contains tabular data in comma-separated values format.
    """
    description: str = "CSV data file"
    ext: str = 'csv'
    package: bool = False

@dataclass
class PNGImageFileArtifact(FileArtifact):
    """
    Represents a PNG image file artifact.
    This class is used to store a PNG image file that may be generated during the execution of tasks.
    """
    description: str = "PNG image file"
    ext: str = 'png'
    package: bool = False
    mime_type: str = 'image/png'