# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
A class for handling artifacts in the Pestifer core. Artifacts are files or 
data generated during the execution of tasks and are managed by the 
:class:`~pestifer.core.pipeline.PipelineContext`.
"""
from __future__ import annotations
from abc import ABC, abstractmethod
from collections import UserList, UserDict
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
    """ The task or process that generated the artifact. """    
    provenance: list[object] = field(default_factory=list)
    """ A list of objects that contributed to the creation of the artifact, ordered chronologically. """

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
            raise ValueError("Owner must not be None.")
    
        if self.has_stamp():
            if self.produced_by is owner:
                logger.debug(f"Artifact {self.key} is already stamped by {repr(self.produced_by)}. No need to override.")
                return self
            logger.debug(f"Artifact {self.key} is already stamped by {repr(self.produced_by)}. Overriding with new owner {repr(owner)}.")
        self.provenance.append(self.produced_by)
        self.produced_by = owner
        for attr_name, attr_value in self.__dict__.items():
            if isinstance(attr_value, Artifact):
                logger.debug(f"Stamping member Artifact {attr_name} {type(attr_value)} of artifact {self.key} ")
                attr_value.stamp(owner)
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
    A list of Artifacts.
    """
    data: list[Artifact] = field(default_factory=list)

    @classmethod
    def from_dict(cls, artifact_dict: ArtifactDict) -> ArtifactList:
        """
        Create an ArtifactList from an ArtifactDict.
        """
        artifact_list = cls()
        for artifact in artifact_dict.values():
            artifact_list.append(artifact)
        artifact_list.produced_by = artifact_dict.produced_by
        artifact_list.key = artifact_dict.key
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
        if owner is None:
            raise ValueError("Owner must not be None.")
        self.produced_by = owner
        for artifact in self.data:
            if not isinstance(artifact, Artifact):
                raise TypeError(f'Expected Artifact, got {type(artifact)}')
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
    Dictionary of Artifacts.
    """
    data: dict[str, Artifact] = field(default_factory=dict)

    @classmethod
    def from_list(cls, artifact_list: ArtifactList) -> ArtifactDict:
        """
        Create an ArtifactDict from an ArtifactList.
        """
        artifact_dict = cls()
        for artifact in artifact_list:
            artifact_dict[artifact.key] = artifact
        artifact_dict.produced_by = artifact_list.produced_by
        artifact_dict.key = artifact_list.key
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
        if owner is None:
            raise ValueError("Owner must not be None.")
        self.produced_by = owner
        for artifact in self.data.values():
            if not isinstance(artifact, Artifact):
                raise TypeError(f'Expected Artifact, got {type(artifact)}')
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

    def exists(self) -> bool:
        """
        Check if the file artifact exists.
        """
        return self.path.exists()

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

        if isinstance(other, Path):
            other_path = other
        else:
            other_path = other.path

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

    def append(self, item: FileArtifact) -> None:
        """
        Append a FileArtifact to the list.

        Raises
        ------
        TypeError
            If the item is not a FileArtifact.
        """
        if not isinstance(item, FileArtifact):
            raise TypeError(f"Expected FileArtifact, got {type(item)}")
        super().append(item)

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

    def unique_paths(self) -> list[Path]:
        """
        Get a list of unique file paths from the artifact files.  If two matching artifacts are found, the one that is kept is the one with pytestable set, if any.
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
            logger.debug(f"Removing files: {remove_files}")
            for f in remove_files:
                try:
                    os.remove(f)
                except Exception as e:
                    logger.warning(f"Failed to remove {f}: {e}")
            for a in remove_artifacts:
                self.remove(a)

@dataclass
class StateArtifacts(Artifact):
    """
    Compound artifact holding congruent PSF/COOR/PDB/VEL files and XSC files 
    corresponding to the same state.
    """
    key: str = 'state'
    """ The key identifying the state artifacts; 'state' by default. """
    description: str = "Set of congruent topology/coordinate/system files"
    """ A brief description of the state artifacts. """
    psf: PSFFileArtifact | None = None
    """ The PSF file artifact. """
    pdb: PDBFileArtifact | None = None
    """ The PDB file artifact. """
    coor: NAMDCoorFileArtifact | None = None
    """ The coordinate file artifact. """
    vel: NAMDVelFileArtifact | None = None
    """ The velocity file artifact. """
    xsc: NAMDXscFileArtifact | None = None
    """ The XSC file artifact. """

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
    A text file artifact.
    """
    description: str = "Text file"
    ext: str = 'txt'
    mime_type: str = 'text/plain'

@dataclass
class CharmmffTopFileArtifact(TXTFileArtifact):
    """
    A toplevel CHARMM force field topology file artifact.
    """
    description: str = "CHARMM force field topology file"
    ext: str = 'rtf'
    
@dataclass
class CharmmffParFileArtifact(TXTFileArtifact):
    """
    A toplevel CHARMM force field parameter file artifact.
    """
    description: str = "CHARMM force field parameter file"
    ext: str = 'prm'

@dataclass
class CharmmffStreamFileArtifact(TXTFileArtifact):
    """
    A CHARMM force field stream file artifact.
    """
    description: str = "CHARMM force field stream file"
    ext: str = 'str'
    # stream: str = ''

@dataclass
class CharmmffTopFileArtifacts(FileArtifactList):
    """
    A collection of CHARMM force field topology file artifacts.
    """
    description: str = "CHARMM force field topology files"
    
    def append(self, item: CharmmffTopFileArtifact) -> None:
        """
        Append a CHARMM force field topology file artifact to the collection.
        """
        if not isinstance(item, CharmmffTopFileArtifact):
            raise TypeError(f"Expected CharmmffTopFileArtifact, got {type(item)}")
        super().append(item)

@dataclass
class CharmmffParFileArtifacts(FileArtifactList):
    """
    A collection of CHARMM force field parameter file artifacts.
    """
    description: str = "CHARMM force field parameter files"

    def append(self, item: CharmmffParFileArtifact) -> None:
        """
        Append a CHARMM force field parameter file artifact to the collection.
        """
        if not isinstance(item, CharmmffParFileArtifact):
            raise TypeError(f"Expected CharmmffParFileArtifact, got {type(item)}")
        super().append(item)

@dataclass
class CharmmffStreamFileArtifacts(FileArtifactList):
    """
    A collection of CHARMM force field stream file artifacts.
    """
    description: str = "CHARMM force field stream files"

    def append(self, item: CharmmffStreamFileArtifact) -> None:
        """
        Append a CHARMM force field stream file artifact to the collection.
        """
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
class NAMDCoorFileArtifact(NAMDOutputFileArtifact):
    """
    A NAMD binary coordinate file artifact.
    """
    description: str = "Binary coordinate file"
    ext: str = 'coor'

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
        if not isinstance(other, PDBFileArtifact | Path) or not self.pytestable:
            raise TypeError(f"Expected PDBFileArtifact, got {type(other)}")
        if isinstance(other, Path):
            other_name = other.name
        else:
            other_name = other.name
        # Compare the two PDB file artifacts
        my_struct = PDBParser(filepath=self.name).parse().parsed
        other_struct = PDBParser(filepath=other_name).parse().parsed
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

    def append(self, item: PDBFileArtifact) -> None:
        """
        Append a PDB file artifact to the list.
        """
        if not isinstance(item, PDBFileArtifact):
            raise TypeError(f"Expected PDBFileArtifact, got {type(item)}")
        super().append(item)

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
        if not isinstance(other, PSFFileArtifact | Path):
            raise TypeError(f"Expected PSFFileArtifact, got {type(other)}")
        # Compare the two PSF file artifacts
        my_struct = PSFContents(self.name)
        other_struct = PSFContents(other.name)
        my_non_waters = filter(lambda x: x.segtype != 'water', my_struct.atoms)
        other_non_waters = filter(lambda x: x.segtype != 'water', other_struct.atoms)
        return all(x.atomname == y.atomname for x, y in zip(my_non_waters, other_non_waters))

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
class CSVDataFileArtifact(DataFileArtifact):
    """
    A CSV data file artifact.
    """
    description: str = "CSV data file"
    ext: str = 'csv'

@dataclass
class PNGImageFileArtifact(FileArtifact):
    """
    A PNG image file artifact.
    """
    description: str = "PNG image file"
    ext: str = 'png'
    mime_type: str = 'image/png'