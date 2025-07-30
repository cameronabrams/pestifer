# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
A class for handling artifacts in the Pestifer core. Artifacts are files or data generated during the execution of tasks.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
import logging
from typing import Optional, List, Any, Iterable, Iterator

logger=logging.getLogger(__name__)

@dataclass
class Artifact(ABC):
    """
    Represents an artifact in the Pestifer core.
    This class is used to store information about an artifact, including its key, value, type, and the task that produced it.
    The ``produced_by`` attribute indicates which task generated this artifact.
    The ``package`` attribute indicates whether the most recently generated artifact should be included in the product package, if it is a file.
    The ``value`` attribute can hold any data type, and the ``key`` attribute is used to identify the artifact.
    The ``describe`` method should be implemented by subclasses to provide a description of the artifact.
    """
    value: Any = None
    key: Optional[str] = None
    produced_by: object = None
    package: bool = False

    @abstractmethod
    def describe(self) -> str:
        pass

    def stamp(self, owner: Any = None) -> None:
        """
        Stamp the artifact with the owner information.
        This method is used to set the owner of the artifact, which can be useful for tracking who created or modified it.  If owner is None, do nothing.

        Parameters
        ----------
        owner : Any, optional
            The owner of the artifact, which can be any object that produced this artifact. If not provided, the artifact will not be stamped.
        """
        if owner:
            self.produced_by = owner
        return self

@dataclass
class ArtifactList:
    """
    Represents a list of artifacts in the Pestifer core.
    This class is used to store a collection of artifacts, allowing for easy management and retrieval.
    """
    value: List[Artifact] = field(default_factory=list)

    def append(self, item: Artifact) -> None:
        self.value.append(item)

    def __len__(self) -> int:
        return len(self.value)

    def extend(self, items: Iterable[Artifact]) -> None:
        for item in items:
            self.append(item)

    def __iter__(self) -> Iterator[Artifact]:
        return iter(self.value)
    
    def __add__(self, other: "ArtifactList") -> "ArtifactList":
        if type(self) is not type(other):
            raise TypeError(f"Cannot add {type(self)} to {type(other)}")

        new = type(self)()  # assumes default constructor works
        new.extend(self)
        new.extend(other)
        return new

    def __getitem__(self, index):
        return self.value[index]

    def stamp(self, owner: Any = None) -> None:
        """
        Stamp the artifact list with the owner information.
        This method is used to set the owner of the artifact, which can be useful for tracking who created or modified it.
        """
        if owner:
            self.produced_by = owner
        else:
            logger.warning("No owner provided for artifact stamping.")
        return self

@dataclass
class ArtifactFile(Artifact):
    """
    Represents a file artifact in the Pestifer core.
    This class is used to store information about a file artifact, including its path and an optional description.
    """
    description: Optional[str] = None
    mime_type: Optional[str] = 'application/octet-stream'

    @property
    @abstractmethod
    def ext(self) -> str:
        """Return a string representing this file's 'type key' fallback."""
        pass

    def __post_init__(self):
        if self.key is None:
            self.key = self.ext

    @property
    def name(self) -> str:
        """Return the name of the file without the extension."""
        return self.value+r'.'+self.ext

    @property
    def path(self) -> Path:
        return Path(self.name)

    def exists(self) -> bool:
        return self.path.exists()

    def validate(self):
        if not self.path.exists():
            raise FileNotFoundError(f"{self.path} not found")
    
    def describe(self) -> str:
        return f"ArtifactFile(key={self.key}, value={self.value}, ext={self.ext}, mime_type={self.mime_type}, description={self.description})"

@dataclass
class TXTFile(ArtifactFile):
    """
    Represents a text file artifact.
    This class is used to store a text file that may be generated during the execution of tasks.
    The text file can contain any textual data, such as logs or configuration settings.
    """
    description: str = "Text file"
    ext: str = 'txt'
    package: bool = False
    mime_type: str = 'text/plain'

@dataclass
class ArtifactData(Artifact):
    """
    Represents a data artifact in the Pestifer core.
    This class is used to store information about a data artifact, including its value and the task that produced it.
    """
    description: Optional[str] = None

    def describe(self) -> str:
        return f"ArtifactData(key={self.key}, value={self.value}, description={self.description})"

@dataclass
class ArtifactFileList(ArtifactList):
    value: List[ArtifactFile] = field(default_factory=list)

    def make_tarball(self, basename: str) -> None:
        """
        Create a tarball from the list of artifact files.
        
        Parameters
        ----------
        basename : str
            The base name for the tarball file.
        """
        import tarfile
        with tarfile.open(f"{basename}.tar.gz", "w:gz") as tar:
            for artifact in self.value:
                if artifact.exists():
                    tar.add(artifact.path, arcname=artifact.path.name)
                else:
                    logger.warning(f"Artifact {artifact.path} does not exist and cannot be added to tarball.")

@dataclass
class CharmmffTopFile(TXTFile):
    """
    Represents a toplevel CHARMM force field topology file.
    This class is used to store the topology file for CHARMM force fields that are used in any particular psfgen invocation.
    """
    description: str = "CHARMM force field topology file"
    ext: str = 'rtf'
    package: bool = False

@dataclass
class CharmmffParFile(TXTFile):
    """
    Represents a toplevel CHARMM force field parameter file.
    This class is used to store the parameter file for CHARMM force fields that are used in any particular psfgen invocation.
    """
    description: str = "CHARMM force field parameter file"
    ext: str = 'prm'
    package: bool = False

@dataclass
class CharmmffStreamFile(TXTFile):
    """
    Represents a CHARMM force field stream file.
    This class is used to store the stream file for CHARMM force fields that are used in any particular psfgen invocation.
    """
    description: str = "CHARMM force field stream file"
    ext: str = 'str'
    package: bool = False
    stream: str = ''

@dataclass
class CharmmffTopFiles(ArtifactFileList):
    """
    Represents CHARMM force field topology files.
    This class is used to store the topology files for CHARMM force fields that are used in any particular psfgen invocation.
    """
    description: str = "CHARMM force field topology files"
    
    def append(self, item: CharmmffTopFile) -> None:
        if not isinstance(item, CharmmffTopFile):
            raise TypeError(f"Expected CharmmffTopFile, got {type(item)}")
        self.value.append(item)

@dataclass
class CharmmffParFiles(ArtifactFileList):
    """
    Represents CHARMM force field parameter files.
    This class is used to store the parameter files for CHARMM force fields that are used in any particular psfgen invocation.
    """
    description: str = "CHARMM force field parameter files"

    def append(self, item: CharmmffParFile) -> None:
        if not isinstance(item, CharmmffParFile):
            raise TypeError(f"Expected CharmmffParFile, got {type(item)}")
        self.value.append(item)

@dataclass
class CharmmffStreamFiles(ArtifactFileList):
    """
    Represents CHARMM force field stream files.
    This class is used to store the stream files for CHARMM force fields that are used in any particular psfgen invocation.
    """
    description: str = "CHARMM force field stream files"

    def append(self, item: CharmmffStreamFile) -> None:
        if not isinstance(item, CharmmffStreamFile):
            raise TypeError(f"Expected CharmmffStreamFile, got {type(item)}")
        self.value.append(item)

@dataclass
class YAMLFile(TXTFile):
    """
    Represents a YAML file artifact.
    This class is used to store a YAML file that may be generated during the execution of tasks.
    The YAML file can contain various configurations or data in a structured format.
    """
    description: str = "YAML file"
    ext: str = 'yaml'
    package: bool = False

@dataclass
class JSONFile(TXTFile):
    """
    Represents a JSON file artifact.
    This class is used to store a JSON file that may be generated during the execution of tasks.
    The JSON file can contain structured data in JavaScript Object Notation format.
    """
    description: str = "JSON file"
    ext: str = 'json'
    package: bool = False

@dataclass
class TclScript(TXTFile):
    """
    Represents a Tcl script file artifact.
    This class is used to store a Tcl script file that may be generated during the execution of tasks.
    The Tcl script file can contain various commands and configurations for molecular simulations.
    """
    description: str = "Tcl script file"
    ext: str = 'tcl'
    package: bool = False

@dataclass
class LogFile(TXTFile):
    """
    Represents a log file artifact.
    This class is used to store the log file generated during the execution of tasks.
    The log file contains information about the execution of tasks and any errors that may have occurred.
    """
    description: str = "Log file for task execution"
    ext: str = 'log'
    package: bool = False

@dataclass
class PackmolLogFile(LogFile):
    """
    Represents a Packmol log file artifact.
    This class is used to store the log file generated by Packmol during the execution of tasks.
    The Packmol log file contains information about the packing of molecules and any errors that may have occurred.
    """
    description: str = "Log file for Packmol execution"

@dataclass
class NAMDLogFile(LogFile):
    """
    Represents a NAMD log file artifact.
    This class is used to store the log file generated by NAMD during the execution of tasks.
    The NAMD log file contains information about the molecular dynamics simulation and any errors that may have occurred.
    """
    description: str = "Log file for NAMD execution"

class VMDLogFile(LogFile):
    """
    Represents a VMD log file artifact.
    This class is used to store the log file generated by VMD during the execution of tasks.
    The VMD log file contains information about the visualization of molecular structures and any errors that may have occurred.
    """
    description: str = "Log file for VMD execution"

@dataclass
class PsfgenLogFile(VMDLogFile):
    """
    Represents a psfgen log file artifact.
    This class is used to store the log file generated by psfgen during the execution of tasks.
    The psfgen log file contains information about the generation of PSF files and any errors that may have occurred.
    """
    description: str = "Log file for psfgen execution"

@dataclass
class NAMDOutputFile(ArtifactFile):
    """
    Represents a NAMD output file artifact.
    This class is used to store the output file generated by NAMD during the execution of tasks.
    The NAMD output file contains information about the results of the molecular dynamics simulation.
    """
    description: str = "Output file for NAMD execution"

@dataclass
class NAMDCoorFile(NAMDOutputFile):
    """
    Represents a NAMD coordinate file artifact.
    This class is used to store the coordinate file generated by NAMD during the execution of tasks.
    The NAMD coordinate file contains the atomic coordinates after the simulation.
    """
    description: str = "Binary coordinate file"
    ext: str = 'coor'
    package: bool = True

@dataclass
class NAMDVelFile(NAMDOutputFile):
    """
    Represents a NAMD velocity file artifact.
    This class is used to store the velocity file generated by NAMD during the execution of tasks.
    The NAMD velocity file contains the atomic velocities after the simulation.
    """
    description: str = "Binary velocity file"
    ext: str = 'vel'
    package: bool = True

@dataclass
class NAMDXscFile(NAMDOutputFile):
    """
    Represents a NAMD XSC file artifact.
    This class is used to store the XSC file generated by NAMD during the execution of tasks.
    The NAMD XSC file contains the simulation box dimensions and periodic boundary conditions.
    """
    description: str = "XSC file"
    ext: str = 'xsc'
    package: bool = True

@dataclass
class NAMDXstFile(NAMDOutputFile):
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
class NAMDDcdFile(NAMDOutputFile):
    """
    Represents a NAMD DCD file artifact.
    This class is used to store the DCD file generated by NAMD during the execution of tasks.
    The NAMD DCD file contains the trajectory information of the molecular dynamics simulation.
    """
    description: str = "Binary DCD file"
    ext: str = 'dcd'
    package: bool = False

@dataclass
class PDBFile(TXTFile):
    """
    Represents a base PDB file artifact.
    This class is used to store any PDB file generated during the execution of tasks.
    The base PDB file contains the initial atomic coordinates before any simulation.
    """
    description: str = "PDB file"
    ext: str = 'pdb'
    package: bool = True

@dataclass
class PDBFileList(ArtifactFileList):
    """
    Represents a list of PDB files.
    This class is used to store a collection of PDB files, allowing for easy management and retrieval.
    """
    description: str = "List of PDB files"

    def append(self, item: PDBFile) -> None:
        if not isinstance(item, PDBFile):
            raise TypeError(f"Expected PDBFile, got {type(item)}")
        self.value.append(item)

@dataclass
class CIFFile(TXTFile):
    """
    Represents a CIF file artifact.
    This class is used to store a CIF file that may be fetched by a ``pestifer.tasks.fetch.FetchTask``.
    """
    description: str = "CIF file"
    ext: str = 'cif'
    package: bool = False

@dataclass
class PQRFile(TXTFile):
    """
    Represents a PQR file artifact.
    This class is used to store a PQR file that may be generated during the execution of ``pestifer.tasks.pdb2pqr.PDB2PQRTask``.
    """
    description: str = "PQR file"
    ext: str = 'pqr'
    package: bool = False

@dataclass
class PSFFile(TXTFile):
    """
    Represents a PSF file artifact.
    This class is used to store the PSF file generated during the execution of tasks.
    The PSF file contains the topology information of the molecular system.
    """
    description: str = "PSF file"
    ext: str = 'psf'
    package: bool = True

@dataclass
class VMDScript(TclScript):
    """
    Represents a VMD script file artifact.
    This class is used to store the VMD script file that may be generated during the execution of tasks.
    The VMD script file contains commands for visualizing molecular structures in VMD.
    """
    description: str = "VMD script file"
    package: bool = False

@dataclass
class PsfgenInputScript(VMDScript):
    """
    Represents a PSFgen input script artifact.
    This class is used to store the input script for PSFgen, which contains the specifications for generating PSF files.
    The PSFgen input script is typically a text file with commands and parameters for the PSFgen tool.
    """
    description: str = "PSFgen input script"
    package: bool = False

@dataclass
class NAMDConfigFile(TclScript):
    """
    Represents a NAMD configuration file artifact.
    This class is used to store the configuration file for NAMD, which contains the parameters and settings for the molecular dynamics simulation.
    The NAMD configuration file is typically a text file with commands and parameters for the NAMD tool.
    """
    description: str = "NAMD configuration file"
    ext: str = 'namd'
    package: bool = True

@dataclass
class NAMDColvarsConfig(JSONFile):
    """
    Represents a NAMD Colvars configuration file artifact.
    This class is used to store the configuration file for NAMD Colvars, which contains the parameters and settings for collective variables in the molecular dynamics simulation.
    The NAMD Colvars configuration file is typically a text file with commands and parameters for the NAMD Colvars tool.
    """
    description: str = "NAMD Colvars configuration file"
    ext: str = 'in'
    package: bool = False

@dataclass
class PackmolInputScript(TXTFile):
    """
    Represents a Packmol input script artifact.
    This class is used to store the input script for Packmol, which contains the specifications for packing molecules into a defined space.
    The Packmol input script is typically a text file with commands and parameters for the Packmol tool.
    """
    description: str = "Packmol input script"
    ext: str = 'inp'
    package: bool = False

@dataclass
class NAMDColvarsOutput(TXTFile):
    """
    Represents a NAMD Colvars output file artifact.
    This class is used to store the output file generated by NAMD Colvars, which contains the results of the collective variables during the molecular dynamics simulation.
    The NAMD Colvars output file is typically a text file with data on the collective variables.
    """
    description: str = "NAMD Colvars output file"
    ext: str = 'colvars.traj'
    package: bool = False

@dataclass
class InputFile(TXTFile):
    """
    Represents a generic input file artifact.
    This class is used to store a data file that may be generated during the execution of tasks.
    The data file can contain various types of information, such as simulation parameters or results.
    """
    description: str = "Generic input file"
    ext: str = 'inp'
    package: bool = False
    
@dataclass
class DataFile(TXTFile):
    """
    Represents a data file artifact.
    This class is used to store a data file that may be generated during the execution of tasks.
    The data file can contain various types of information, such as simulation parameters or results.
    """
    description: str = "Data file"
    ext: str = 'dat'
    package: bool = False

@dataclass
class CSVDataFile(TXTFile):
    """
    Represents a CSV data file artifact.
    This class is used to store a CSV file that may be generated during the execution of tasks.
    The CSV data file contains tabular data in comma-separated values format.
    """
    description: str = "CSV data file"
    ext: str = 'csv'
    package: bool = False

@dataclass
class PNGImageFile(ArtifactFile):
    """
    Represents a PNG image file artifact.
    This class is used to store a PNG image file that may be generated during the execution of tasks.
    """
    description: str = "PNG image file"
    ext: str = 'png'
    package: bool = False
    mime_type: str = 'image/png'