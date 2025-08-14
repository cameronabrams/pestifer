"""
This subpackage defines several scripters for writing different types of input files used
by pestifer. All scripters used in pestifer are subclasses of
:class:`pestifer.scripters.genericscripter.GenericScripter`.

.. mermaid::
  :caption: Scripter dependencies.
  
  graph TD;
    GenericScripter --> TcLScripter;
    GenericScripter --> NAMDColvarInputScripter;
    GenericScripter --> PackmolScripter;
    TcLScripter --> VMDScripter;
    TcLScripter --> NAMDScripter;
    VMDScripter --> PsfgenScripter;

"""

from .genericscripter import GenericScripter
from .namdscripter import NAMDScripter
from .namdcolvarinputscripter import NAMDColvarInputScripter
from .packmolscripter import PackmolScripter
from .psfgenscripter import PsfgenScripter
from .tclscripter import TcLScripter
from .vmdscripter import VMDScripter

scripters = {
    'psfgen': PsfgenScripter,
    'vmd': VMDScripter,
    'tcl': TcLScripter,
    'packmol': PackmolScripter,
    'data': GenericScripter,
    'namd': NAMDScripter,
    'namd_colvar': NAMDColvarInputScripter
}
