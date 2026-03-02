"""
This subpackage defines several scripters for writing different types of input files used
by pestifer. All scripters used in pestifer are subclasses of
:class:`pestifer.scripters.generic.GenericScripter`.

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

from .generic import GenericScripter
from .namd import NAMDScripter
from .namd_colvar_input import NAMDColvarInputScripter
from .packmol import PackmolScripter
from .psfgen import PsfgenScripter
from .tcl import TcLScripter
from .vmd import VMDScripter

scripters = {
    'psfgen': PsfgenScripter,
    'vmd': VMDScripter,
    'tcl': TcLScripter,
    'packmol': PackmolScripter,
    'data': GenericScripter,
    'namd': NAMDScripter,
    'namd_colvar': NAMDColvarInputScripter
}
