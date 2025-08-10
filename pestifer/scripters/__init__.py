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
