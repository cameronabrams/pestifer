from .filewriter import Filewriter
from .namdscripter import NAMDScripter
from .namdcolvarinputscripter import NAMDColvarInputScripter
from .packmolscripter import PackmolScripter
from .psfgenscripter import PsfgenScripter
from .tclscripters import TcLScripter, VMDScripter

scripters = {
    'psfgen': PsfgenScripter,
    'vmd': VMDScripter,
    'tcl': TcLScripter,
    'packmol': PackmolScripter,
    'data': Filewriter,
    'namd': NAMDScripter,
    'namd_colvar': NAMDColvarInputScripter
}
