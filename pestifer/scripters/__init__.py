from .filewriter import Filewriter
from .namdscripter import NAMDScripter
from .packmolscripter import PackmolScripter
from .psfgenscripter import PsfgenScripter
from .tclscripters import TcLScripter, VMDScripter

scripters = {
    'psfgen': PsfgenScripter,
    'vmd': VMDScripter,
    'tcl': TcLScripter,
    'packmol': PackmolScripter,
    'data': Filewriter,
    'namd': NAMDScripter
}
