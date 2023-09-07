from argparse import Namespace
from .util import inspect_classes
from .mods import ModTypes
import logging
logger=logging.getLogger(__name__)
class ModContainer(Namespace):
    mod_classes,modlist_classes=inspect_classes('pestifer.mods','List')
    def __init__(self,input_specs={},**kwargs):
        prep_dict={f'{a}s':kwargs.get(a,Namespace()) for a in ModTypes}
        # logger.debug(f'modcontainer prepdict {prep_dict}')
        # prep_dict['modclasses']=mod_classes
        if input_specs:
            # we are parsing yaml input directly from user config file
            for name,Cls in self.mod_classes.items():
                modcat=f'{Cls.modtype}s'
                if Cls.yaml_header in input_specs:
                    LCls=self.modlist_classes.get(f'{name}List',list)
                    prep_dict[modcat].__dict__[Cls.yaml_header]=LCls([])
                    for entry in input_specs[Cls.yaml_header]:
                        assert type(entry) in [str,dict]
                        prep_dict[modcat].__dict__[Cls.yaml_header].append(Cls(entry))
        # already parsed mods can be incorporated
        for yaml_name,parsed_modlist in kwargs.items():
            modtype=[Cls.modtype for name,Cls in self.mod_classes.items() if Cls.yaml_header==yaml_name][0]
            prep_dict[modtype].__dict__[yaml_name]=parsed_modlist
        super().__init__(**prep_dict) 
        # logger.debug(f'modcontainer dict {self.__dict__}')
