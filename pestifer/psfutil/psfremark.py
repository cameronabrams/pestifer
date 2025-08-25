# Author: Cameron F. Abrams <cfa22@drexel.edu>

import logging

from collections import UserList
from dataclasses import dataclass, field

from .psfpatch import PSFDISUPatch, PSFLinkPatch

logger = logging.getLogger(__name__)

from ..objs.link import Link

@dataclass
class PSFTopoRemark:
    topofilename: str = ''

    @classmethod
    def _from_remarkline(cls, remarkline: str):
        tokens = [x.strip() for x in remarkline.split()]
        assert tokens[0] == 'REMARKS'
        assert tokens[1] == 'topology'
        topofilename = tokens[2]
        return cls(topofilename=topofilename)
    
class PSFTopoRemarkList(UserList[PSFTopoRemark]):
    pass

@dataclass
class PSFPatchRemark:
    patchname: str = ''
    patchspecs: list[str] = field(default_factory=list)
    data: object = None

    @classmethod
    def _from_remarkline(cls, remarkline: str):
        tokens = [x.strip() for x in remarkline.split()]
        assert tokens[0] == 'REMARKS'
        assert tokens[1] == 'patch'
        patchname = tokens[2]
        patchspecs = tokens[3:]
        data = None
        if patchname == 'DISU':
            data = PSFDISUPatch(patchspecs)
        elif patchname in Link._patch_atomnames:
            data = PSFLinkPatch([patchname]+patchspecs)
        return cls(patchname=patchname, patchspecs=patchspecs, data=data)

class PSFPatchRemarkList(UserList[PSFPatchRemark]):
    pass

@dataclass
class PSFSegmentRemark:
    segname: str = ''
    segdata: str = ''

    @classmethod
    def _from_remarkline(cls, remarkline: str):
        tokens = [x.strip() for x in remarkline.split()]
        assert tokens[0] == 'REMARKS'
        assert tokens[1] == 'segment'
        segname = tokens[2]
        segdata = ' '.join(tokens[3:])
        return cls(segname=segname, segdata=segdata)

class PSFSegmentRemarkList(UserList[PSFSegmentRemark]):
    pass

@dataclass
class PSFRemark:
    remarkline: str = ''
    data: object = None

    @classmethod
    def from_remarkline(cls, remarkline: str):
        tokens = [x.strip() for x in remarkline.split()]
        assert tokens[0] == 'REMARKS'
        remarktype = tokens[1]
        # logger.debug(f'Parsing PSF remark of type {remarktype}')
        data = None
        if remarktype == 'patch':
            data = PSFPatchRemark._from_remarkline(remarkline)
            # logger.debug(f'Parsed patch remark: {type(data)}')
        elif remarktype == 'segment':
            data = PSFSegmentRemark._from_remarkline(remarkline)
        elif remarktype == 'topology':
            data = PSFTopoRemark._from_remarkline(remarkline)
        return cls(remarkline=remarkline, data=data)
    
class PSFRemarkList(UserList[PSFRemark]):
    def get_segmentremarks(self) -> PSFSegmentRemarkList:
        return PSFSegmentRemarkList([r.data for r in self if isinstance(r.data, PSFSegmentRemark)])
    def get_patchremarks(self) -> PSFPatchRemarkList:
        return PSFPatchRemarkList([r.data for r in self if isinstance(r.data, PSFPatchRemark)])
    def get_toporemarks(self) -> PSFTopoRemarkList:
        return PSFTopoRemarkList([r.data for r in self if isinstance(r.data, PSFTopoRemark)])