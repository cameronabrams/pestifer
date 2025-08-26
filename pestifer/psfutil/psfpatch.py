from ..objs.resid import ResID

class PSFDISUPatch:
    def __init__(self, pl: list[str]):
        """
        Initialize a PSFDISUPatch object from a list of strings.
        
        Parameters
        ----------
        pl : list[str]
            A list of strings representing the disulfide bond patch, typically in the format:
            ['segid1:resid1', 'segid2:resid2'].
        """
        self.seg1, self.res1 = pl[0].split(':')
        self.seg2, self.res2 = pl[1].split(':')
        self.resid1 = ResID(self.res1)
        self.resid2 = ResID(self.res2)

class PSFLinkPatch:
    def __init__(self, pl: list[str]):
        """
        Initialize a PSF Link Patch object from a list of strings.
        
        Parameters
        ----------
        pl : list[str]
            A list of strings representing the link patch, typically in the format:
            ['linkname', 'segid1:resid1', 'segid2:resid2'].
        """
        self.patchname = pl[0]
        self.seg1, self.res1 = pl[1].split(':')
        self.seg2, self.res2 = pl[2].split(':')
        self.resid1 = ResID(self.res1)
        self.resid2 = ResID(self.res2)