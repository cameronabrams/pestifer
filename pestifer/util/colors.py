# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Definition of the :class:`PestiferColor` class for handling color definitions in Pestifer.
This class is used to create a list of colors based on a specified matplotlib colormap.
It inherits from :class:`collections.UserList` and provides a way to access colors in a structured manner.
"""
from collections import UserList
from colorist import ColorHex, ColorRGB
import matplotlib as mpl
import numpy as np

__colormapnames__=list(mpl.colormaps.keys())

class PestiferColor(UserList):
    """
    A class for handling color definitions in Pestifer, based on a specified matplotlib colormap.
    """
    def __init__(self,mpl_colormapname):
        if mpl_colormapname not in __colormapnames__:
            return None
        cmap=mpl.colormaps[mpl_colormapname]
        self.data=[]
        for c in cmap.colors:
            cc=np.array(np.array(c)*255,dtype=int)
            self.data.append(ColorRGB(*cc))

__plasma__=PestiferColor('plasma')

lavender=ColorHex("#e39ff6")
cobalt=ColorHex("#1338be")
rose=ColorHex("#e4234b")
watermelon_red=ColorHex("#ff5733")
seafoam=ColorHex("#3ded97")
