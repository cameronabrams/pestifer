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

__colormapnames__ = list(mpl.colormaps.keys())

class PestiferColorMap(UserList):
    """
    A class for handling color definitions in Pestifer, based on a specified matplotlib colormap.
    """
    def __init__(self, mpl_colormapname):
        if mpl_colormapname not in __colormapnames__:
            return None
        cmap = mpl.colormaps[mpl_colormapname]
        self.data = []
        for c in cmap.colors:
            cc = np.array(np.array(c) * 255, dtype=int)
            self.data.append(ColorRGB(*cc))

    def hue(self, x: float):
        """
        Get a color from the colormap based on a float input in the range [0, 1].
        
        Parameters
        ----------
        x : float
            A float value in the range [0, 1] to select a color from the colormap.
        
        Returns
        -------
        ColorRGB
            The corresponding ColorRGB object from the colormap.
        """
        if x < 0.0 or x > 1.0:
            raise ValueError("Input must be in the range [0, 1]")
        index = int(x * (len(self.data) - 1))
        return self.data[index]


__plasma__ = PestiferColorMap('plasma')

PestiferColors = dict(
    adobe = ColorHex("#c04737"),
    alabaster = ColorHex("#f2f0e6"),
    alizarin = ColorHex("#e74c3c"),
    aluminum = ColorHex("#b2b2b2"),
    amber = ColorHex("#ffbf00"),
    amethyst = ColorHex("#9966cc"),
    apricot = ColorHex("#fbceb1"),
    aqua = ColorHex("#00ffff"),
    arctic = ColorHex("#b0e0e6"),
    ash = ColorHex("#b2beb5"),
    azure = ColorHex("#f0ffff"),
    bamboo = ColorHex("#c7d295"),
    basalt = ColorHex("#4b5563"),
    beige = ColorHex("#f5f5dc"),
    bisque = ColorHex("#ffe4c4"),
    black = ColorHex("#000000"),
    blue = ColorHex("#0000ff"),
    bronze = ColorHex("#cd7f32"),
    brown = ColorHex("#a52a2a"),
    burgundy = ColorHex("#800020"),
    cactus = ColorHex("#7d7c84"),
    caramel = ColorHex("#af6f09"),
    carbon = ColorHex("#000000"),
    celadon = ColorHex("#ace1af"),
    cerulean = ColorHex("#007ba7"),
    champagne = ColorHex("#f7e7ce"),
    charcoal = ColorHex("#36454f"),
    chartreuse = ColorHex("#7fff00"),
    chromium = ColorHex("#8a8a8a"),
    cinnabar = ColorHex("#e34234"),
    cobalt = ColorHex("#1338be"),
    copper = ColorHex("#b87333"),
    coral = ColorHex("#ff7f50"),
    cosmos = ColorHex("#ffd1dc"),
    cream = ColorHex("#fffdd0"),
    crimson = ColorHex("#dc143c"),
    cyan = ColorHex("#00ffff"),
    dandelion = ColorHex("#fed85d"),
    darkgray = ColorHex("#a9a9a9"),
    denim = ColorHex("#1560bd"),
    dove = ColorHex("#6b6b47"),
    driftwood = ColorHex("#af8751"),
    ebony = ColorHex("#555d50"),
    eggplant = ColorHex("#614051"),
    elderberry = ColorHex("#6e3f7e"),
    emerald = ColorHex("#50c878"),
    eucalyptus = ColorHex("#278a5b"),
    fig = ColorHex("#715c87"),
    flamingo = ColorHex("#fc8eac"),
    fluorescent_blue = ColorHex("#1f51ff"),
    fluorescent_green = ColorHex("#39ff14"),
    fluorescent_orange = ColorHex("#ff5f1f"),
    fluorescent_pink = ColorHex("#ff6ec7"),
    fluorescent_purple = ColorHex("#bf00ff"),
    fluorescent_red = ColorHex("#ff073a"),
    fluorescent_yellow = ColorHex("#ccff00"),
    forest = ColorHex("#228b22"),
    frost = ColorHex("#f7f7f7"),
    fuchsia = ColorHex("#ff00ff"),
    garnet = ColorHex("#733635"),
    ginger = ColorHex("#b06500"),
    glacier = ColorHex("#9bb0ff"),
    gold = ColorHex("#ffbf00"),
    granite = ColorHex("#676767"),
    gray = ColorHex("#808080"),
    green = ColorHex("#00ff00"),
    harlequin = ColorHex("#3fff00"),
    haze = ColorHex("#8b7d6b"),
    hibiscus = ColorHex("#c7375f"),
    honey = ColorHex("#ffc30b"),
    iceberg = ColorHex("#b9d9eb"),
    indigo = ColorHex("#4b0082"),
    ink = ColorHex("#1c1c1c"),
    iris = ColorHex("#5a4fcf"),
    iron = ColorHex("#b7410e"),
    ivory = ColorHex("#fffff0"),
    jade = ColorHex("#00a86b"),
    jasmine = ColorHex("#f8de7e"),
    jasper = ColorHex("#d73502"),
    juniper = ColorHex("#6b9b73"),
    kale = ColorHex("#4d7c0f"),
    kelp = ColorHex("#454b1b"),
    khaki = ColorHex("#f0e68c"),
    kiwi = ColorHex("#8ee53f"),
    lapis = ColorHex("#26619c"),
    lavender = ColorHex("#e39ff6"),
    lead = ColorHex("#575961"),
    lemon = ColorHex("#fff700"),
    lightgray = ColorHex("#d3d3d3"),
    lightning = ColorHex("#e6e6fa"),
    lilac = ColorHex("#c8a2c8"),
    lime = ColorHex("#bfff00"),
    magenta = ColorHex("#ff00ff"),
    mahogany = ColorHex("#c04000"),
    mango = ColorHex("#fdbe02"),
    manganese = ColorHex("#9c7c5d"),
    maroon = ColorHex("#800000"),
    mauve = ColorHex("#e0b0ff"),
    mercury = ColorHex("#d3d3d3"),
    mint = ColorHex("#98ff98"),
    moonstone = ColorHex("#3aa8c1"),
    moss = ColorHex("#addfad"),
    navy = ColorHex("#000080"),
    nebula = ColorHex("#663399"),
    nectarine = ColorHex("#ffb347"),
    neon_blue = ColorHex("#1b03a3"),
    nickel = ColorHex("#727472"),
    nutmeg = ColorHex("#81613c"),
    obsidian = ColorHex("#0b1426"),
    ochre = ColorHex("#cc7722"),
    olive = ColorHex("#808000"),
    onyx = ColorHex("#353839"),
    opal = ColorHex("#a8c3bc"),
    orange = ColorHex("#ff8000"),
    ozone = ColorHex("#6699cc"),
    papaya = ColorHex("#ffefd5"),
    peach = ColorHex("#ffe5b4"),
    periwinkle = ColorHex("#c5c5ff"),
    phosphorus = ColorHex("#ff8000"),
    pine = ColorHex("#01796f"),
    pink = ColorHex("#ffc0cb"),
    platinum = ColorHex("#e5e4e2"),
    prism = ColorHex("#ff6ec7"),
    purple = ColorHex("#800080"),
    quail = ColorHex("#b6aa9c"),
    quasar = ColorHex("#6a5acd"),
    quartz = ColorHex("#51484f"),
    quince = ColorHex("#ee8a00"),
    radiance = ColorHex("#fff68f"),
    red = ColorHex("#ff0000"),
    rhubarb = ColorHex("#c21807"),
    rose = ColorHex("#e4234b"),
    ruby = ColorHex("#e0115f"),
    rust = ColorHex("#b7410e"),
    saffron = ColorHex("#f4c430"),
    sage = ColorHex("#9caf88"),
    salmon = ColorHex("#fa8072"),
    sapphire = ColorHex("#0f52ba"),
    seafoam = ColorHex("#3ded97"),
    silicon = ColorHex("#c2b280"),
    silver = ColorHex("#c0c0c0"),
    stardust = ColorHex("#9bb3c0"),
    steel = ColorHex("#4682b4"),
    sulfur = ColorHex("#ffff30"),
    tan = ColorHex("#d2b48c"),
    tangerine = ColorHex("#f28500"),
    taupe = ColorHex("#483c32"),
    teal = ColorHex("#008080"),
    thistle = ColorHex("#d8bfd8"),
    thunder = ColorHex("#708090"),
    tin = ColorHex("#9ea0a1"),
    titanium = ColorHex("#8f9779"),
    topaz = ColorHex("#ffc87c"),
    turquoise = ColorHex("#40e0d0"),
    ultramarine = ColorHex("#120a8f"),
    umber = ColorHex("#635147"),
    universe = ColorHex("#2e2d88"),
    urchin = ColorHex("#404040"),
    vanilla = ColorHex("#f3e5ab"),
    verbena = ColorHex("#da70d6"),
    vermillion = ColorHex("#e34234"),
    violet = ColorHex("#ee82ee"),
    vortex = ColorHex("#da70d6"),
    walnut = ColorHex("#773f1a"),
    watermelon_red = ColorHex("#ff5733"),
    wheat = ColorHex("#f5deb3"),
    whirlpool = ColorHex("#00868b"),
    white = ColorHex("#ffffff"),
    willow = ColorHex("#9cbb58"),
    wisteria = ColorHex("#c9a0dc"),
    xanthe = ColorHex("#f1e788"),
    xerus = ColorHex("#cd853f"),
    xylem = ColorHex("#8fbc8f"),
    yam = ColorHex("#d19fe8"),
    yarrow = ColorHex("#ffe784"),
    yellow = ColorHex("#ffff00"),
    yucca = ColorHex("#faf0e6"),
    zaffre = ColorHex("#0014a8"),
    zephyr = ColorHex("#4f97a3"),
    zinc = ColorHex("#7f7f7f"),
    zucchini = ColorHex("#4d5d53"))
