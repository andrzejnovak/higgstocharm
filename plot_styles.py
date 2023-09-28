from collections import OrderedDict
import copy

def scale_lightness(rgb, scale_l):
    import colorsys
    import matplotlib
    if isinstance(rgb, str):
        rgb = matplotlib.colors.to_rgba(rgb)[:-1]
    # convert rgb to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as rgb
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s=s)


def colorize(snames, mmap=None, cmap=None):
    import matplotlib
    if mmap is None:
        from itertools import groupby
        grouped = [list(g) for k, g in groupby(snames, key=lambda x: x.split("_")[0])]
    else:
        grouped = [[g for g in group if g in snames] for group in mmap.values()]
        grouped = sorted([group for group in grouped if len(group) > 0])
    if cmap is None:
        cmap = matplotlib.cm.get_cmap('tab10')
    elif isinstance(cmap, list):
        cmap = matplotlib.colors.ListedColormap(cmap)
    cdict = {}
    for i, group in enumerate(grouped):
        for samp, scale in zip(sorted(group),
                               [1, 1.2, 1.4, 1.6, 1.8, 0.9, 1.1, 1.3, 1.5]):
            cdict[samp] = scale_lightness(cmap(i)[:-1], scale)
    return cdict


merge_dict = {
    'top': ['stqq', 'tqq'],
    'vlep': ['zll', 'wln'],
    'hcc': ['hcc', 'vbfhcc', 'whcc', 'zhcc'],
    'hbb': ['hbb', 'vbfhbb', 'whbb', 'zhbb', 'tthbb'],
    'bkg': ['stqq', 'tqq', 'zll', 'wln', "wcq", "wqq", "zbb", "zqq"],
    'bkgNoW': ['stqq', 'tqq', 'zll', 'wln', "wqq", "zbb", "zqq"],
    'nonV': ['stqq', 'tqq', 'zll', 'wln'],
}

style_set_A = {
    'color_dict': {
        'hqq': '#003049',
        'hbb': '#4E5A7E',
        'hcc': '#3a86ff',
        'wqq': "#84a98c",
        'wcq': "#588157",
        'qcd': 'gray',
        'tqq': 'plum',
        'stqq': 'lightblue',
        'top': 'gray',
        'vlep': '#DEC1FF',
        'zbb': '#f77f00',
        'zcc': '#d62828',
        'zqq': '#fcbf49',
        'vvqq': '#C5C9A4',
        'nonV': 'gray',
        'bkg': 'gray',
        'bkgNoW': 'gray',
    },
    'hatch_dict': {
        'hqq': None,
        'hbb': None,
        'hcc': None,
        'wqq': None,
        'wcq': None,
        'qcd': None,
        'tqq': '\\\\\\\\',
        'stqq': '\\\\\\\\',
        'top': None,
        'vlep': '///',
        'zbb': None,
        'zcc': None,
        'zqq': None,
        'vvqq': None,
        'bkg': '\\\\\\\\',
        'nonV': '\\\\\\\\',
        'bkgNoW': '\\\\\\\\',
    },
    # Sequence of tuples because python2 is stupid
    'label_dict':
    OrderedDict([
        ('Data', 'Data'),
        ('MC', 'MC'),
        ('Toys', 'PostFit\nToys'),
        ('zcc', "$\mathrm{Z(c\\bar{c})}$"),
        ('zbb', "$\mathrm{Z(b\\bar{b})}$"),
        ('zqq', "$\mathrm{Z(q\\bar{q})}$"),
        ('wcq', "$\mathrm{W(c\\bar{q})}$"),
        ('wqq', "$\mathrm{W(q\\bar{q})}$"),
        ('hbb', "$\mathrm{H(b\\bar{b})}$"),
        ('hcc', "$\mathrm{H(c\\bar{c})}$"),
        ('vvqq', "$\mathrm{VV}$"),
        ('top', "$\mathrm{Top}$"),
        ('tqq', "$\mathrm{t\\bar{t}}$"),
        ('zll', "$\mathrm{DY}$"),
        ('vlep', "$\mathrm{V(lep.)}$"),
        ('stqq', "$\mathrm{single-t}$"),
        ('bkg', "Other"),
        ('nonV', "Other"),
        ('bkgNoW', "Other"),
        ('qcd', "QCD"),
    ])
}

style_set_C = copy.deepcopy(style_set_A)
style_set_C['hatch_dict']['wcq'] = "////"

style_set_B = {
    'color_dict': {
        'hqq': '#6479B9',
        'hbb': '#6479B9',
        'hcc': '#EE2F36',
        'wqq': '#6CAE75',
        'wcq': '#007A7A',
        'qcd': 'gray',
        'tqq': 'plum',
        'stqq': 'lightblue',
        'top': 'gray',
        'vlep': '#DEC1FF',
        'zbb': '#2C497F',
        'zcc': '#A4243B',
        'zqq': '#E09F3E',
        'vvqq': '#C5C9A4',
        'bkg': 'gray',
        'nonV': 'gray',
        'bkgNoW': 'gray',
    },
    'hatch_dict': {
        'hqq': None,
        'hbb': None,
        'hcc': None,
        'wqq': '///',
        'wcq': '///',
        'qcd': None,
        'tqq': '///',
        'stqq': '///',
        'top': '///',
        'vlep': '///',
        'zbb': None,
        'zcc': None,
        'zqq': None,
        'vvqq': '///',
        'bkg': '///',
        'nonV': '///',
        'bkgNoW': '///',
    },
    # Sequence of tuples because python2 is stupid
    'label_dict':
    OrderedDict([
        ('Data', 'Data'),
        ('MC', 'MC'),
        ('Toys', 'PostFit\nToys'),
        ('hbb', "$\mathrm{H(b\\bar{b})}$"),
        ('hqq', "$\mathrm{H(b\\bar{b})}$"),
        ('zbb', "$\mathrm{Z(b\\bar{b})}$"),
        ('zcc', "$\mathrm{Z(c\\bar{c})}$"),
        ('zqq', "$\mathrm{Z(q\\bar{q})}$"),
        ('wcq', "$\mathrm{W(c\\bar{q})}$"),
        ('wqq', "$\mathrm{W(q\\bar{q})}$"),
        ('vvqq', "$\mathrm{VV}$"),
        ('top', "$\mathrm{Top}$"),
        ('tqq', "$\mathrm{t\\bar{t}}$"),
        ('zll', "$\mathrm{DY}$"),
        ('vlep', "$\mathrm{V(lep.)}$"),
        ('stqq', "$\mathrm{single-t}$"),
        ('hcc', "$\mathrm{H(c\\bar{c})}$"),
        ('bkg', "Other"),
        ('nonV', "Other"),
        ('bkgNoW', "Other"),
        ('qcd', "QCD"),
    ])
}