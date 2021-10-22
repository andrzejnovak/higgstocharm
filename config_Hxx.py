import numpy as np

# Define Bins
ptbins = np.array([450, 500, 550, 600, 675, 800, 1200])
npt = len(ptbins) - 1
msdbins = np.linspace(40, 201, 24)

# Define pt/msd/rho grids
ptpts, msdpts = np.meshgrid(ptbins[:-1] + 0.3 * np.diff(ptbins),
                            msdbins[:-1] + 0.5 * np.diff(msdbins),
                            indexing='ij')
rhopts = 2*np.log(msdpts/ptpts)
ptscaled = (ptpts - 450.) / (1200. - 450.)
rhoscaled = (rhopts - (-6)) / ((-2.1) - (-6))
validbins = (rhoscaled >= 0) & (rhoscaled <= 1)
rhoscaled[~validbins] = 1  # we will mask these out later

# Define fine bins for smooth TF plots
fptbins = np.arange(450, 1202, 2)
fmsdbins = np.arange(40, 201.5, .5)

fptpts, fmsdpts = np.meshgrid(fptbins[:-1] + 0.3 * np.diff(fptbins),
                              fmsdbins[:-1] + 0.5 * np.diff(fmsdbins),
                              indexing='ij')
frhopts = 2*np.log(fmsdpts/fptpts)
fptscaled = (fptpts - 450.) / (1200. - 450.)
frhoscaled = (frhopts - (-6)) / ((-2.1) - (-6))
fvalidbins = (frhoscaled >= 0) & (frhoscaled <= 1)
frhoscaled[~fvalidbins] = 1  # we will mask these out later


def TF_params(xparlist, xparnames=None, nrho=None, npt=None):
    # TF map from param/name lists
    if xparnames is not None:
        from operator import methodcaller

        def _get(s):
            return s[-1][0]

        ptdeg = max(
            list(
                map(
                    int,
                    list(
                        map(_get, list(map(methodcaller("split", 'pt_par'),
                                           xparnames)))))))
        rhodeg = max(
            list(
                map(
                    int,
                    list(
                        map(_get, list(map(methodcaller("split", 'rho_par'),
                                           xparnames)))))))
    else:
        rhodeg, ptdeg = nrho, npt

    TF_cf_map = np.array(xparlist).reshape(rhodeg + 1, ptdeg + 1)

    return TF_cf_map, rhodeg, ptdeg


# SFs
# def ddxSF(pbin, flav):
#     # ptbins = np.array([450, 500, 550, 600, 675, 800, 1200])
#     _SFdict = {
#         'cc': {
#             "SF": [0.899, 1.152, 0.692],
#             "up": [0.254, 0.428, 0.309],
#             "down": [0.254, 0.426, 0.278],
#         }
#     }
#     if flav not in ['bb', 'cc', 'qq']:
#         raise ValueError("``flav`` has be one of ['bb', 'cc', 'qq'].")
#     if pbin in [0, 1, 2]:
#         return _SFdict[flav]['SF'][0], _SFdict[flav]['up'][0], _SFdict[flav]['down'][0]
#     elif pbin in [3, 4]:
#         return _SFdict[flav]['SF'][1], _SFdict[flav]['up'][1], _SFdict[flav]['down'][1]
#     elif pbin in [5]:
#         return _SFdict[flav]['SF'][2], _SFdict[flav]['up'][2], _SFdict[flav]['down'][2]
#     else:
#         raise RuntimeError()

def ddxSF(pbin, flav, year=2017):
    # ptbins = np.array([450, 500, 550, 600, 675, 800, 1200])
    _SFdict = {
        2016: {
            'cc': {
                "SF": [0.77],
                "up": [0.15],
                "down": [0.15],
            }
        },
        2017: {
            'cc': {
                "SF": [1.24],
                "up": [0.2],
                "down": [0.2],
            }
        },
        2018: {
            'cc': {
                "SF": [0.94],
                "up": [0.15],
                "down": [0.15],
            }
        },
    }
    if flav not in ['bb', 'cc', 'qq']:
        raise ValueError("``flav`` has be one of ['bb', 'cc', 'qq'].")
    if pbin in [0, 1, 2, 3, 4, 5]:
        return _SFdict[year][flav]['SF'][0], _SFdict[year][flav]['up'][0], _SFdict[year][flav]['down'][0]
    else:
        raise RuntimeError()
