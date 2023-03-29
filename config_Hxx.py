import numpy as np

# Define Bins
def get_bins(year="2017"):
    if year == "2017":
        ptbins = np.array([475, 500, 550, 600, 675, 800, 1200])
    else:
        ptbins = np.array([450, 500, 550, 600, 675, 800, 1200])
    msdbins = np.linspace(40, 201, 24)

    npt = len(ptbins) - 1

    # Define pt/msd/rho grids
    ptpts, msdpts = np.meshgrid(ptbins[:-1] + 0.3 * np.diff(ptbins),
                                msdbins[:-1] + 0.5 * np.diff(msdbins),
                                indexing='ij')
    rhopts = 2*np.log(msdpts/ptpts)
    ptscaled = (ptpts - ptbins[0]) / (ptbins[-1] - ptbins[0])
    rhoscaled = (rhopts - (-6)) / ((-2.1) - (-6))
    validbins = (rhoscaled >= 0) & (rhoscaled <= 1)
    rhoscaled[~validbins] = 1  # we will mask these out later

    # Define fine bins for smooth TF plots
    fptbins = np.arange(ptbins[0], ptbins[0]+2, 2)
    fmsdbins = np.arange(40, 201.5, .5)

    fptpts, fmsdpts = np.meshgrid(fptbins[:-1] + 0.3 * np.diff(fptbins),
                                fmsdbins[:-1] + 0.5 * np.diff(fmsdbins),
                                indexing='ij')
    frhopts = 2*np.log(fmsdpts/fptpts)
    fptscaled = (fptpts - ptbins[0]) / (ptbins[-1] - ptbins[0])
    frhoscaled = (frhopts - (-6)) / ((-2.1) - (-6))
    fvalidbins = (frhoscaled >= 0) & (frhoscaled <= 1)
    frhoscaled[~fvalidbins] = 1  # we will mask these out later

    return msdbins, ptbins, ptscaled, rhoscaled, validbins

msdbins, ptbins, ptscaled, rhoscaled, validbins = get_bins(year="2016")


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
        # 2016: {
        #     # 'cc': {
        #     #     "SF": [0.77],
        #     #     "up": [0.15],
        #     #     "down": [0.15],
        #     # }
        #     'cc': {
        #         "SF": [1.09],
        #         "up": [0.1],
        #         "down": [0.1],
        #     }
        # },
        # 2017: {
        #     'cc': {
        #         "SF": [1.24],
        #         "up": [0.2],
        #         "down": [0.2],
        #     }
        # },
        # 2018: {
        #     'cc': {
        #         "SF": [0.94],
        #         "up": [0.15],
        #         "down": [0.15],
        #     }
        # },
        2016: {
            'cc': {
                "SF": [1.15],
                "up": [0.25],
                "down": [0.25],
            }
        },
        2017: {
            'cc': {
                "SF": [0.85],
                "up": [0.16],
                "down": [0.16],
            }
        },
        2018: {
            'cc': {
                "SF": [0.74],
                "up": [0.2],
                "down": [0.2],
            }
        },
    }
    if flav not in ['bb', 'cc', 'qq']:
        raise ValueError("``flav`` has be one of ['bb', 'cc', 'qq'].")
    if pbin in [0, 1, 2, 3, 4, 5]:
        return _SFdict[year][flav]['SF'][0], _SFdict[year][flav]['up'][0], _SFdict[year][flav]['down'][0]
    else:
        raise RuntimeError()


def ddxSFUL(pbin, flav, year=2017):
    # ptbins = np.array([450, 500, 550, 600, 675, 800, 1200])
    _SFdict = {
        2016: {
            'cc': {
                "SF": [1., 0.97, 0.624],
                "up": [0.15, 0.12, 0.4],
                "down": [0.17, 0.16, 0.25],
            }
        },
        2017: {
            'cc': {
                "SF": [1.18, 1.04, 1.18],
                "up": [0.17, 0.13, 0.18],
                "down": [0.17, 0.13, 0.18],
            }
        },
        2018: {
            'cc': {
                "SF": [0.99, 0.98, 0.92],
                "up": [0.12, 0.14, 0.23],
                "down": [0.12, 0.14, 0.23],
            }
        },
    }
    if flav not in ['bb', 'cc', 'qq']:
        raise ValueError("``flav`` has be one of ['bb', 'cc', 'qq'].")
    if pbin in [0, 1, 2]:
        return _SFdict[year][flav]['SF'][0], _SFdict[year][flav]['up'][0], _SFdict[year][flav]['down'][0]
    elif pbin in [3, 4]:
        return _SFdict[year][flav]['SF'][1], _SFdict[year][flav]['up'][1], _SFdict[year][flav]['down'][1]
    elif pbin in [5]:
        return _SFdict[year][flav]['SF'][2], _SFdict[year][flav]['up'][2], _SFdict[year][flav]['down'][2]
    else:
        raise RuntimeError()