import argparse
import os
import json
import ROOT as r
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import uproot3
from matplotlib.offsetbox import AnchoredText
from util import make_dirs
import rhalphalib as rl


def prepareTF(
    rl_poly,
    xaxis=(40, 200),
    yaxis=(450, 1200),
    raxis=(-6, -2.1),
    ):
    fptbins, fmsdbins = np.linspace(yaxis[0], yaxis[-1], 500), np.linspace(xaxis[0], xaxis[-1], 500)
    frhobins = 2*np.log(fmsdbins/fptbins)

    fptpts, fmsdpts = np.meshgrid(fptbins, fmsdbins, indexing='ij')
    frhopts = 2*np.log(fmsdpts/fptpts)

    fptscaled = (fptpts - yaxis[0]) / (yaxis[-1] - yaxis[0])
    frhoscaled = (frhopts - (raxis[0])) / (raxis[1] - raxis[0])
    fvalidbins = (frhoscaled >= 0) & (frhoscaled <= 1)
    frhoscaled = np.clip(frhoscaled, 0, 1)

    TF = rl_poly(fptscaled, frhoscaled, nominal=True)
    TF = np.ma.array(TF, mask=~fvalidbins)
    return TF, fptpts, fmsdpts, frhopts

def base_plot(TF, fptpts, fmsdpts, frhopts, ax=None, rho=False):
    if ax is None:
        ax = plt.gca()
    coff = 0.3
    zmin, zmax = np.floor(10 * np.min(TF)) / 10, np.ceil(10 * np.max(TF)) / 10
    zmin, zmax = zmin + 0.001, zmax - 0.001
    clim = np.round(np.min([abs(zmin - 1), abs(zmax - 1)]), 1)
    clim = np.clip(clim, 0.2, 0.4)
    levels = np.linspace(1 - clim, 1 + clim, 500)
    
    if np.min(TF) < 1 - clim and np.max(TF) > 1 + clim:
        _extend = 'both'
    elif np.max(TF) > 1 + clim:
        _extend = 'max'
    elif np.min(TF) < 1 - clim:
        _extend = 'min'
    else:
        _extend = 'neither'
    if abs(1 - zmin) > coff and abs(1 - zmax) > coff:
        c_extend = 'both'
    elif abs(1 - zmin) > coff:
        c_extend = 'min'
    elif abs(1 - zmax) > coff:
        c_extend = 'max'
    else:
        c_extend = 'neither'
    
    pmesh = plt.contourf(frhopts if rho else fmsdpts, fptpts, TF, cmap='RdBu_r', levels=levels, extend=_extend, corner_mask=False)
    cax = hep.make_square_add_cbar(ax, pad=0.2, size=0.5)
    cbar = plt.colorbar(pmesh, cax=cax, extend=c_extend)
    cbar.set_ticks([np.arange(1 - clim, 1 + clim + 0.1, 0.1)])

    ax.set_xlim(-6, -2.1) if rho else ax.set_xlim(round(np.min(fmsdpts)), round(np.max(fmsdpts)))
    ax.set_ylim(round(np.min(fptpts)), round(np.max(fptpts)))
    ax.invert_yaxis()
    cax.set_ylabel("TF", y=1, ha='right')
    return ax


def plotTF_styling(ax, rho=False):
    def rho_bound(ms, rho):
        fpt = ms * np.e**(-rho / 2)
        return fpt

    x = np.arange(40, 70)
    ax.plot(x, rho_bound(x, -6), 'black', lw=3)
    ax.fill_between(x,
                    rho_bound(x, -6),
                    1200,
                    facecolor="none",
                    hatch="xx",
                    edgecolor="black",
                    linewidth=0.0)
    x = np.arange(150, 201)
    ax.plot(x, rho_bound(x, -2.1) + 5, 'black', lw=3)
    ax.fill_between(x,
                    rho_bound(x, -2.1),
                    facecolor="none",
                    hatch="xx",
                    edgecolor="black",
                    linewidth=0.0)

    _mbins = np.linspace(40, 201, 24)
    _pbins = np.array([450, 500, 550, 600, 675, 800, 1200])
    sampling_m, sampling_pt = np.meshgrid(_mbins[:-1] + 0.5 * np.diff(_mbins),
                                          _pbins[:-1] + 0.3 * np.diff(_pbins))
    sampling_rho = 2*np.log(sampling_m/sampling_pt)
    valmask = (sampling_pt > rho_bound(sampling_m, -2.1)) & (sampling_pt < rho_bound(
        sampling_m, -6))
    x_sampling = sampling_rho[valmask] if rho else sampling_m[valmask]
    ax.scatter(
        x_sampling,
        sampling_pt[valmask],
        marker='x',
        color='black',
        s=40,
        alpha=.4,
    )
    if rho:
        ax.set_xlabel(r'Jet $\rho$', ha='right', x=1)
    else:
        ax.set_xlabel(r'Jet $\mathrm{m_{SD}}$', ha='right', x=1)
    ax.set_ylabel(r'Jet $\mathrm{p_{T}}$', ha='right', y=1)
    return ax
    
def singleTF(tf_poly,  rho = False):
    TF, fptpts, fmsdpts, frhopts = prepareTF(tf_poly)    
    fig, ax = plt.subplots()
    base_plot(TF, fptpts, fmsdpts, frhopts, ax=ax, rho=rho)
    plotTF_styling(ax, rho=rho)

    _latexify = {'rho': '\\rho', 'pt': 'p_T'}
    _deg_str = "$({})$".format(",".join([_latexify[n] +"=" + str(o) for n, o in zip(tf_poly.dim_names, tf_poly.order)]))
    _label =  "{} {}".format(tf_poly.basis, _deg_str)

    at = AnchoredText(_label, loc='lower right', frameon=False, prop=dict(size=20))
    ax.add_artist(at);
    return ax

def combinedTF(tf1, tf2, rho=False):
    TF1, fptpts, fmsdpts, frhopts = prepareTF(tf1)
    TF2, _, _, _ = prepareTF(tf2)    
    fig, ax = plt.subplots()
    base_plot(TF1 * TF2, fptpts, fmsdpts, frhopts, ax=ax, rho=rho)
    plotTF_styling(ax, rho=rho)

    _latexify = {'rho': '\\rho', 'pt': 'p_T'}
    _deg_str = "$({})$".format(",".join([_latexify[n] +"=" + str(o) for n, o in zip(tf1.dim_names, tf1.order)]))
    _label =  "{} {}".format(tf1.basis, _deg_str)
    _deg_str = "$({})$".format(",".join([_latexify[n] +"=" + str(o) for n, o in zip(tf2.dim_names, tf2.order)]))
    _label += "\nx\n{} {}".format(tf2.basis, _deg_str)

    at = AnchoredText(_label, loc='lower right', frameon=False, prop=dict(size=20, multialignment='center'))
    ax.add_artist(at);
    return ax


def plotTF_ratio(in_ratio, mask, region, args=None, zrange=None):
    fig, ax = plt.subplots()

    H = np.ma.masked_where(in_ratio * mask <= 0.01, in_ratio * mask)
    zmin, zmax = np.nanmin(H), np.nanmax(H)
    if zrange is None:
        # Scale clim to fit range up to a max of 0.6
        clim = np.max([.3, np.min([0.6, 1 - zmin, zmax - 1])])
    else:
        clim = zrange
    ptbins = np.array([450, 500, 550, 600, 675, 800, 1200])
    msdbins = np.linspace(40, 201, 24)
    hep.hist2dplot(H.T,
                   msdbins,
                   ptbins,
                   vmin=1 - clim,
                   vmax=1 + clim,
                   cmap='RdBu_r',
                   cbar=False)
    cax = hep.make_square_add_cbar(ax, pad=0.2, size=0.5)
    if abs(1 - zmin) > .3 and abs(1 - zmax) > .3:
        c_extend = 'both'
    elif abs(1 - zmin) > .3:
        c_extend = 'min'
    elif abs(1 - zmax) > .3:
        c_extend = 'max'
    else:
        c_extend = 'neither'
    cbar = fig.colorbar(ax.get_children()[0], cax=cax, extend=c_extend)

    ax.set_xticks(np.arange(40, 220, 20))
    ax.tick_params(axis='y', which='minor', left=False, right=False)
    ax.invert_yaxis()

    ax.set_title('{} QCD Ratio'.format(region), pad=15, fontsize=26)
    ax.set_xlabel(r'Jet $\mathrm{m_{SD}}$', ha='right', x=1)
    ax.set_ylabel(r'Jet $\mathrm{p_{T}}$', ha='right', y=1)
    cbar.set_label(r'(Pass QCD) / (Fail QCD * eff)', ha='right', y=1)
    return ax
    

if __name__ == '__main__':
    plt.style.use([hep.cms.style.ROOT])
    plt.switch_backend('agg')
    np.seterr(divide='ignore', invalid='ignore')

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", default='', help="Model/Fit dir")
    parser.add_argument("-f",
                        "--fit",
                        default='fitDiagnostics.root',
                        dest='fit',
                        help="fitDiagnostics file")
    parser.add_argument("-o",
                        "--output-folder",
                        default='plots',
                        dest='output_folder',
                        help="Folder to store plots - will be created ? doesn't exist.")
    parser.add_argument("--year",
                        default=None,
                        choices={"2016", "2017", "2018"},
                        type=str,
                        help="year label")
    parser.add_argument("--MC",
                        action='store_true',
                        dest='isMC',
                        help="Use 'simulation' label")
    parser.add_argument("--rho",
                        action='store_true',
                        help="Plot rho as axis instead of msd")
    args = parser.parse_args()

    if args.output_folder.split("/")[0] != args.dir:
        args.output_folder = os.path.join(args.dir, args.output_folder)
    make_dirs(args.output_folder)

    configs = json.load(open(os.path.join(args.dir, "config.json")))
    if args.year is None:
        args.year = str(configs['year'])

    # Get fitDiagnostics File
    rf = r.TFile.Open(os.path.join(args.dir, args.fit))

    # Define TFs
    basis1, basis2 = configs['basis'].split(",")
    degs = tuple([int(s) for s in configs['degs'].split(',')])
    degsMC = tuple([int(s) for s in configs['degsMC'].split(',')])
    tf_MC = rl.BasisPoly("tf{}_MCtempl".format(args.year), degsMC, ['pt', 'rho'], basis=basis1)
    tf_res = rl.BasisPoly("tf{}_dataResidual".format(args.year), degs, ['pt', 'rho'], basis=basis2)

    # Set to values from fit
    tf_res.update_from_roofit(rf.Get('fit_s'))
    par_names = rf.Get('fit_s').floatParsFinal().contentsString().split(',')
    MC_nuis = [round(rf.Get('fit_s').floatParsFinal().find(p).getVal(), 3) for p in par_names if 'tf{}_MCtempl'.format(args.year) in p]
    _vect = np.load(os.path.join(args.dir, 'decoVector.npy'))
    _MCTF_nominal = np.load(os.path.join(args.dir, 'MCTF.npy'))
    _values = _vect.dot(np.array(MC_nuis)) + _MCTF_nominal
    tf_MC.set_parvalues(_values)

    ax = singleTF(tf_MC, rho=args.rho)
    ax.set_title("Tagger Response TF", x=0, ha='left', fontsize='small')
    hep.cms.label(ax=ax, loc=2, year=args.year)
    ax.figure.savefig('{}/TF_MC.png'.format(args.output_folder), dpi=300, bbox_inches="tight")
    ax.figure.savefig('{}/TF_MC.pdf'.format(args.output_folder), transparent=True, bbox_inches="tight")

    ax = singleTF(tf_res, rho=args.rho)
    ax.set_title("Residual (Data/MC) TF", x=0, ha='left', fontsize='small')
    hep.cms.label(ax=ax, loc=2, year=args.year)
    ax.figure.savefig('{}/TF_data.png'.format(args.output_folder), dpi=300, bbox_inches="tight")
    ax.figure.savefig('{}/TF_data.pdf'.format(args.output_folder), transparent=True, bbox_inches="tight")

    ax = combinedTF(tf_MC, tf_res, rho=args.rho)
    ax.set_title("Effective Transfer Factor", x=0, ha='left', fontsize='small')
    hep.cms.label(ax=ax, loc=2, year=args.year)
    ax.figure.savefig('{}/TF_eff.png'.format(args.output_folder), dpi=300, bbox_inches="tight")
    ax.figure.savefig('{}/TF_eff.pdf'.format(args.output_folder), transparent=True, bbox_inches="tight")


    f = uproot3.open(os.path.join(args.dir, args.fit))
    region = 'prefit'
    fail_qcd, pass_qcd = [], []
    bins = []
    for ipt in range(6):
        fail_qcd.append(f['shapes_{}/ptbin{}{}{}/qcd;1'.format(
            'prefit', ipt, 'fail', args.year)].values)
        pass_qcd.append(f['shapes_{}/ptbin{}{}{}/qcd;1'.format(
            'prefit', ipt, 'pass', args.year)].values)

    fail_qcd = np.array(fail_qcd)
    pass_qcd = np.array(pass_qcd)
    mask = ~np.isclose(pass_qcd, np.zeros_like(pass_qcd))
    mask *= ~np.isclose(fail_qcd, np.zeros_like(fail_qcd))
    q = np.sum(pass_qcd[mask]) / np.sum(fail_qcd[mask])
    in_data_rat = (pass_qcd / (fail_qcd * q))

    ax = plotTF_ratio(in_data_rat, mask, region="Prefit", args=args)
    ax.figure.savefig('{}/{}{}.png'.format(args.output_folder, "TF_ratio_", region),
                      bbox_inches="tight", dpi=300)
    ax.figure.savefig('{}/{}{}.pdf'.format(args.output_folder, "TF_ratio_", region),
                      bbox_inches="tight", transparent=True)

    region = 'postfit'
    fail_qcd, pass_qcd = [], []
    bins = []
    for ipt in range(6):
        fail_qcd.append(f['shapes_{}/ptbin{}{}{}/qcd;1'.format(
            'fit_s', ipt, 'fail', args.year)].values)
        pass_qcd.append(f['shapes_{}/ptbin{}{}{}/qcd;1'.format(
            'fit_s', ipt, 'pass', args.year)].values)

    fail_qcd = np.array(fail_qcd)
    pass_qcd = np.array(pass_qcd)

    mask = ~np.isclose(pass_qcd, np.zeros_like(pass_qcd))
    mask *= ~np.isclose(fail_qcd, np.zeros_like(fail_qcd))
    q = np.sum(pass_qcd[mask]) / np.sum(fail_qcd[mask])
    in_data_rat = (pass_qcd / (fail_qcd * q))

    ax = plotTF_ratio(in_data_rat, mask, region="Postfit", args=args)
    ax.figure.savefig('{}/{}{}.png'.format(args.output_folder, "TF_ratio_", region),
                      bbox_inches="tight", dpi=300)
    ax.figure.savefig('{}/{}{}.pdf'.format(args.output_folder, "TF_ratio_", region),
                      bbox_inches="tight", transparent=True)
