import os
import sys
import time
import argparse
import json

import ROOT as r
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import mplhep as hep

plt.style.use(hep.style.ROOT)

BASE = os.getcwd()

def get_vals(fname):
    rfile = r.TFile.Open(fname)
    rtree = rfile.Get("limit")
    vals = []
    for i in range(rtree.GetEntries()):
        rtree.GetEntry(i)
        mu = rtree.limit
        vals.append(mu)
    return vals

def fval(lambda1, lambda2, p1, p2, nbins):

    numerator = -2.0 * np.log(np.array(lambda1) / np.array(lambda2)) / (p2 - p1)
    denominator = -2.0 * np.log(np.array(lambda2)) / (nbins - p2)

    return numerator / denominator

def fplot(dname, ref, alt, year=2017, savename='fplotX', nbins=130):
    ref_pt, ref_rho = ref
    alt_pt, alt_rho = alt
    p1 = (ref_pt + 1) * (ref_rho + 1)
    p2 = (alt_pt + 1) * (alt_rho + 1)
    
    alt = get_vals('{dname}/bkgtest_{ref_pt}-{ref_rho}_{alt_pt}-{alt_rho}/gofsalt.root'
                   .format(dname=dname, ref_pt=ref_pt, ref_rho=ref_rho, alt_pt=alt_pt, alt_rho=alt_rho))
    base = get_vals('{dname}/bkgtest_{ref_pt}-{ref_rho}_{alt_pt}-{alt_rho}/gofsbase.root'
                    .format(dname=dname, ref_pt=ref_pt, ref_rho=ref_rho, alt_pt=alt_pt, alt_rho=alt_rho))
    if len(alt) != len(base):
        raise ValueError("Number of toys for base and ref does not match.")
    fvals = fval(base, alt, 2, 3, nbins)
    f_data = fval(get_vals('{dname}/bkgtest_{ref_pt}-{ref_rho}_{alt_pt}-{alt_rho}/refbase.root'
                           .format(dname=dname, ref_pt=ref_pt, ref_rho=ref_rho, alt_pt=alt_pt, alt_rho=alt_rho)),
                  get_vals('{dname}/bkgtest_{ref_pt}-{ref_rho}_{alt_pt}-{alt_rho}/refalt.root'
                           .format(dname=dname, ref_pt=ref_pt, ref_rho=ref_rho, alt_pt=alt_pt, alt_rho=alt_rho)),
                  2, 3, nbins)[0]

    from scipy.stats import f
    x_lim = max(np.percentile(fvals, 90), f_data*1.2)
    x = np.linspace(0, x_lim, 200)
    bins = np.linspace(0, x_lim, 30)
    width = bins[1] - bins[0]

    fig, ax = plt.subplots()
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.plot(x, len(base) * width * f.pdf(x, p2 - p1, nbins - p2), color='red', label='F-dist, ndf({},{})'.format(p2-p1, nbins-p2))
    ax.hist(fvals, bins, facecolor='none', edgecolor='black', histtype='stepfilled', lw=2,
            label="Toys, N = {}".format(len(fvals)))
    ax.hist(fvals[fvals > f_data], bins, facecolor='steelblue', edgecolor='gray', histtype='stepfilled', alpha=0.3,
            label='p-value = {}'.format(round(float(len(fvals[fvals > f_data]))/len(fvals), 3)));
    ax.annotate("", xy=(f_data, 0), xycoords=trans,
                xytext=(f_data, 0.25), textcoords=trans,
                arrowprops=dict(lw='4', color='b', arrowstyle="->,head_length=1.5,head_width=0.5"),
                )
    ax.plot([], [], color='blue', lw=2, label="Observed = {:.3f}".format(f_data))

    title = "TF({},{}) x TF({},{})".format(ref_pt, ref_rho, alt_pt, alt_rho)
    ax.legend(title=title)
    hep.cms.label(data=True, year=year, ax=ax)
    ax.set_xlim(0, x_lim)

    fig.savefig('{}.pdf'.format(savename), dpi=300, transparent=True, bbox_inches='tight')
    fig.savefig('{}.png'.format(savename), dpi=300, transparent=True, bbox_inches='tight')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make plots for ftest/gof files.')
    parser.add_argument('dname', nargs='?', default=os.getcwd())
    parser.add_argument("--degs",
                        type=str,
                        default=None,
                        help="Polynomial degrees in the shape 'pt,rho' e.g. '2,2'")
    parser.add_argument('-d', '--savedir', type=str, default=None, help="Plots out")
    parser.add_argument('--year', type=str, default=None, help="Year to display on plots.")
    args = parser.parse_args()

    configs = json.load(open(os.path.join(BASE, args.dname, "config.json")))
    if args.year is None:
        args.year = str(configs['year'])
    if args.degs is None:
        args.degs = str(configs['degs'])
    ref_pt, ref_rho = tuple([int(s) for s in args.degs.split(',')])
    nbins = int(configs['NBINS']) + int(configs['NBINSMU'])

    sname = "{}_{}_{}".format(os.path.join(args.savedir, args.dname),
                              str(ref_pt) + str(ref_rho),
                              str(ref_pt + 1) + str(ref_rho),
                              )
    fplot(os.path.join(BASE, args.dname),
          (ref_pt, ref_rho), (ref_pt + 1, ref_rho), args.year, sname)
    sname = "{}_{}_{}".format(os.path.join(args.savedir, args.dname),
                              str(ref_pt) + str(ref_rho),
                              str(ref_pt) + str(ref_rho + 1),
                              )
    fplot(os.path.join(BASE, args.dname),
          (ref_pt, ref_rho), (ref_pt, ref_rho + 1), args.year, sname)
    sname = "{}_{}_{}".format(os.path.join(args.savedir, args.dname),
                              str(ref_pt) + str(ref_rho),
                              str(ref_pt + 1) + str(ref_rho + 1),
                              )
    fplot(os.path.join(BASE, args.dname),
          (ref_pt, ref_rho), (ref_pt + 1, ref_rho + 1), args.year, sname)
