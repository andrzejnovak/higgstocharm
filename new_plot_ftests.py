import os
import argparse
import ROOT as r
import numpy as np
import matplotlib
matplotlib.use("Agg")
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.ROOT)
import matplotlib.transforms as transforms


def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.makedirs(file_path)


def skim_gofs(file_names):
    out = {}
    for i, fname in enumerate(file_names):
        try:
            rfile = r.TFile.Open(fname)
            rtree = rfile.Get("limit")
            for j in range(rtree.GetEntries()):
                rtree.GetEntry(j)
                mu = rtree.limit
                _fname = fname.split("/")[-1].replace("BaseToys", "Toys").replace("AltToys", "Toys")
                out[_fname+":"+str(j)] = mu
        except:
            print("        Skipping: {}".format(fname.split("/")[-1]))
            pass
    return out


def collate_gofs(base_dict, alt_dict):
    bar, altr = [], []
    for key in base_dict.keys():
        bar.append(base_dict[key])
        altr.append(alt_dict[key])
    return np.array(bar), np.array(altr)


@np.vectorize
def fval(lambda1, lambda2, p1, p2, nbins):
    with np.errstate(divide='ignore'):
        numerator = -2.0 * np.log(float(lambda1) / float(lambda2)) / float(p2 - p1)
        denominator = -2.0 * np.log(float(lambda2)) / float(nbins - p2)
        return numerator / denominator


def fgofs(gofs_base, gofs_alt, data_ref, data_alt, ref, alt, year="2016", mc=False, savename=None):
    ref_pt, ref_rho = ref
    alt_pt, alt_rho = alt
    
    fig, ax = plt.subplots()
    ax.hist(gofs_base, alpha=0.5, label="Base - toys", bins=30, color='r')
    ax.hist(gofs_alt, alpha=0.5, label="Alt - toys", bins=30, color='blue');
    ax.axvline(data_ref, label="Base - ref", color='red')
    ax.axvline(data_alt, label="Alt - ref", color='blue')
    ax.legend()

    title = "TF({},{}) x TF({},{})".format(ref_pt, ref_rho, alt_pt, alt_rho)
    ax.legend(title=title, loc="best")
    hep.cms.label(data=not mc, year=year, ax=ax)
    x_lim = max(np.percentile(gofs_base, 90), np.percentile(gofs_alt, 90), max(data_ref, data_alt)*1.05)
    ax.set_xlim(0, x_lim)
    xlab = r"$-2log(\lambda)$"
    ax.set_xlabel(xlab , x=1, ha='right')
    ax.set_ylabel("Pseudoexperiments", y=1, ha='right')

    if savename is not None:
        fig.savefig('{}.pdf'.format(savename), dpi=300, transparent=True, bbox_inches='tight')
        fig.savefig('{}.png'.format(savename), dpi=300, transparent=True, bbox_inches='tight')


def fplot(fvals, f_data, ref, alt, year=2017, nbins=130, savename=None, mc=False):
    ref_pt, ref_rho = ref
    alt_pt, alt_rho = alt
    p1 = (ref_pt + 1) * (ref_rho + 1)
    p2 = (alt_pt + 1) * (alt_rho + 1)
    
    from scipy.stats import f
    x_lim = max(np.percentile(fvals, 95), f_data*1.05, np.median(fvals) * 3)
    x = np.linspace(0, x_lim, 200)
    bins = np.linspace(0, x_lim, 30)
    width = bins[1]-bins[0]
    
    goodvals = fvals[fvals > 0]

    fig, ax = plt.subplots()
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    with np.errstate(divide='ignore'):
        ax.plot(x, len(goodvals)*width*f.pdf(x, p2-p1, nbins-p2), color='red', label='F-dist, ndf({},{})'.format(p2-p1, nbins-p2))
    ax.hist(fvals, bins, facecolor='none', edgecolor='black', histtype='stepfilled', lw=2,
            label="Toys > 0, N = {}".format(len(goodvals)));
    ax.hist(goodvals[goodvals > f_data], bins, facecolor='steelblue', edgecolor='gray', histtype='stepfilled', alpha=0.3,
            label='p-value = {}'.format(round(float(len(goodvals[goodvals > f_data]))/len(goodvals), 3)));
    ax.annotate("", xy=(f_data, 0), xycoords=trans,
                xytext=(f_data, 0.25), textcoords=trans,
                arrowprops=dict(lw='4', color='b', arrowstyle="->,head_length=1.5,head_width=0.5"),
               )
    ax.plot([], [], color='blue', lw=2, label="Observed = {:.3f}".format(f_data))

    title = "TF({},{}) x TF({},{})".format(ref_pt, ref_rho, alt_pt, alt_rho)
    ax.legend(title=title)
    hep.cms.label(data=not mc, year=year, ax=ax)
    ax.set_xlim(0, x_lim);
    xlab = r"$\frac{-2log(\lambda_1/\lambda_2)/(p_2-p_1)}{-2log\lambda_2/(n-p_2)}$"
    ax.set_xlabel(xlab , x=1, ha='right')
    ax.set_ylabel("Pseudoexperiments", y=1, ha='right')

    if savename is not None: 
        fig.savefig('{}.pdf'.format(savename), dpi=300, transparent=True, bbox_inches='tight')
        fig.savefig('{}.png'.format(savename), dpi=300, transparent=True, bbox_inches='tight')
    

def prep_data(infile, reg):
    prefit_qcd = infile['shapes_prefit/{}/qcd'.format(reg)].values
    postfit_qcd = infile['shapes_fit_s/{}/qcd'.format(reg)].values
    data = infile['shapes_fit_s/{}/data'.format(reg)].yvalues
    bins = infile['shapes_prefit/{}/qcd'.format(reg)].edges
    centers = infile['shapes_fit_s/{}/data'.format(reg)].xvalues
    return data, prefit_qcd, postfit_qcd, bins, centers


def prepostplot(data, prefit_qcd, postfit_qcd, bins, centers, degs=(1,2), year='2016', chi2=None, reg=None, savename=None):
    width = (bins[1] - bins[0]) / 2
    
    fig, (ax, subax) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
    fig.subplots_adjust(hspace=0)
    ax.set_xlim(40, 201)

    if chi2 is None:
        with np.errstate(divide='ignore'):
            chi2_pre = np.sum(np.nan_to_num((data-prefit_qcd)**2/prefit_qcd, 0))
            chi2_post = np.sum(np.nan_to_num((data-postfit_qcd)**2/prefit_qcd, 0))
    else:
        chi2_pre, chi2_post = chi2
    
    # Plot
    hep.histplot(prefit_qcd, bins, color='red', linestyle="--", ax=ax,
                 label='PreFit, $\chi^2 = {:.1f}$'.format(chi2_pre))
    hep.histplot(postfit_qcd, bins, color='blue', ax=ax, 
                 label='PostFit, $\chi^2 = {:.1f}$'.format(chi2_post))
    ax.errorbar(centers, data, fmt='o', xerr=width, color='black', label='True QCD')
    
    leg = ax.legend(title=(reg+"\n" if reg is not None else "") + "TF({},{})".format(degs[0], degs[1]))
    plt.setp(leg.get_title(), multialignment='center')
    hep.cms.label(data=False, year=year, ax=ax)
    ax.set_ylabel("Events / 7GeV", y=1, ha='right')
    ax.set_yticks(ax.get_yticks()[1:])

    #Subplots
    subax.axhline(1, color='grey', ls='--')
    subax.errorbar(centers, (postfit_qcd/data),
                   fmt='o', xerr=width, color='blue')
    subax.axhline(1, color='grey', ls='--')
    eb = subax.errorbar(centers, (prefit_qcd/data),
                   fmt='o', xerr=width, color='red')
    eb[-1][0].set_linestyle('--')
    subax.set_ylim(0.5, 1.5)
    subax.set_ylabel('Pred/True')
    subax.set_xlabel('jet $m_{SD}$ [GeV]', x=1, ha='right')

    if savename is not None: 
        fig.savefig('{}.pdf'.format(savename), dpi=300, transparent=True, bbox_inches='tight')
        fig.savefig('{}.png'.format(savename), dpi=300, transparent=True, bbox_inches='tight')

def prepostplotall(infile, year="2016", degs=(2,4), savename=None):
    col = []
    for i in range(6):
        reg = 'ptbin'+str(i)+'pass'+str(year)
        col.append(prep_data(infile, reg))

    data = np.sum(np.array(col)[:, 0], axis=0)
    prefit_qcd = np.sum(np.array(col)[:, 1], axis=0)
    postfit_qcd = np.sum(np.array(col)[:, 2], axis=0)
    bins = np.array(col)[0, 3]
    centers = np.array(col)[0, 4]

    with np.errstate(divide='ignore'):
        chi2_pre = np.sum(np.nan_to_num(np.vstack((np.array(col)[:, 0] - np.array(col)[:, 1])**2 / np.array(col)[:, 1])))
        chi2_post = np.sum(np.nan_to_num(np.vstack((np.array(col)[:, 0] - np.array(col)[:, 2])**2 / np.array(col)[:, 2])))

    prepostplot(data, prefit_qcd, postfit_qcd, bins, centers, degs=degs, chi2=(chi2_pre, chi2_post), year=year, savename=savename)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make plots for ftest/gof files.')
    # parser.add_argument('rundir', nargs='?', default=os.getcwd())
    parser.add_argument('-d', '--rundir', type=str, required=True, help="Run dir")
    parser.add_argument("--degs",
                        type=str,
                        default=None,
                        help="Polynomial degrees in the shape 'pt,rho' e.g. '2,2'")
    parser.add_argument("--degsalt",
                        type=str,
                        default=None,
                        help="Polynomial degrees in the shape 'pt,rho' e.g. '2,2'")
    parser.add_argument('-o', '--outdir', type=str, default=None, help="Plots out")
    parser.add_argument('--year', type=str, choices=["2016", "2017", "2018"], required=True, help="Year to display on plots.")
    parser.add_argument('--mc', action='store_true', dest='mc')
    parser.add_argument('--qplots', action='store_true')
    parser.add_argument('--nbins', type=int, default=130, help="NBINS")
    args = parser.parse_args()


    ref_pt, ref_rho = tuple([int(s) for s in args.degs.split(',')])
    alt_pt, alt_rho = tuple([int(s) for s in args.degsalt.split(',')])
    p1 = (ref_pt + 1) * (ref_rho + 1)
    p2 = (alt_pt + 1) * (alt_rho + 1)

    dirname = args.rundir
    if args.outdir is None:
        args.outdir = "plots"
    ensure_dir(args.outdir)

    # Extract
    names = sorted([k for k in os.listdir(dirname) if k.endswith('root') and "BaseToys.Good" in k])
    names_alt = sorted([k for k in os.listdir(dirname) if k.endswith('root') and "AltToys.Good" in k])

    base_dict = skim_gofs(sorted([os.path.join(dirname, n) for n in names]))
    alt_dict = skim_gofs(sorted([os.path.join(dirname, n) for n in names_alt]))
    gofs_base, gofs_alt = collate_gofs(base_dict, alt_dict)

    data_ref = skim_gofs([os.path.join(dirname, 'higgsCombine.Base.GoodnessOfFit.mH120.root')]).values()[0]
    data_alt = skim_gofs([os.path.join(dirname, 'higgsCombine.Alt.GoodnessOfFit.mH120.root')]).values()[0]


    # Plot
    savepath = os.path.join(args.outdir, "fplot_gofs_{}{}_{}{}".format(ref_pt, ref_rho, alt_pt, alt_rho))
    fgofs(gofs_base, gofs_alt, data_ref, data_alt, (ref_pt, ref_rho), (alt_pt, alt_rho), mc=args.mc, year=args.year, savename=savepath)

    f_vals = fval(gofs_base, gofs_alt, p1, p2, args.nbins)
    f_ref = fval(data_ref, data_alt, p1, p2, args.nbins)

    savepath = os.path.join(args.outdir, "fplot_{}{}_{}{}".format(ref_pt, ref_rho, alt_pt, alt_rho))
    fplot(f_vals, f_ref, (ref_pt, ref_rho), (alt_pt, alt_rho), args.year, mc=args.mc, savename=savepath)


    if args.qplots:
        ensure_dir(os.path.join(args.outdir, "qcd_{}{}".format(ref_pt, ref_rho)))
        ensure_dir(os.path.join(args.outdir, "qcd_{}{}".format(alt_pt, alt_rho)))
    
        fi_base = uproot.open(os.path.join(dirname,'fitDiagnostics.Base.root'))
        fi_alt = uproot.open(os.path.join(dirname,'fitDiagnostics.Alt.root'))

        for i in range(6):
            reg = 'ptbin'+str(i)+'pass'+str(args.year)
            savepath = os.path.join(args.outdir, "qcd_{}{}".format(ref_pt, ref_rho), "qcd_pt{}".format(i))
            prepostplot(*prep_data(fi_base, reg), degs=(ref_pt, ref_rho), reg=reg, year=args.year, savename=savepath)   
        savepath = os.path.join(args.outdir, "qcd_{}{}".format(ref_pt, ref_rho), "qcd_allpt") 
        prepostplotall(fi_base, savename=savepath, degs=(ref_pt, ref_rho), year=args.year)

        for i in range(6):
            reg = 'ptbin'+str(i)+'pass'+str(args.year)
            savepath = os.path.join(args.outdir, "qcd_{}{}".format(alt_pt, alt_rho), "qcd_pt{}".format(i))
            prepostplot(*prep_data(fi_alt, reg), degs=(alt_pt, alt_rho), reg=reg, year=args.year, savename=savepath)   
        savepath = os.path.join(args.outdir, "qcd_{}{}".format(alt_pt, alt_rho), "qcd_allpt") 
        prepostplotall(fi_alt, savename=savepath, degs=(alt_pt, alt_rho), year=args.year)
