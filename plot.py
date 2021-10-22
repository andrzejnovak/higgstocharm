from __future__ import print_function
import argparse
import os
from collections import OrderedDict
import json

import uproot

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

from util import make_dirs

from utils.plot_fractions import plot_fractions
from rhalphalib.plot.plot_cov import plot_cov

import mplhep as hep
plt.style.use([hep.cms.style.ROOT, {'font.size': 24}])
plt.switch_backend('agg')

pbins = [450, 500, 550, 600, 675, 800, 1200]

cdict = {
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
    # 'vvqq': '#DEC1FF',
    'vvqq': '#C5C9A4'
}

# Sequence of tuples because python2 is stupid
label_dict = OrderedDict([
    ('Data', 'Data'),
    ('MC', 'MC'),
    ('Toys', 'PostFit\nToys'),
    ('hbb', "$\mathrm{H(b\\bar{b})}$"),
    ('hqq', "$\mathrm{H(b\\bar{b})}$"),
    ('zbb', "$\mathrm{Z(b\\bar{b})}$"),
    ('zcc', "$\mathrm{Z(c\\bar{c})}$"),
    ('zqq', "$\mathrm{Z(q\\bar{q})}$"),
    ('hcc', "$\mathrm{H(c\\bar{c})}$"),
    ('wcq', "$\mathrm{W(c\\bar{q})}$"),
    ('wqq', "$\mathrm{W(q\\bar{q})}$"),
    ('vvqq', "$\mathrm{VV}$"),
    ('top', "$\mathrm{Top}$"),
    ('tqq', "$\mathrm{t\\bar{t}}$"),
    ('zll', "$\mathrm{DY}$"),
    ('vlep', "$\mathrm{V(lep.)}$"),
    ('stqq', "$\mathrm{single-t}$"),
    ('qcd', "QCD"),
])

mergedict = {'top': ['stqq', 'tqq'],
             'vlep': ['zll', 'wln'],
             'hcc': ['hcc', 'vbfhcc', 'whcc', 'zhcc'],
             'hbb': ['hbb', 'vbfhbb', 'whbb', 'zhbb', 'tthbb']
}


def full_plot(
        cats,
        pseudo=True,
        fittype="",
        mask=False,
        toys=False,
        sqrtnerr=True,
        filled=True,
        format='png',
        scaleH=True,
        run2=False,
        year='',
):

    # Determine:
    if "pass" in str(cats[0].name) or "fail" in str(cats[0].name):
        regs = "pf"
    elif "pqq" in str(cats[0].name) or "pcc" in str(cats[0].name) or "pbb" in str(
            cats[0].name):
        regs = "3reg"
    else:
        print("Unknown regions")
        return

    # For masking 0 bins (don't want to show them)
    class Ugh():
        def __init__(self):
            self.plot_bins = None

    ugh = Ugh()

    def tgasym_to_err(tgasym):
        # https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/wiki/nonstandard
        # Rescale density by binwidth for actual value
        _binwidth = tgasym._fEXlow + tgasym._fEXhigh
        _x = tgasym._fX
        _y = tgasym._fY * _binwidth
        _xerrlo, _xerrhi = tgasym._fEXlow, tgasym._fEXhigh
        _yerrlo, _yerrhi = tgasym._fEYlow * _binwidth, tgasym._fEYhigh * _binwidth
        return _x, _y, [_yerrlo, _yerrhi], [_xerrlo, _xerrhi]

    def plot_data(x, y, yerr, xerr, ax=None, pseudo=pseudo, ugh=None, **kwargs):
        if ugh is None:
            ugh = Ugh()
        data_err_opts = {
            'linestyle': 'none',
            'marker': '.',
            'markersize': 12.,
            'color': 'k',
            'elinewidth': 2,
        }
        for k in kwargs:
            data_err_opts.setdefault(k, kwargs[k])
        if np.sum([y != 0][0]) > 0:
            if ugh.plot_bins is None:
                ugh.plot_bins = [y != 0][0]
            else:
                ugh.plot_bins = (ugh.plot_bins & [y != 0][0])

        x = np.array(x)[ugh.plot_bins]
        y = np.array(y)[ugh.plot_bins]

        yerr = [np.array(yerr[0])[ugh.plot_bins], np.array(yerr[1])[ugh.plot_bins]]
        xerr = [np.array(xerr)[0][ugh.plot_bins], np.array(xerr)[1][ugh.plot_bins]]

        if mask and not pseudo:
            _y = y
            _y[10:14] = np.nan
            _y[6:9] = np.nan
        else:
            _y = y

        _d_label = "MC" if pseudo else "Data"
        if toys: _d_label = "Toys"
        ax.errorbar(x, y, yerr, xerr, fmt='+', label=_d_label, **data_err_opts)

    def th1_to_step(th1, restoreNorm=True):
        _h, _bins = th1.numpy()
        if restoreNorm:
            _h = _h * np.diff(_bins)
        return _bins, np.r_[_h, _h[-1]], th1.variances

    def th1_to_err(th1, restoreNorm=True):
        _h, _bins = th1.numpy()
        _x = _bins[:-1] + np.diff(_bins) / 2
        _xerr = [abs(_bins[:-1] - _x), _bins[1:] - _x]
        _var = th1.variances
        if restoreNorm:
            _h = _h * np.diff(_bins)
            _var = _var * np.diff(_bins)

        return _x, _h, _var, [_xerr[0], _xerr[1]]

    def plot_step(bins, h, ax=None, label=None, nozeros=True, **kwargs):
        ax.step(bins, h, where='post', label=label, c=cdict[label], lw=2, **kwargs)

    def plot_filled(bins, h, h0=0, ax=None, label=None, nozeros=True, **kwargs):
        if np.sum(h0) == 0:
            h0 = np.zeros_like(h)
        if 'hatch' not in kwargs and 'color' not in kwargs:
            kwargs['color'] = cdict[label]
        else:
            kwargs['edgecolor'] = cdict[label]
        ax.fill_between(bins, h, h0, step='post', label=label, **kwargs)

    def from_cats(fcn, name):
        out = []
        if name in mergedict.keys():
            samples = mergedict[name]
        else:
            samples = [name]
        for _name in samples:
            _count_missing = 0
            for cat in cats:
                try:
                    out.append(fcn(cat[_name]))
                except:
                    _count_missing += 1
            if _count_missing > 0:
                print('Missing {} in {}/{} categories'.format(_name, _count_missing, len(cats)))
        return np.array(out)

    # Sample proofing
    by_cat_samples = []
    for _cat in cats:
        cat_samples = [
            k.decode(encoding="utf-8").split(';')[0] for k in _cat.keys()
            if b'total' not in k
        ]
        by_cat_samples.append(cat_samples)

    # Plotting
    fig, (ax, rax) = plt.subplots(2,
                                  1,
                                  gridspec_kw={'height_ratios': (3, 1)},
                                  sharex=True)
    plt.subplots_adjust(hspace=0)

    #  Main
    if hasattr(cats[0]['data'], "_fEXhigh"):
        res = np.array(list(map(tgasym_to_err, [cat['data'] for cat in cats])))
    else:
        res = np.array(list(map(th1_to_err, [cat['data'] for cat in cats])))
    _x, _h = res[:, 0][0], np.sum(res[:, 1], axis=0)
    _xerr = res[:, -1][0]
    if sqrtnerr:
        _yerr = np.sqrt(_h)
        plot_data(_x, _h, yerr=[_yerr, _yerr], xerr=_xerr, ax=ax, ugh=ugh)
    else:
        # FIXME
        # _yerrlo = np.sqrt(np.sum(res[:, 2][0]**2, axis=0))
        # _yerrhi = np.sqrt(np.sum(res[:, 2][1]**2, axis=0))
        # print(_yerrlo, _yerrhi)
        #plot_data(_x, _h, yerr=[_yerrlo, _yerrhi], xerr=_xerr, ax=ax, ugh=ugh)
        plot_data(_x, _h, yerr=[np.sqrt(_h), np.sqrt(_h)], xerr=_xerr, ax=ax, ugh=ugh)

    # Unc
    res = from_cats(th1_to_step, 'total')

    # Stack qcd/ttbar
    tot_h = None
    for mc, zo in zip(['qcd'], [1]):
        res = from_cats(th1_to_step, mc)
        if len(res) == 0:
            tot_h = np.zeros(len(_x)+1)
            continue
        bins, h = res[:, 0][0], np.sum(res[:, 1], axis=0)
        if tot_h is None:
            plot_step(bins, h, ax=ax, label=mc, zorder=zo)
            tot_h = h
        else:
            plot_step(bins, h + tot_h, label=mc, ax=ax, zorder=zo)
            tot_h += h

    for mc in ['vlep', 'top', 'vvqq', 'wcq', 'wqq']:
        res = from_cats(th1_to_step, mc)
        if len(res) == 0:
            continue
        bins, h = res[:, 0][0], np.sum(np.nan_to_num(res[:, 1], 0), axis=0)
        plot_filled(bins,
                    h + tot_h,
                    h0=tot_h,
                    ax=ax,
                    label=mc,
                    facecolor='white',
                    hatch='///')
        tot_h += h

    # Separate scaled signal
    if scaleH:
        for mc in ['hcc']:
            res = from_cats(th1_to_step, mc)
            if len(res) == 0:
                continue
            bins, h = res[:, 0][0], np.sum(res[:, 1], axis=0)
            if np.sum(abs(h)) > 30:
                scaleH = False
            else:
                h = h * 100
            plot_step(bins, h, ax=ax, label=mc, linestyle='--')

    # Stack plots
    tot_h, bins = None, None
    #stack_samples = ['zcc', 'zbb', 'zqq', 'wcq', 'wqq']
    stack_samples = ['hbb', 'zcc', 'zbb', 'zqq']
    if not scaleH:
        stack_samples = ['hcc'] + stack_samples
    for mc in stack_samples:
        res = from_cats(th1_to_step, mc)
        if len(res) == 0:
            continue
        bins, h = res[:, 0][0], np.sum(res[:, 1], axis=0)
        if tot_h is None:
            if filled:
                plot_filled(bins, h, h0=0, ax=ax, label=mc)
            else:
                plot_step(bins, h, ax=ax, label=mc)
            tot_h = h

        else:
            if filled:
                plot_filled(bins, h + tot_h, h0=tot_h, ax=ax, label=mc)
            else:
                plot_step(bins, h + tot_h, label=mc, ax=ax)
            tot_h += h


    #######
    # Ratio plot
    rax.axhline(0, c='gray', ls='--')

    # Caculate diff
    if hasattr(cats[0]['data'], "_fEXhigh"):
        res = np.array(list(map(tgasym_to_err, [cat['data'] for cat in cats])))
        variance = np.sum([((res[:, 2][i][0] + res[:, 2][i][1])/2)**2 for i in range(len(res[:, 2]))], axis=0)
    else:
        res = np.array(list(map(th1_to_err, [cat['data'] for cat in cats])))
        variance = np.sum(res[:, 1], axis=0)
    _x, _y = res[:, 0][0], np.sum(res[:, 1], axis=0)
    _xerr = res[:, -1][0]
    _yerr = np.sqrt(variance)
    _yerr += 0.00000001  # pad zeros

    # y = np.copy(_y)
    bkgs = np.zeros_like(_y)
    bkg_vars = np.zeros_like(_y)
    for mc in ['qcd', 'top', 'vvqq', 'wcq', 'wqq', 'zbb', 'zqq', 'hbb', 'vlep']:
        res = from_cats(th1_to_step, mc)
        if len(res) == 0:
            continue
        bins, h = res[:, 0][0], np.sum(res[:, 1], axis=0)
        bkgs += h[:-1]
        bkg_vars += np.sum(from_cats(th1_to_err, mc)[:, 2], axis=0) # Should this be squared?

    y = _y - bkgs

    err =_yerr/_yerr
    y /= _yerr
    _scale_for_mc = np.r_[_yerr, _yerr[-1]]
    plot_data(_x, y, yerr=[err, err], xerr=_xerr, ax=rax, ugh=ugh, zorder=10)

    # print(_yerr)
    # print(np.sqrt(bkg_vars))
    band_size = np.sqrt(bkg_vars)/_yerr
    # print(band_size)
    band_size = np.r_[band_size, band_size[-1]]
    rax.fill_between(bins, band_size, -band_size, label=None, color='gray', alpha=0.2, hatch='xxx',
                     step='post')

    # Stack plots
    tot_h, bins = None, None
    #stack_samples = ['hbb', 'zbb', 'zcc', 'zqq', 'wcq', 'wqq']
    #stack_samples = ['hbb', 'zcc', 'zbb', 'zqq', ]
    stack_samples = ['zcc']
    if not scaleH:
        stack_samples = ['hcc'] + stack_samples
    for mc in stack_samples:
        res = from_cats(th1_to_step, mc)
        if len(res) == 0:
            continue
        bins, h = res[:, 0][0], np.sum(res[:, 1], axis=0)
        if tot_h is None:
            if filled:
                plot_filled(bins, h / _scale_for_mc, ax=rax, label=mc)
            else:
                plot_step(bins, h / _scale_for_mc, ax=rax, label=mc)
            tot_h = h
        else:
            if filled:
                plot_filled(bins, (h + tot_h) / _scale_for_mc,
                            tot_h / _scale_for_mc,
                            ax=rax,
                            label=mc)
            else:
                plot_step(bins, (h + tot_h) / _scale_for_mc, label=mc, ax=rax)
            tot_h += h

    # Separate scaled signal
    if scaleH:
        for mc in ['hcc']:
            res = from_cats(th1_to_step, mc)
            if len(res) == 0:
                continue
            bins, h = res[:, 0][0], np.sum(res[:, 1], axis=0)
            plot_step(bins, 100 * h / _scale_for_mc, ax=rax, label=mc, linestyle='--')

    ############
    # Style
    lumi = {
        "jet": {
            "2016": 35.5,
            "2017": 41.5,
            "2018": 59.2,
        },
        "mu": {
            "2016": 35.2,
            "2017": 41.1,
            "2018": 59.0,
        }
    }
    if b'muon' in cats[0].name:
        lumi_t = "mu"
    else:
        lumi_t = "jet"
    if run2:
        ax = hep.cms.cmslabel(ax=ax,
                              data=((not pseudo) | toys),
                              year='',
                              lumi=np.sum([float(v) for k, v in lumi['jet'].items()]),
                              fontsize=22)
    else:
        ax = hep.cms.cmslabel(ax=ax,
                              data=((not pseudo) | toys),
                              year=year,
                              lumi=lumi[lumi_t][str(year)],
                              fontsize=22)
    ax.legend(ncol=2)

    ax.set_ylabel('Events / 7GeV', ha='right', y=1)
    rax.set_xlabel('jet $\mathrm{m_{SD}}$ [GeV]', ha='right', x=1)
    rax.set_ylabel(
        #r'$\mathrm{\frac{Data-(MultiJet+t\bar{t})}{\sigma_{Data}}}$')
        r'$\mathrm{\frac{Data-Bkg}{\sigma_{Data}}}$')

    if b'muon' in cats[0].name:
        ax.set_xlim(0, 1)
        rax.set_xticks([0, 1])
        rax.set_xticklabels([40, 200])
    else:
        ax.set_xlim(40, 200)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.4)
    # ax.ticklabel_format(axis='y', style='sci', scilimits=(0,3), useOffset=False)
    # ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
    f = mtick.ScalarFormatter(useOffset=False, useMathText=True)

    # g = lambda x, pos: "${}$".format(f._formatSciNotation('%1.10e' % x))

    def g(x, pos):
        return "${}$".format(f._formatSciNotation('%1.10e' % x))

    ax.yaxis.set_major_formatter(mtick.FuncFormatter(g))
    rax.set_ylim(rax.get_ylim()[0] * 1.3, rax.get_ylim()[1] * 1.3)

    ipt = int(str(cats[0].name).split('ptbin')[1][0]) if b'ptbin' in cats[0].name else 0
    if len(cats) == 1:
        pt_range = str(pbins[ipt]) + "$< \mathrm{p_T} <$" + str(pbins[ipt + 1]) + " GeV"
    else:
        pt_range = str(pbins[0]) + "$< \mathrm{p_T} <$" + str(pbins[-1]) + " GeV"
    if b'muon' in cats[0].name:
        pt_range = str(pbins[0]) + "$< \mathrm{p_T} <$" + str(pbins[-1]) + " GeV"

    lab_mu = ", MuonCR" if b'muon' in cats[0].name else ""
    if regs == "pf":
        lab_reg = "Passing" if "pass" in str(cats[0].name) else "Failing"
    else:
        if "pqq" in str(cats[0].name):
            lab_reg = "Light"
        elif "pcc" in str(cats[0].name):
            lab_reg = "Charm"
        elif "pbb" in str(cats[0].name):
            lab_reg = "Bottom"

    annot = pt_range + '\nDeepDoubleX{}'.format(lab_mu) + '\n{} Region'.format(lab_reg)

    ax.annotate(annot,
                linespacing=1.7,
                xy=(0.04, 0.94),
                xycoords='axes fraction',
                ha='left',
                va='top',
                ma='center',
                fontsize='small',
                bbox={
                    'facecolor': 'white',
                    'edgecolor': 'white',
                    'alpha': 0,
                    'pad': 13
                },
                annotation_clip=False)

    # Leg sort
    if scaleH:
        label_dict['hcc'] = "$\mathrm{H(c\\bar{c})}$ x 100"
        #label_dict['hqq'] = "$\mathrm{H(b\\bar{b})}$ x 500"
        #label_dict['hbb'] = "$\mathrm{H(b\\bar{b})}$ x 500"
    else:
        label_dict['hcc'] = "$\mathrm{H(c\\bar{c})}$"

    sorted_handles_labels = hep.plot.sort_legend(ax, label_dict)
    # Insert dummy to uneven legend to align right
    if len(sorted_handles_labels[0]) % 2 != 0:
        _insert_ix = len(sorted_handles_labels[0]) / 2
        sorted_handles_labels[0].insert(
            _insert_ix, plt.Line2D([], [], linestyle='none', marker=None))
        sorted_handles_labels[1].insert(_insert_ix, '')
    leg = ax.legend(*sorted_handles_labels, ncol=2, columnspacing=0.8)
    if fittype == 'fit_s':
        fittype = 'postfit'
    leg.set_title(title=fittype.capitalize(), prop={'size': "smaller"})

    if b'muon' in cats[0].name:
        _iptname = "_MuonCR"
    else:
        _iptname = str(str(ipt) if len(cats) == 1 else "")
    # name = str("pass" if "pass" in str(cats[0].name) else "fail"
    #            ) + _iptname
    name = str(lab_reg) + _iptname

    if args:
        fig.savefig('{}/{}.{}'.format(args.output_folder, fittype + "_" + name, format),
                bbox_inches="tight")


if __name__ == '__main__':
    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')


    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", default='', help="Model/Fit dir")
    parser.add_argument("-i",
                        "--input",
                        default='fitDiagnostics.root',
                        help="Input shapes file")
    parser.add_argument("--fit",
                        default=None,
                        choices={"prefit", "fit_s"},
                        dest='fit',
                        help="Shapes to plot")
    parser.add_argument("--3reg",
                        action='store_true',
                        dest='three_regions',
                        help="By default plots pass/fail region. Set to plot pqq/pcc/pbb")
    parser.add_argument("--unmask",
                        action='store_false',
                        dest='mask',
                        help="Mask Higgs bins")
    parser.add_argument("--all",
                        action='store_true',
                        dest='run_all',
                        help="Include split pT bin plots")
    parser.add_argument("--run2", action='store_true', dest='run2', help="Stack all years")
    parser.add_argument("-o",
                        "--output-folder",
                        default='plots',
                        dest='output_folder',
                        help="Folder to store plots - will be created if it doesn't exist.")
    parser.add_argument("--year",
                        default=None,
                        choices={"2016", "2017", "2018"},
                        type=str,
                        help="year label")

    parser.add_argument("--scaleH",
                        type=str2bool,
                        default='True',
                        choices={True, False},
                        help="Scale Higgs signal in plots by 100")

    parser.add_argument("--filled",
                        type=str2bool,
                        default='True',
                        choices={True, False},
                        help="Use filled stack plots")

    parser.add_argument('-f',
                        "--format",
                        type=str,
                        default='png',
                        choices={'png', 'pdf'},
                        help="Plot format")

    pseudo = parser.add_mutually_exclusive_group(required=True)
    pseudo.add_argument('--data', action='store_false', dest='pseudo')
    pseudo.add_argument('--MC', action='store_true', dest='pseudo')
    pseudo.add_argument('--toys', action='store_true', dest='toys')

    parser.add_argument('--inputonly', action='store_true', dest='inputonly')

    args = parser.parse_args()
    if args.output_folder.split("/")[0] != args.dir:
        args.output_folder = os.path.join(args.dir, args.output_folder)

    if not args.run2:
        configs = json.load(open("config.json"))
        if args.year is None:
            args.year = str(configs['year'])

    make_dirs(args.output_folder)

    if args.fit is None:
        shape_types = ['prefit', 'fit_s']
    else:
        shape_types = [args.fit]
    if args.three_regions:
        regions = ['pqq', 'pcc', 'pbb']
    else:
        regions = ['pass', 'fail']

    # f = uproot.open(os.path.join(args.dir, args.input))
    f = uproot.open(args.input)
    for shape_type in shape_types:
        pbins = [450, 500, 550, 600, 675, 800, 1200]
        if args.inputonly:
            continue
        for region in regions:
            print("Plotting {} region".format(region), shape_type)
            mask = (args.mask &
                    (region == "pass")) | (args.mask &
                                        (region == "pcc")) | (args.mask &
                                                                (region == "pbb"))
            for i in range(0, 6):
                if not args.run_all: continue
                cat_name = 'shapes_{}/ptbin{}{}{};1'.format(shape_type, i, region,
                                                            args.year)
                print([cat_name])
                try:
                    cat = f[cat_name]
                except Exception:
                    raise ValueError("Namespace {} is not available, only following"
                                    "namespaces were found in the file: {}".format(
                                        args.fit, f.keys()))

                fig = full_plot([cat],
                                pseudo=args.pseudo,
                                fittype=shape_type,
                                mask=mask,
                                toys=args.toys,
                                sqrtnerr=True,
                                format=args.format,
                                year=args.year)

            if args.run2:
                cat_list = [
                    f['shapes_{}/ptbin{}{}{};1'.format(shape_type, i, region, '2016')]
                    for i in range(0, 6)
                ]
                cat_list += [
                    f['shapes_{}/ptbin{}{}{};1'.format(shape_type, i, region, '2017')]
                    for i in range(0, 6)
                ]
                cat_list += [
                    f['shapes_{}/ptbin{}{}{};1'.format(shape_type, i, region, '2018')]
                    for i in range(0, 6)
                ]
            else:
                cat_list = [
                    f['shapes_{}/ptbin{}{}{};1'.format(shape_type, i, region, args.year)]
                    for i in range(0, 6)
                ]
            full_plot(cat_list,
                    pseudo=args.pseudo,
                    fittype=shape_type,
                    mask=mask,
                    toys=args.toys,
                    sqrtnerr=True,
                    format=args.format,
                    year="" if args.run2 else args.year,
                    run2=args.run2, 
                    )

            # MuonCR if included
            if any(['muonCR' in s for s in f['shapes_{}'.format(shape_type)].keys()]):
                if args.run2:
                    cat_list = [f['shapes_{}/muonCR{}{};1'.format(shape_type, region, '2016')]]
                    cat_list += [f['shapes_{}/muonCR{}{};1'.format(shape_type, region, '2017')]]
                    cat_list += [f['shapes_{}/muonCR{}{};1'.format(shape_type, region, '2018')]]
                    full_plot(cat_list,
                              args.pseudo,
                              fittype=shape_type,
                              toys=args.toys,
                              sqrtnerr=True,
                              format=args.format,
                              year="",
                              scaleH=args.scaleH,
                              run2=args.run2, 
                              )
                else:
                    cat = f['shapes_{}/muonCR{}{};1'.format(shape_type, region, args.year)]
                    full_plot(
                        [cat],
                        args.pseudo,
                        fittype=shape_type,
                        toys=args.toys,
                        sqrtnerr=True,
                        format=args.format,
                        year=args.year,
                        scaleH=args.scaleH,
                    )
                print("Plotted muCR", region, shape_type)
            else:
                print("Muon region not found")

    ##### Input shape plotter
    # Take sqrt N err for data
    # Mock QCD while unavailable as template in rhalpha
    import os
    from rhalphalib.plot.input_shapes import input_dict_maker

    # try:
    mockd = input_dict_maker(os.getcwd() + ".pkl")
    input_pseudo = True
    if args.toys or not args.pseudo:
        input_pseudo = False
    for shape_type in ["inputs"]:
        pbins = [450, 500, 550, 600, 675, 800, 1200]
        for region in regions:
            print("Plotting inputs", region)
            _mask = not input_pseudo
            mask = (_mask &
                    (region == "pass")) | (_mask &
                                            (region == "pcc")) | (_mask &
                                                                    (region == "pbb"))
            full_plot([
                mockd['ptbin{}{}{}_{}'.format(i, region, args.year, shape_type)]
                for i in range(0, 6)
                        ],
                        pseudo=input_pseudo,
                        fittype=shape_type,
                        mask=mask,
                        sqrtnerr=True,
                        toys=False,
                        format=args.format,
                        year=args.year,
                        scaleH=args.scaleH,
                        )
            # Per bin plots
            for i in range(0, 6):
                if not args.run_all: continue
                full_plot(
                    [mockd['ptbin{}{}{}_{}'.format(i, region, args.year, shape_type)]],
                    pseudo=input_pseudo,
                    fittype=shape_type,
                    mask=mask,
                    sqrtnerr=True,
                    toys=False,
                    format=args.format,
                    year=args.year,
                    scaleH=args.scaleH,
                    )

            # MuonCR if included
            try:
                cat = mockd['muonCR{}{}_{}'.format(region, args.year, shape_type)]
                full_plot([cat],
                            fittype=shape_type,
                            pseudo=input_pseudo,
                            mask=False,
                            sqrtnerr=True,
                            toys=False,
                            format=args.format,
                            year=args.year,
                            scaleH=args.scaleH
                )
                print("Plotted input, muCR", region, shape_type)
            except Exception:
                print("Muon region not found")
                pass
    # except:
    #     print("Input pkl file not found")

    if args.three_regions:
        plot_fractions(
            os.path.join(args.dir, 'fitDiagnostics.root'),
            os.path.join(args.dir, 'model_combined.root'),
            out='{}/{}.{}'.format(args.output_folder, 'fractions', args.format),
            data=((not args.pseudo) | args.toys),
            year=args.year,
        )

    plot_cov(
        os.path.join(args.dir, 'fitDiagnostics.root'),
        out='{}/{}.{}'.format(args.output_folder, 'covariances', args.format),
        data=((not args.pseudo) | args.toys),
        year=args.year,
    )

    plot_cov(
        os.path.join(args.dir, 'fitDiagnostics.root'),
        out='{}/{}_wTF.{}'.format(args.output_folder, 'covariances', args.format),
        data=((not args.pseudo) | args.toys),
        year=args.year,
        include='tf',
    )
