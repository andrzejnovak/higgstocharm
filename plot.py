from __future__ import print_function
import argparse
import os
from collections import OrderedDict
import json, copy
from re import sub

import uproot
import ROOT as r

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

import logging

from rhalphalib.plot.plot_cov import plot_cov

from util import make_dirs
from utils.plot_fractions import plot_fractions

from plot_utils import format_legend


import mplhep as hep
plt.style.use([hep.cms.style.ROOT, {'font.size': 24}])
plt.switch_backend('agg')


from plot_styles import scale_lightness, colorize, merge_dict
from plot_styles import style_set_A, style_set_B
from plot_utils import format_legend


def get_fit_val(fitDiag, val, fittype='fit_s', substitute=1.):
    if fitDiag is None:
        return substitute
    if val in fitDiag.Get(fittype).floatParsFinal().contentsString().split(','):
        return fitDiag.Get(fittype).floatParsFinal().find(val).getVal()
    else:
        return substitute


def lite_plot(
        cats,
        pseudo=True,
        fittype="",
        mask=False,
        toys=False,
        filled=True,  # 
        ratio_samples=['zcc', 'hcc'],
        scaleH=None,
        stack_style=0,
        style_set=style_set_A,
        prelim=True,
        format='png',
        fitDiag=None,  # Add signal strength labels
        bkgUnc=None,  # External bkg uncertainty
        args=None,  # Configs
        run2=False,  # For label
        paper=False,  
        pbins=[450, 500, 550, 600, 675, 800, 1200],  # For label
        year='',  # For label
):
    color_dict = style_set['color_dict']
    hatch_dict = style_set['hatch_dict']
    label_dict = style_set['label_dict']

    # Determine:
    if "pass" in str(cats[0].name) or "fail" in str(cats[0].name):
        regs = "pf"
    elif "pqq" in str(cats[0].name) or "pcc" in str(cats[0].name) or "pbb" in str(
            cats[0].name):
        regs = "3reg"
    else:
        print("Unknown regions")
        return

    if fitDiag is not None and not b'muon' in cats[0].name:
        _r = get_fit_val(fitDiag, 'r', fittype='fit_s', substitute=1.),
        _z = get_fit_val(fitDiag, 'z', fittype='fit_s', substitute=1.)
    if fittype == 'prefit':
        _r = 1.
        _z = 1.

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

    def plot_data(x,
                  y,
                  yerr,
                  xerr,
                  ax=None,
                  pseudo=pseudo,
                  label=None,
                  ugh=None,
                  **kwargs):
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
            # _y[6:9] = np.nan
        else:
            _y = y

        if label is None:
            _d_label = "MC" if pseudo else "Data"
            if toys: _d_label = "Toys"
        elif label == "":
            _d_label = None
        else:
            _d_label = label
        ax.errorbar(x, y, yerr, xerr, fmt='+', label=_d_label, **data_err_opts)

    def th1_to_step(th1, restoreNorm=True):
        _h, _bins = th1.numpy()
        _variances = th1.variances
        if restoreNorm:
            _h = _h * np.diff(_bins)
            _variances = _variances * np.diff(_bins)
        return _bins, np.r_[_h, _h[-1]], _variances

    def th1_to_err(th1, restoreNorm=True):
        _h, _bins = th1.numpy()
        _x = _bins[:-1] + np.diff(_bins) / 2
        _xerr = [abs(_bins[:-1] - _x), _bins[1:] - _x]
        _var = th1.variances
        if restoreNorm:
            _h = _h * np.diff(_bins)
            _var = _var * np.diff(_bins)

        return _x, _h, _var, [_xerr[0], _xerr[1]]

    def plot_step(bins, h, ax=None, sample=None, label=None, 
                  nozeros=True, linewidth=2, **kwargs):
        try:
            ax.step(bins,
                    h,
                    where='post',
                    label=label,
                    c=color_dict[sample],
                    linewidth=linewidth,
                    **kwargs)
        except Exception as e:
            print("Some issue with {}, {}".format(label, sample))
            raise e

    def plot_filled(bins,
                    h,
                    h0=0,
                    ax=None,
                    sample=None,
                    label=None,
                    nozeros=True,
                    **kwargs):
        if np.sum(h0) == 0:
            h0 = np.zeros_like(h)
        if 'hatch' not in kwargs and 'color' not in kwargs:
            kwargs['facecolor'] = color_dict[sample]
        else:
            kwargs['edgecolor'] = color_dict[sample]
        if 'hatch' in kwargs and kwargs['hatch'] is None:
            kwargs['edgecolor'] = "none"
        if np.sum(h) < 1:
            print(label)
            kwargs['edgecolor'] = 'none'
            kwargs['linewidth'] = 0
        ax.fill_between(bins, h, h0, step='post', label=label, **kwargs)

    def from_cats(fcn, name):
        out = []
        if name in merge_dict.keys():
            samples = merge_dict[name]
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
                # logging.warning('Watch out!')
                logging.info('Missing {} in {}/{} categories'.format(
                    _name, _count_missing, len(cats)))
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
    _yerr = np.sqrt(_h)
    plot_data(_x, _h, yerr=[_yerr, _yerr], xerr=_xerr, ax=ax, ugh=ugh)

    # Unc
    res = from_cats(th1_to_step, 'total')

    # Stack qcd/ttbar
    tot_h = None
    for mc, zo in zip(['qcd'], [1]):
        res = from_cats(th1_to_step, mc)
        if len(res) == 0:
            tot_h = np.zeros(len(_x) + 1)
            continue
        bins, h = res[:, 0][0], np.sum(res[:, 1], axis=0)
        if tot_h is None:
            plot_step(bins, h, ax=ax, label=mc, sample=mc, zorder=zo)
            tot_h = h
        else:
            plot_step(bins, h + tot_h, label=mc, sample=mc, ax=ax, zorder=zo)
            tot_h += h

    # Stack peaky stuff
    # 'vlep' not visible
    if stack_style == 0:
        _stack_samples = ['top', 'vlep', 'vvqq', 'wqq', 'wcq']
    elif stack_style == 1:
        _stack_samples = [
            'top', 'vlep', 'vvqq', 'wqq', 'wcq', "zqq", "zbb", "zcc", "hbb", 'hcc'
        ]
    elif stack_style == 2:
        _stack_samples = ['nonV', "wqq", "wcq", "zqq", "zbb", 'zcc', 'hbb', 'hcc']
    elif stack_style == 3:
        _stack_samples = ['bkgNoW', "wcq", 'zcc', 'hbb', 'hcc']
    elif stack_style == 4:
        _stack_samples = ['bkg', 'zcc', 'hcc']
    else:
        _stack_samples = [
            'vlep', 'top', 'vvqq', 'wqq', 'wcq', "zqq", "zbb", "zcc", "hbb", 'hcc'
        ]


#     if scaleH is None:
#         _stack_samples = ['hcc'] + _stack_samples
    for mc in _stack_samples:
        res = from_cats(th1_to_step, mc)
        if len(res) == 0:
            continue
        bins, h = res[:, 0][0], np.sum(np.nan_to_num(res[:, 1], 0), axis=0)
        plot_filled(bins,
                    h + tot_h,
                    h0=tot_h,
                    ax=ax,
                    label=mc,
                    sample=mc,
                    facecolor='white' if hatch_dict[mc] is not None else color_dict[mc],
                    hatch=hatch_dict[mc])
        tot_h += h

    # Separate scaled signal
    if scaleH is not None:
        for mc in ['hcc']:
            res = from_cats(th1_to_step, mc)
            if len(res) == 0:
                continue
            bins, h = res[:, 0][0], np.sum(res[:, 1], axis=0)
            h = h * scaleH / _r
            plot_step(bins, h, ax=ax, label="hccs", sample=mc, linestyle='--', linewidth=3)

    # Stack project onto x-axis
    if stack_style == 0:
        tot_h, bins = None, None
        stack_samples = ['hbb', 'zcc', 'zbb', 'zqq', 'hcc']
        for mc in stack_samples:
            res = from_cats(th1_to_step, mc)
            if len(res) == 0:
                continue
            bins, h = res[:, 0][0], np.sum(res[:, 1], axis=0)
            if tot_h is None:
                if not filled:
                    plot_step(bins, h, ax=ax, label=mc, sample=mc)
                elif filled:
                    plot_filled(bins, h, h0=0, ax=ax, label=mc, sample=mc)
                tot_h = h

            else:
                if not filled:
                    plot_step(bins, h + tot_h, label=mc, sample=mc, ax=ax)
                elif filled:
                    plot_filled(bins, h + tot_h, h0=tot_h, label=mc, ax=ax, sample=mc)
                tot_h += h

    #######
    # Ratio plot
    rax.axhline(0, c='gray', ls='--')

    # Caculate diff
    if hasattr(cats[0]['data'], "_fEXhigh"):
        res = np.array(list(map(tgasym_to_err, [cat['data'] for cat in cats])))
        variance = np.sum([((res[:, 2][i][0] + res[:, 2][i][1]) / 2)**2
                           for i in range(len(res[:, 2]))],
                          axis=0)
    else:
        res = np.array(list(map(th1_to_err, [cat['data'] for cat in cats])))
        variance = np.sum(res[:, 1], axis=0)
    _x, _y = res[:, 0][0], np.sum(res[:, 1], axis=0)
    _xerr = res[:, -1][0]
    _yerr = np.sqrt(variance)
    _yerr += 0.00000001  # pad zeros

    bkgs = np.zeros_like(_y)
    bkg_vars = np.zeros_like(_y)
    MCLIST = ['qcd', 'top', 'vvqq', 'wcq', 'wqq', 'zbb', 'zqq', 'hbb', 'vlep']
    for mc in [mc for mc in MCLIST if mc not in ratio_samples]:
        res = from_cats(th1_to_step, mc)
        if len(res) == 0:
            continue
        bins, h = res[:, 0][0], np.sum(res[:, 1], axis=0)
        bkgs += h[:-1]
        bkg_vars += np.sum(from_cats(th1_to_err, mc)[:, 2],
                           axis=0)  # Should this be squared?

    y = _y - bkgs

    err = _yerr / _yerr
    y /= _yerr
    _scale_for_mc = np.r_[_yerr, _yerr[-1]]
    plot_data(_x, y, yerr=[err, err], xerr=_xerr, ax=rax, ugh=ugh, zorder=10, label="")

    band_size = np.sqrt(bkg_vars) / _yerr
    if bkgUnc is not None:
        band_size = bkgUnc / _yerr
    band_size = np.r_[band_size, band_size[-1]]
    rax.fill_between(bins,
                     band_size,
                     -band_size,
                     color='gray',
                     alpha=0.3,
                     hatch='xxx',
                     step='post',
                     label='Bkg. unc.')
    rax.legend(columnspacing=1, loc='upper right', markerscale=0.9, fontsize='small')

    # Stack plots
    tot_h, bins = None, None
    for mc in ratio_samples:
        res = from_cats(th1_to_step, mc)
        if len(res) == 0:
            continue
        bins, h = res[:, 0][0], np.sum(res[:, 1], axis=0)
        if tot_h is None:
            if filled:
                plot_filled(
                    bins,
                    h / _scale_for_mc,
                    ax=rax,
                    sample=mc,
                    facecolor='white' if hatch_dict[mc] is not None else color_dict[mc],
                    hatch=hatch_dict[mc])
            else:
                plot_step(bins, h / _scale_for_mc, ax=rax, label=mc, sample=mc)
            tot_h = h
        else:
            if filled:
                plot_filled(
                    bins, (h + tot_h) / _scale_for_mc,
                    tot_h / _scale_for_mc,
                    ax=rax,
                    sample=mc,
                    facecolor='white' if hatch_dict[mc] is not None else color_dict[mc],
                    hatch=hatch_dict[mc])
            else:
                plot_step(bins, (h + tot_h) / _scale_for_mc,
                          label=mc,
                          sample=mc,
                          ax=rax)
            tot_h += h

    # Separate scaled signal
    if scaleH is not None:
        for mc in ['hcc']:
            res = from_cats(th1_to_step, mc)
            if len(res) == 0:
                continue
            bins, h = res[:, 0][0], np.sum(res[:, 1], axis=0)
            _scaled_h = scaleH * h / _scale_for_mc / _r
            plot_step(bins, _scaled_h, ax=rax, label=mc, sample=mc, linestyle='--', linewidth=3)

    ############
    # Style
    lumi = {
        "jet": {
            "2016": 36.3,
            "2017": 41.5,
            "2018": 59.8,
        },
        "mu": {
            "2016": 36.3,
            "2017": 41.5,
            "2018": 59.8,
        }
    }
    if b'muon' in cats[0].name:
        lumi_t = "mu"
    else:
        lumi_t = "jet"

    _data = ((not pseudo) | toys)
    _label = "" if _data else "Simulation"
    if prelim:
        _label += " Preliminary"
    else:
        if stack_style not in [2]:  # not paper
            _label += " Supplementary"

    if run2 or year is None or year == "":
        ax = hep.cms.cmslabel(ax=ax, llabel=_label, lumi=138, fontsize=22, year="", paper=paper)
    else:
        ax = hep.cms.cmslabel(ax=ax,
                              llabel=_label,
                              paper=paper,
                              year=year,
                              lumi=lumi[lumi_t][str(year)],
                              fontsize=22)
    ax.legend(ncol=2)

    ax.set_ylabel('Events / 7GeV', ha='right', y=1)
    rax.set_xlabel('Jet $\mathrm{m_{SD}}$ [GeV]', ha='right', x=1)
    rax.set_ylabel(r'$\mathrm{\frac{Data-Bkg}{\sigma_{Data}}}$')

    if b'muon' in cats[0].name:
        ax.set_xlim(0, 1)
        rax.set_xticks([0, 1])
        rax.set_xticklabels([40, 200])
    else:
        ax.set_xlim(40, 200)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.4)
    if len(ax.get_legend_handles_labels()[0]) > 9:
        ax.set_ylim(0, ax.get_ylim()[1] * 1.2)
    f = mtick.ScalarFormatter(useOffset=False, useMathText=True)

    # g = lambda x, pos: "${}$".format(f._formatSciNotation('%1.10e' % x))

    def g(x, pos):
        return "${}$".format(f._formatSciNotation('%1.1e' % x))

    ax.yaxis.set_major_formatter(mtick.FuncFormatter(g))
    rax.set_ylim(rax.get_ylim()[0] * 1.3, rax.get_ylim()[1] * 1.3)

    cat_ixes = [int(str(cat.name).split('ptbin')[1][0]) if b'ptbin' in cat.name else 0 for cat in cats]
    if len(set(cat_ixes)) == 1:
        pt_range = str(pbins[cat_ixes[0]]) + "$< \mathrm{p_T} <$" + str(pbins[cat_ixes[0] + 1]) + " GeV"
    else:
        pt_range = str(pbins[min(cat_ixes)]) + "$< \mathrm{p_T} <$" + str(pbins[max(cat_ixes) + 1]) + " GeV"
    if b'muon' in cats[0].name:
        pt_range = str(pbins[0]) + "$< \mathrm{p_T} <$" + str(pbins[-1]) + " GeV"

    # ipt = int(str(cats[0].name).split('ptbin')[1][0]) if b'ptbin' in cats[0].name else 0
    # if len(cats) == 1:
    #     pt_range = str(pbins[ipt]) + "$< \mathrm{p_T} <$" + str(pbins[ipt + 1]) + " GeV"
    # else:
    #     pt_range = str(pbins[0]) + "$< \mathrm{p_T} <$" + str(pbins[-1]) + " GeV"
    # if b'muon' in cats[0].name:
    #     pt_range = str(pbins[0]) + "$< \mathrm{p_T} <$" + str(pbins[-1]) + " GeV"

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
    if scaleH is not None:
        # label_dict['hccs'] = "$\mathrm{H(c\\bar{c})}_{\mu=%s}$" % scaleH
        label_dict['hccs'] = "$\mathrm{H(c\\bar{c})} \\times %s$" % scaleH
    else:
        label_dict['hccs'] = "$\mathrm{H(c\\bar{c})}$"

    order = np.argsort(
        [list(label_dict.keys()).index(i) for i in ax.get_legend_handles_labels()[-1]])
    handles, labels = ax.get_legend_handles_labels()
    handles = [handles[i] for i in order]
    labels = [label_dict[labels[i]] for i in order]
    leg = format_legend(ax,
                        ncols=2,
                        handles_labels=(handles, labels),
                        bbox_to_anchor=(1, 1),
                        markerscale=0.8, fontsize='small',
                        labelspacing=0.4,
                        columnspacing=1.5)
    if fittype == 'fit_s':
        fit_name = 'postfit'
    else:
        fit_name = copy.deepcopy(fittype)
    if fittype != 'fit_s' or b'muon' in cats[0].name:
        leg.set_title(title=fit_name.capitalize(), prop={'size': "small"})
    _r, _z = 1., 1.
    if fitDiag is not None and (b'muon' not in cats[0].name and 'prefit' not in fittype):
        fig.canvas.draw()
        _r = get_fit_val(fitDiag, 'r', fittype='fit_s', substitute=1.)
        _z = get_fit_val(fitDiag, 'z', fittype='fit_s', substitute=1.)
        # print("X", _r, _z)
        # print(fittype, fit_name, cats[0].name)
        val_str = "$\mu_{Z(c\\bar{c})}$ = " + "{:.1f}".format(_z) + "   " + "$\mu_{H(c\\bar{c})}$ = " + "{:.1f}".format(_r)
        # print(val_str)
        fig.canvas.draw()
        bbox = leg.get_window_extent().inverse_transformed(ax.transAxes)
        x, y = bbox.xmin + bbox.width / 2, bbox.ymin
        ax.text(x,
                y,
                val_str,
                fontsize='small',
                ha='center',
                va='top',
                transform=ax.transAxes)

    if b'muon' in cats[0].name:
        _iptname = "_MuonCR"
    else:
        _iptname = str(str(cat_ixes[0]) if len(cats) == 1 else "")
    name = str(lab_reg) + _iptname

    if args:
        if len(set(cat_ixes)) == 1:
            if not run2: # Single cat
                save_name = [year]
            else: # pt bin across years
                if b'muon' in cats[0].name:
                    save_name = []
                else:
                    save_name = ["pt{}".format(cat_ixes[0])]
        else:
            if not run2:
                save_name = [year]
            else:
                save_name = []
        save_name += [fit_name, name]
        save_name = "_".join(save_name)
        prel_name = "prelim" if prelim else "paper"
        dir_path = '{}/{}/style{}'.format(args.output_folder, prel_name, str(stack_style))
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        fig.savefig('{}/{}/style{}/{}.{}'.format(args.output_folder, prel_name, str(stack_style),
                                             save_name,
                                             format),
                    dpi=300,
                    bbox_inches='tight',
                    transparent=True)
        plt.close(fig)
    return fig, (ax, rax)


def full_plot(*args, **kwargs):
    mod = 'stack_style' not in kwargs
    for st_style in [0, 2, 3]:
        if mod:
            kwargs['stack_style'] = st_style
        lite_plot(*args, **kwargs)


def get_cat_keys(cat, nominal):
    keys = [k for k in nominal.keys() if 'shapeBkg' in k and cat in k]
    return keys

def get_yerr(nominal, varied, keys, ntoys=200):
    nsum = sum(nominal[k][0] for k in keys)
    vsum = sum(varied[k] for k in keys)
    lo, hi = np.quantile(vsum, [0.5 - 0.6827 / 2, 0.5 + 0.6827 / 2], axis=0)
    combinestyle = np.sqrt(np.sum(np.power(vsum - nsum, 2), axis=0) / ntoys)

    mid, width = (hi + lo) / 2, (hi - lo) / 2
    midr, widthr = mid / nsum, width / nsum
    return width


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
                        default='fitDiagnosticsTest.root',
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
    parser.add_argument("--mask",
                        action='store_true',
                        dest='mask',
                        help="Mask Higgs bins")
    parser.add_argument("--covonly",
                        action='store_true',
                        help="only cov plots")
    parser.add_argument("--all",
                        action='store_true',
                        dest='run_all',
                        help="Include split pT bin plots")
    parser.add_argument("--run2", action='store_true', dest='run2', help="Stack all years")
    parser.add_argument("--paper", action='store_true', dest='paper', help="Stack all years")
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

    parser.add_argument("--scale",
                        type=int,
                        default=200,
                        help="Scale Higgs signal in plots by X")

    parser.add_argument("--filled",
                        type=str2bool,
                        default='True',
                        choices={True, False},
                        help="Use filled stack plots")
    parser.add_argument("--resample",
                        action='store_true',
                        help="Resample to get yerr")
    parser.add_argument("--resampled",
                        action='store_true',
                        help="Use existing")

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

    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    parser.add_argument("--debug", "-vv", action="store_true", help="Debug logging")
    args = parser.parse_args()

    log_level = logging.WARNING
    if args.verbose:
        log_level = logging.INFO
    if args.debug:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)

    if args.output_folder.split("/")[0] != args.dir:
        args.output_folder = os.path.join(args.dir, args.output_folder)
    if args.resampled:
        args.resampled = "resamp_{fit}.pkl.gz"

    if not args.run2:
        if args.year is None:
            configs = json.load(open("config.json"))
            args.year = str(configs['year'])

    if args.year == "2017":
        pbins = [475, 500, 550, 600, 675, 800, 1200]
    else:
        pbins = [450, 500, 550, 600, 675, 800, 1200]


    from functools import partial
    full_plot = partial(full_plot, pbins=pbins, args=args, fitDiag=r.TFile.Open(args.input))

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
    if not args.covonly:
        f = uproot.open(args.input)
        for shape_type in shape_types:
            if args.resample:
                # python ../correlations/resampler.py -w model_combined.root:w -f fitDiagnosticsTest.root:fit_s --out here.pkl.gz -o msd
                cmd_str = "python ../resampler.py -w model_combined.root:w -f {inp}:{fit} --out resamp_{fit}.pkl.gz -o msd".format(inp=args.input, fit=shape_type)
                print(cmd_str)
                os.system(cmd_str)
                args.resampled = "resamp_{fit}.pkl.gz".format(fit=shape_type)
            if args.resampled:
                import pickle, gzip
                print("Reading from " +  args.resampled.format(fit=shape_type))
                try:
                    with gzip.open(args.resampled.format(fit=shape_type), 'rb') as fin:
                        nominal, varied, ntoys = pickle.load(fin)
                except:
                    raise ValueError("No resampled pkl found.")
            else:
                nominal, varied, ntoys = None, None, None

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
                    if args.run2:
                        cat_list = ['shapes_{}/ptbin{}{}{};1'.format(shape_type, i, region, year) for year in ['2016', '2017', '2018'] ]
                    else:
                        cat_name = 'shapes_{}/ptbin{}{}{};1'.format(shape_type, i, region,
                                                                    args.year)
                        try:
                            cat_list = [cat_name]
                        except Exception:
                            raise ValueError("Namespace {} is not available, only following"
                                            "namespaces were found in the file: {}".format(
                                                args.fit, f.keys()))
                    template_set = [f[cat] for cat in cat_list]
                    if nominal is not None:
                        all_temp_names = sum([get_cat_keys(cat.split("/")[-1].split(";")[0], nominal) for cat in cat_list], [])
                        _bkgunc = get_yerr(nominal, varied, all_temp_names, ntoys)
                    else:
                        _bkgunc = None

                    fig = full_plot(template_set,
                                    pseudo=args.pseudo,
                                    fittype=shape_type,
                                    mask=mask,
                                    toys=args.toys,
                                    # sqrtnerr=True,
                                    bkgUnc=_bkgunc,
                                    run2=args.run2,
                                    paper=args.paper,
                                    scaleH=args.scale,
                                    format=args.format,
                                    year="" if args.run2 else args.year,
                                    )

                if args.run2:
                    cat_list = [
                        'shapes_{}/ptbin{}{}{};1'.format(shape_type, i, region, '2016')
                        for i in range(0, 6)
                    ]
                    cat_list += [
                        'shapes_{}/ptbin{}{}{};1'.format(shape_type, i, region, '2017')
                        for i in range(0, 6)
                    ]
                    cat_list += [
                        'shapes_{}/ptbin{}{}{};1'.format(shape_type, i, region, '2018')
                        for i in range(0, 6)
                    ]
                else:
                    cat_list = [
                        'shapes_{}/ptbin{}{}{};1'.format(shape_type, i, region,
                                                         args.year)
                        for i in range(0, 6)
                    ]
                template_set = [f[cat] for cat in cat_list]
                if nominal is not None:
                    all_temp_names = sum([get_cat_keys(cat.split("/")[-1].split(";")[0], nominal) for cat in cat_list], [])
                    _bkgunc = get_yerr(nominal, varied, all_temp_names, ntoys)
                else:
                    _bkgunc = None

                full_plot(template_set,
                        pseudo=args.pseudo,
                        fittype=shape_type,
                        mask=mask,
                        toys=args.toys,
                        # sqrtnerr=True,
                        bkgUnc=_bkgunc,
                        format=args.format,
                        year="" if args.run2 else args.year,
                        scaleH=args.scale,
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
                                # sqrtnerr=True,
                                style_set=style_set_B,
                                stack_style=1, # stack al always
                                format=args.format,
                                year="",
                                scaleH=args.scale,
                                run2=args.run2,
                                )
                    else:
                        cat = f['shapes_{}/muonCR{}{};1'.format(shape_type, region, args.year)]
                        full_plot(
                            [cat],
                            args.pseudo,
                            fittype=shape_type,
                            toys=args.toys,
                            # sqrtnerr=True,
                            style_set=style_set_B,
                            format=args.format,
                            year=args.year,
                            scaleH=args.scale,
                        )
                    print("Plotted muCR", region, shape_type)
                else:
                    print("Muon region not found")

    ##### Input shape plotter
    # Take sqrt N err for data
    # Mock QCD while unavailable as template in rhalpha
    import os
    from rhalphalib.plot.input_shapes import input_dict_maker

    if not args.run2:
        mockd = input_dict_maker(os.getcwd() + ".pkl")
        input_pseudo = True
        if args.toys or not args.pseudo:
            input_pseudo = False
        for shape_type in ["inputs"]:
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
                            # sqrtnerr=True,
                            toys=False,
                            format=args.format,
                            year=args.year,
                            scaleH=args.scale,
                            )
                # Per bin plots
                for i in range(0, 6):
                    if not args.run_all: continue
                    full_plot(
                        [mockd['ptbin{}{}{}_{}'.format(i, region, args.year, shape_type)]],
                        pseudo=input_pseudo,
                        fittype=shape_type,
                        mask=mask,
                        # sqrtnerr=True,
                        toys=False,
                        format=args.format,
                        year=args.year,
                        scaleH=args.scale,
                        )

                # MuonCR if included
                try:
                    cat = mockd['muonCR{}{}_{}'.format(region, args.year, shape_type)]
                    full_plot([cat],
                                fittype=shape_type,
                                pseudo=input_pseudo,
                                mask=False,
                                # sqrtnerr=True,
                                toys=False,
                                format=args.format,
                                year=args.year,
                                scaleH=args.scale
                    )
                    print("Plotted input, muCR", region, shape_type)
                except Exception:
                    print("Muon region not found")
                    pass
        # except:
        #     print("Input pkl file not found")

    if args.three_regions:
        plot_fractions(
            os.path.join(args.dir, args.input),
            os.path.join(args.dir, 'model_combined.root'),
            out='{}/{}.{}'.format(args.output_folder, 'fractions', args.format),
            data=((not args.pseudo) | args.toys),
            year=args.year,
        )

    plot_cov(
        os.path.join(args.dir, args.input),
        out='{}/{}.{}'.format(args.output_folder, 'covariances', args.format),
        data=((not args.pseudo) | args.toys),
        year=args.year,
    )

    plot_cov(
        os.path.join(args.dir, args.input),
        out='{}/{}_wTF.{}'.format(args.output_folder, 'covariances', args.format),
        data=((not args.pseudo) | args.toys),
        year=args.year,
        include='tf',
    )

    plot_cov(
        os.path.join(args.dir, args.input),
        out='{}/{}_filtered.{}'.format(args.output_folder, 'covariances', args.format),
        data=((not args.pseudo) | args.toys),
        year=args.year,
        include='r,z,r16,r17,r18,z16,z17,z18,*eff*,*scale_2*,*smear_2*',
    )
