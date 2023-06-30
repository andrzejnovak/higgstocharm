#!/usr/bin/env python
import argparse
import numpy as np
import tqdm
import pickle
import gzip

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gEnv.SetValue("RooFit.Banner=0")
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)


def _RooAbsCollection__iter__(self):
    it = self.iterator()
    obj = it.Next()
    while obj != None:  # noqa: E711
        yield obj
        obj = it.Next()


def _RooAbsCollection_assign(self, other):
    if self == other:
        return
    for el in self:
        if not hasattr(el, 'setVal'):
            continue
        theirs = other.find(el)
        if not theirs:
            continue
        el.setVal(theirs.getVal())
        el.setError(theirs.getError())
        el.setAsymError(theirs.getErrorLo(), theirs.getErrorHi())
        el.setAttribute("Constant", theirs.isConstant())


def _RooArgList_fromiter(cls, iterable, silent=False):
    items = cls()
    for item in iterable:
        items.add(item, silent)
    return items


ROOT.RooAbsCollection.assign = _RooAbsCollection_assign
ROOT.RooAbsCollection.__iter__ = _RooAbsCollection__iter__
ROOT.RooArgList.fromiter = classmethod(_RooArgList_fromiter)


def to_numpy(hinput):
    if "<class 'ROOT.TH1" in str(type(hinput)):
        sumw = np.zeros(hinput.GetNbinsX())
        binning = np.zeros(sumw.size + 1)
        name = hinput.GetName()
        for i in range(1, sumw.size + 1):
            sumw[i-1] = hinput.GetBinContent(i)
            binning[i-1] = hinput.GetXaxis().GetBinLowEdge(i)
        binning[i] = hinput.GetXaxis().GetBinUpEdge(i)
        return (sumw, binning, name)
    else:
        raise TypeError("Don't know how to convert %r to numpy" % type(hinput))


def resample_fit_result(args):
    fws, wsname = args.workspace.split(':')
    fin = ROOT.TFile.Open(fws)
    w = fin.Get(wsname)

    ffit, fitname = args.fit.split(':')
    fin2 = ROOT.TFile.Open(ffit)
    if fitname != "prefit":
        fit = fin2.Get(fitname)
    else:
        fit = fin2.Get("fit_s")

    model = w.obj(args.model)
    obs = model.GetObservables()[args.observable]
    pdf = model.GetPdf()
    params = pdf.getParameters(model.GetObservables())
    # It would seem necessary to do this but sometimes the fit
    # constant parameters are somehow the wrong value w.r.t. the model
    # params.assign(fit.constPars())
    if fitname != "prefit":
        params.assign(fit.floatParsFinal())
    else:
        params.assign(fit.floatParsInit())

    all_processes = []
    cat = pdf.indexCat()
    for iCat in range(cat.numBins('')):
        cat.setBin(iCat)
        catpdf = pdf.getPdf(cat.getLabel())
        addpdf = [p for p in catpdf.pdfList() if p.dependsOn(model.GetObservables())]
        if len(addpdf) != 1:
            raise RuntimeError("I failed to parse the model structure :(")
        addpdf = addpdf[0]
        if not addpdf.dependsOn(obs):
            continue
        for proc, norm in zip(addpdf.pdfList(), addpdf.coefList()):
            all_processes.append((proc, norm))

    nominal = {}
    for proc, norm in all_processes:
        h = proc.createHistogram("", obs)
        h.SetDirectory(0)
        h.Scale(norm.getVal())
        nominal[proc.GetName()] = to_numpy(h)
        del h

    varied = {name: np.zeros((args.samples, ) + h[0].shape) for name, h in nominal.items()}
    for i in tqdm.trange(args.samples):
        params.assignValueOnly(fit.randomizePars())
        for proc, norm in all_processes:
            h = proc.createHistogram("", obs)
            h.Scale(norm.getVal())
            h.SetDirectory(0)
            varied[proc.GetName()][i] = to_numpy(h)[0]
            del h

    return nominal, varied, args.samples


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Resample fit result and save as a pickled dictionary of numpy arrays for subsequent processing.")
    parser.add_argument("-w", "--workspace", metavar="ROOTFILE:WORKSPACE", help="Workspace to load (e.g. card.root:w)", required=True)
    parser.add_argument("-f", "--fit", metavar="ROOTFILE:FIT_NAME", help="Fit result to load (e.g. fitDiagnostics.root:fit_s)", required=True)
    parser.add_argument("-m", "--model", help="Model to load (default: %(default)s)", default="ModelConfig")
    parser.add_argument("-o", "--observable", help="Observable (default: %(default)s)", default="x")
    parser.add_argument("-n", "--samples", help="Number of samples (default: %(default)s)", type=int, default=200)
    parser.add_argument("--out", help="Output filename.  Suggested suffix: .pkl.gz", type=str, required=True)

    args = parser.parse_args()
    arrays = resample_fit_result(args)
    with gzip.open(args.out, 'wb') as fout:
        pickle.dump(arrays, fout)
