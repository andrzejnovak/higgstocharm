from __future__ import print_function, division
import warnings
import rhalphalib as rl
import pickle
import numpy as np
np.set_printoptions(linewidth=1000, precision=2)
import ROOT
import uproot
from rhalphalib import AffineMorphTemplate, MorphHistW2
from util import make_dirs

rl.util.install_roofit_helpers()

SF = {
    "2016": {
        'V_SF': 0.875,
        'V_SF_ERR': 0.017,
        'W_SF': 0.605,
        'W_SF_ERR': 0.088,
        'shift_SF': -1.245,
        'shift_SF_ERR': 0.21,
        'smear_SF': 1.02415,
        'smear_SF_ERR': 0.0174,
    },
    "2017": {
        # Default
        'V_SF': 0.881,
        'V_SF_ERR': 0.022,
        'W_SF': 0.663,
        'W_SF_ERR': 0.091,
        'shift_SF': -1.425,
        'shift_SF_ERR': 0.21,
        'smear_SF': 1.01605,
        'smear_SF_ERR': 0.0131
    },
    "2018": {
        'V_SF': 0.896,
        'V_SF_ERR': 0.013,
        'W_SF': 0.714,
        'W_SF_ERR': 0.079,
        'shift_SF': -0.304,
        'shift_SF_ERR': 0.15,
        'smear_SF': 1.05,
        'smear_SF_ERR': 0.00235,
    },
}


def badtemp(hvalues, eps=0.0000001, mask=None):
    # Need minimum size & more than 1 non-zero bins
    tot = np.sum(hvalues[mask])
    count_nonzeros = np.sum(hvalues[mask] > 0)
    if (tot < eps) or (count_nonzeros < 2):
        return True
    else:
        return False


def smass(sName):
    if 'hbb' in sName or 'hcc' in sName:
        _mass = 125.
    elif sName in ['wqq', 'wcq', 'tqq', 'stqq', 'vvqq']:
        _mass = 80.379
    elif sName in ['zqq', 'zcc', 'zbb']:
        _mass = 91.
    else:
        raise ValueError("DAFUQ is {}".format(sName))
    return _mass


def passfailSF(f, region, sName, ptbin, mask, SF=1, SF_unc=0.1, muon=False):
    """
    Return (SF, SF_unc) for a pass/fail scale factor.
    """
    if region == "pass":
        return SF, 1. + SF_unc / SF
    else:
        _pass = get_templ(f, "pass", sName, ptbin, muon=muon)
        _pass_rate = np.sum(_pass[0] * mask)

        _fail = get_templ(f, "fail", sName, ptbin, muon=muon)
        _fail_rate = np.sum(_fail[0] * mask)

        if _fail_rate > 0:
            _sf = 1 + (1 - SF) * _pass_rate / _fail_rate
            _sfunc = 1. - SF_unc * (_pass_rate / _fail_rate)
            return _sf, _sfunc
        else:
            return 1, 1


def shape_to_num(f, region, sName, ptbin, syst, mask, muon=False, bound=0.5):
    _nom = get_templ(f, region, sName, ptbin, muon=muon)
    if _nom is None:
        return None
    _nom_rate = np.sum(_nom[0] * mask)
    if _nom_rate < .1:
        return 1.0
    _up = get_templ(f, region, sName, ptbin, syst=syst + "Up", muon=muon)
    _down = get_templ(f, region, sName, ptbin, syst=syst + "Down", muon=muon)
    if _up is None and _down is None:
        return None
    else:
        if _down is not None:
            _down_rate = np.sum(_down[0] * mask)
        else:
            _down_rate = _nom_rate
        if _up is not None:
            _up_rate = np.sum(_up[0] * mask)
        else:
            _up_rate = _nom_rate
    _diff = np.abs(_up_rate - _nom_rate) + np.abs(_down_rate - _nom_rate)
    magnitude = _diff / (2. * _nom_rate)
    if bound is not None:
        magnitude = min(magnitude, bound)
    return 1.0 + magnitude

def get_templ(f, region, sample, ptbin, syst=None, muon=False):
    hist_name = '{}_{}'.format(sample, region)
    if syst is not None:
        hist_name += "_" + syst
    else:
        hist_name += "_nominal"
    if not muon:
        hist_name += "_bin{}".format(ptbin)
    try:
        f[hist_name]
    except:
        print("{}Sample {}, {}, {}, {} not found.".format('(Muon) ' if muon else "",
            sample, region, ptbin if not muon else "-", syst))
        return None
    h_vals = f[hist_name].values
    h_edges = f[hist_name].edges
    h_variances = f[hist_name].variances
    if np.any(h_vals < 0):
        print("Sample {}, {}, {}, {}, has {} negative bins. They will be set to 0.".format(
            sample, region, ptbin, syst, np.sum(h_vals < 0)))
        _invalid = h_vals < 0
        h_vals[_invalid] = 0
        h_variances[_invalid] = 0
    if np.any(~np.isfinite(h_vals)):
        print("Sample {}, {}, {}, {}, has {} Nan/Inf bins. They will be set to 0.".format(
            sample, region, ptbin, syst, np.sum(~np.isfinite(h_vals))))
        _invalid = ~np.isfinite(h_vals)
        h_vals[_invalid] = 0
        h_variances[_invalid] = 0
    h_key = 'msd'
    return (h_vals, h_edges, h_key, h_variances)


def one_bin(template):
    N = 5
    try:
        h_vals, h_edges, h_key, h_variances = template
        return (np.array([np.sum(h_vals)]), np.array([0., 1.]), "onebin", np.array([np.sum(h_variances)]))
    except:
        h_vals, h_edges, h_key = template
        return (np.array([np.sum(h_vals)]), np.array([0., 1.]), "onebin")

def dummy_rhalphabet(pseudo,
                     throwPoisson,
                     MCTF,
                     justZ=False,
                     scale_syst=True,
                     smear_syst=True,
                     systs=True,
                     blind=True,
                     runhiggs=False,
                     fitTF=True,
                     muonCR=True,
                     runboth=False,
                     year=2017,
                     opts=None):

    # Default lumi (needs at least one systematics for prefit)
    sys_lumi = rl.NuisanceParameter('CMS_lumi_13TeV_{}'.format(year), 'lnN')
    sys_lumi_correlated = rl.NuisanceParameter('CMS_lumi_13TeV_correlated', 'lnN')
    sys_lumi_1718 = rl.NuisanceParameter('CMS_lumi_13TeV_1718', 'lnN')
    lumi_dict = {
        "2016": 1.01,
        "2017": 1.02,
        "2018": 1.015,
    }
    lumi_correlated_dict = {
        "2016": 1.006,
        "2017": 1.009,
        "2018": 1.02,
    }
    lumi_1718_dict = {
        "2017": 1.006,
        "2018": 1.002,
    }
    # TT params
    tqqeffSF = rl.IndependentParameter('tqqeffSF_{}'.format(year), 1., 0, 10)
    tqqnormSF = rl.IndependentParameter('tqqnormSF_{}'.format(year), 1., 0, 10)
    # Systematics
    sys_shape_dict = {}
    if opts.fast == 0:
        _shape_type = 'lnN'
    else:
        _shape_type = 'shape'
    sys_shape_dict['JES'] = rl.NuisanceParameter('CMS_scale_j_{}'.format(year), _shape_type)
    sys_shape_dict['JER'] = rl.NuisanceParameter('CMS_res_j_{}'.format(year), _shape_type)
    sys_shape_dict['UES'] = rl.NuisanceParameter('CMS_ues_j_{}'.format(year), _shape_type)
    sys_shape_dict['jet_trigger'] = rl.NuisanceParameter('CMS_gghcc_trigger_{}'.format(year), _shape_type)
    sys_shape_dict['mu_trigger'] = rl.NuisanceParameter('CMS_mu_trigger_{}'.format(year), _shape_type)
    sys_shape_dict['mu_isoweight'] = rl.NuisanceParameter('CMS_mu_isoweight_{}'.format(year), _shape_type)
    sys_shape_dict['mu_idweight'] = rl.NuisanceParameter('CMS_mu_idweight_{}'.format(year), _shape_type)
    sys_shape_dict['pileup_weight'] = rl.NuisanceParameter('CMS_gghcc_PU_{}'.format(year), _shape_type)
    sys_shape_dict['HEM18'] = rl.NuisanceParameter('CMS_gghcc_HEM_{}'.format(year), _shape_type)
    sys_shape_dict['L1Prefiring'] = rl.NuisanceParameter('CMS_gghcc_L1prefire_{}'.format(year), _shape_type)
    sys_shape_dict['scalevar_7pt'] = rl.NuisanceParameter('CMS_gghcc_th_scale7pt', _shape_type)
    sys_shape_dict['scalevar_3pt'] = rl.NuisanceParameter('CMS_gghcc_th_scale3pt', _shape_type)

    for sys in ['btagEffStat', 'btagWeight', 'd1kappa_EW', 'Z_d2kappa_EW', 'Z_d3kappa_EW', 'd1K_NLO', 'd2K_NLO', 'd3K_NLO']:
        sys_shape_dict[sys] = rl.NuisanceParameter('CMS_gghcc_{}_{}'.format(sys, year), _shape_type)

    sys_ddxeff = rl.NuisanceParameter('CMS_eff_cc_{}'.format(year), 'lnN')
    sys_ddxeffbb = rl.NuisanceParameter('CMS_eff_bb_{}'.format(year), 'lnN')
    sys_ddxeffw = rl.NuisanceParameter('CMS_eff_w_{}'.format(year), 'lnN')

    sys_eleveto = rl.NuisanceParameter('CMS_gghcc_e_veto_{}'.format(year), 'lnN')
    sys_muveto = rl.NuisanceParameter('CMS_gghcc_m_veto_{}'.format(year), 'lnN')
    sys_tauveto = rl.NuisanceParameter('CMS_gghcc_tau_veto_{}'.format(year), 'lnN')

    sys_wznormEW = rl.NuisanceParameter('CMS_gghcc_wznormEW', 'lnN')
    sys_znormEW = rl.NuisanceParameter('CMS_gghcc_znormEW', 'lnN')
    sys_znormQ = rl.NuisanceParameter('CMS_gghcc_znormQ', 'lnN')

    sys_veff = rl.NuisanceParameter('CMS_gghcc_veff_{}'.format(year), 'lnN')
    sys_scale = rl.NuisanceParameter('CMS_gghcc_scale_{}'.format(year), 'shape')
    sys_smear = rl.NuisanceParameter('CMS_gghcc_smear_{}'.format(year), 'shape')

    # sys_ggHpt = rl.NuisanceParameter('CMS_gghcc_th_scale7pt', 'lnN')
    # sys_VBFpt = rl.NuisanceParameter('CMS_gghcc_th_scale3pt', 'lnN')

    # Import binnings
    # Hidden away to be available to other functions
    from config_Hxx import ptbins, msdbins  # ptpts, msdpts, rhopts
    from config_Hxx import ptscaled, rhoscaled, validbins
    from config_Hxx import ddxSF
    msd = rl.Observable('msd', msdbins)
    npt = len(ptbins) - 1

    # Year setup
    if opts.templates is not None:
        print("Custom templates:", opts.templates)
        model_name = "Nano17Model"
        #f = uproot.open('hxx/hist_1DZcc_pt_scalesmear.root')
        f = uproot.open(opts.templates)
        if muonCR:
            f_mu = uproot.open(opts.mutemplates)
        #muonCR = False
    # elif year == "2018":
    #     print("Year: 2018")
    #     model_name = "temp18Model"
    #     #f = uproot.open('hxx18/hist_1DZcc_pt_scalesmear.root')
    #     f = uproot.open('2018v2/hist_1DZcc_pt_scalesmear.root')
    #     f_mu = uproot.open('2018v2/hist_1DZcc_muonCR.root')
    # elif year == "2017":
    #     print("Year: 2017")
    #     model_name = "temp17Model"
    #     #f = uproot.open('hxx/hist_1DZcc_pt_scalesmear.root')
    #     f = uproot.open('2017v2/hist_1DZcc_pt_scalesmear.root')
    #     f_mu = uproot.open('2017v2/hist_1DZcc_muonCR.root')
    # elif year == "2016":
    #     print("Year: 2016")
    #     model_name = "temp16Model"
    #     f = uproot.open('2016v2/hist_1DZcc_pt_scalesmear.root')
    #     f_mu = uproot.open('2016v2/hist_1DZcc_muonCR.root')
    # else:
    #     raise ValueError("Invalid Year")
    if opts.model is not None:
        model_name = opts.model

    # Prepare dir structure
    make_dirs('{}/plots'.format(model_name))

    # Get QCD efficiency
    if MCTF:
        qcdmodel = rl.Model("qcdmodel")

    qcdpass, qcdfail = 0., 0.
    for ptbin in range(npt):
        failCh = rl.Channel("ptbin%d%s" % (ptbin, 'fail'))
        passCh = rl.Channel("ptbin%d%s" % (ptbin, 'pass'))

        passTempl = get_templ(f, "pass", "qcd", ptbin)
        failTempl = get_templ(f, "fail", "qcd", ptbin)

        failCh.setObservation(failTempl, read_sumw2=True)
        passCh.setObservation(passTempl, read_sumw2=True)
        qcdfail += failCh.getObservation()[0].sum()
        qcdpass += passCh.getObservation()[0].sum()

        if MCTF:
            qcdmodel.addChannel(failCh)
            qcdmodel.addChannel(passCh)

    qcdeff = qcdpass / qcdfail

    # Separate out QCD to QCD fit
    if MCTF:
        degsMC = tuple([int(s) for s in opts.degsMC.split(',')])
        _basisMC = opts.basis.split(",")[0]
        if _basisMC == 'Bernstein':
            _inits = np.ones(tuple(n + 1 for n in degsMC))
        elif _basisMC == 'Chebyshev':
            _inits = np.zeros(tuple(n + 1 for n in degsMC))
            _inits[0,0] = 1
        else:
            raise ValueError("Basis ``{}`` not understood.".format(_basisMC))
        tf_MCtempl = rl.BasisPoly("tf{}_MCtempl".format(year),
                                      degsMC, ['pt', 'rho'], basis=_basisMC,
                                      init_params = _inits,
                                      limits=(-10, 10), coefficient_transform=None)
        tf_MCtempl_params = qcdeff * tf_MCtempl(ptscaled, rhoscaled)

        for ptbin in range(npt):
            failCh = qcdmodel['ptbin%dfail' % ptbin]
            passCh = qcdmodel['ptbin%dpass' % ptbin]
            failObs = failCh.getObservation()[0]
            qcdparams = np.array([
                rl.IndependentParameter('qcdparam_ptbin%d_msdbin%d' % (ptbin, i), 0)
                for i in range(msd.nbins)
            ])
            sigmascale = 10.
            scaledparams = failObs * (
                1 + sigmascale / np.maximum(1., np.sqrt(failObs)))**qcdparams
            fail_qcd = rl.ParametericSample('ptbin%dfail_qcd' % ptbin,
                                            rl.Sample.BACKGROUND, msd, scaledparams)
            failCh.addSample(fail_qcd)
            pass_qcd = rl.TransferFactorSample('ptbin%dpass_qcd' % ptbin,
                                               rl.Sample.BACKGROUND,
                                               tf_MCtempl_params[ptbin, :], fail_qcd)
            passCh.addSample(pass_qcd)

            failCh.mask = validbins[ptbin]
            passCh.mask = validbins[ptbin]

        qcdfit_ws = ROOT.RooWorkspace('qcdfit_ws')
        simpdf, obs = qcdmodel.renderRoofit(qcdfit_ws)
        qcdfit = simpdf.fitTo(
            obs,
            ROOT.RooFit.Extended(True),
            ROOT.RooFit.SumW2Error(True),
            ROOT.RooFit.Strategy(2),
            ROOT.RooFit.Save(),
            ROOT.RooFit.Minimizer('Minuit2', 'migrad'),
            ROOT.RooFit.Offset(True),
            ROOT.RooFit.PrintLevel(-1),
        )
        qcdfit_ws.add(qcdfit)
        qcdfit_ws.writeToFile('{}/qcdfit.root'.format(model_name))
        if qcdfit.status() != 0:
            qcdfit.Print()
            raise RuntimeError('Could not fit qcd')

        qcdmodel.readRooFitResult(qcdfit)

        # # Plot it
        from rhalphalib.plot.plot_TF import TF_smooth_plot, TF_params
        from rhalphalib.plot.plot_TF import plotTF as plotMCTF
        _values = [par.value for par in tf_MCtempl.parameters.flatten()]
        _names = [par.name for par in tf_MCtempl.parameters.flatten()]
        np.save('{}/MCTF'.format(model_name), _values)
        print('ptdeg', degsMC[0], 'rhodeg', degsMC[1])
        plotMCTF(*TF_smooth_plot(*TF_params(_values, _names)),
                 MC=True,
                 raw=True,
                 ptdeg=degsMC[0],
                 rhodeg=degsMC[1],
                 year=args.year,
                 out='{}/plots/TF_MC_only'.format(model_name))

        param_names = [p.name for p in tf_MCtempl.parameters.reshape(-1)]
        decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(
            tf_MCtempl.name + '_deco', qcdfit, param_names)
        np.save('{}/decoVector'.format(model_name), decoVector._transform)
        tf_MCtempl.parameters = decoVector.correlated_params.reshape(
            tf_MCtempl.parameters.shape)
        tf_MCtempl_params_final = tf_MCtempl(ptscaled, rhoscaled)

    # build actual fit model now
    model = rl.Model("shapes" + year)
    if not justZ:
        model.t2w_config = ("-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose "
                            "--PO 'map=.*/*hcc*:r[1,-500,500]' --PO 'map=.*/zcc:z[1,-5,5]'")

    for ptbin in range(npt):
        for region in ['pass', 'fail']:
            #ch = rl.Channel("ptbin%d%s" % (ptbin, region))
            ch = rl.Channel("ptbin{}{}{}".format(ptbin, region, year))
            model.addChannel(ch)
            if justZ:
                include_samples = ['zcc']
            elif opts.justHZ is True:
                include_samples = ['zcc', "hcc"]
            else:
                include_samples = [
                    'zbb', 'zcc', 'zqq',
                    'wcq', 'wqq',
                    'tqq', 'stqq', 'vvqq', 'zll', 'wln',
                    'hcc', 'vbfhcc', 'whcc', 'zhcc',
                    'hbb', 'vbfhbb', 'whbb', 'zhbb', # hbb signals
                ]
            # Remove unavailable samples
            _available = sorted(
                list(set([key.split("_pass")[0] for key in f.keys() if "pass" in key])))
            found_samples = []
            for sName in include_samples:
                if sName not in _available:
                    print(
                        'Sample `{}` is not available in templates file.'.format(sName))
                else:
                    found_samples.append(sName)
            include_samples = found_samples
            #print("x", found_samples)

            # Define mask
            mask = validbins[ptbin].copy()
            if not pseudo and region == 'pass':
                if blind:
                    #mask[4:7] = False
                    mask[6:9] = False
                    mask[10:14] = False

            from functools import partial
            # badtemp_ma = partial(badtemp, mask=mask)
            
            def badtemp_ma(hvalues, eps=0.0000001, mask=None):
                # Need minimum size & more than 1 non-zero bins
                tot = np.sum(hvalues[mask])
                count_nonzeros = np.sum(hvalues[mask] > 0)
                if (tot < eps) or (count_nonzeros < 3):
                    return True
                else:
                    return False
            # Remove empty samples
            for sName in include_samples:
                templ = get_templ(f, region, sName, ptbin)
                if templ is None:
                    print(
                        'Sample {} in region = {}, ptbin = {}, not found in template file.'
                        .format(sName, region, ptbin))
                    include_samples.remove(sName)
                elif np.sum(templ[0][mask]) < 0.0000001:
                    print(
                        'Sample {} in region = {}, ptbin = {}, would be empty, so it will be removed'
                        .format(sName, region, ptbin))
                    include_samples.remove(sName)

            if not fitTF:  # Add QCD sample when not running TF fit
                include_samples.append('qcd')
            for sName in include_samples:
                templ = get_templ(f, region, sName, ptbin)
                if opts.runbb:
                    _signals = ["zbb"]
                    raise NotImplementedError
                else:
                    if sName == "zcc":
                        stype = 0
                    elif "hcc" in sName:
                        stype = -1
                    else:
                        stype = rl.Sample.BACKGROUND

                MORPHNOMINAL = True
                def smorph(templ):      
                    if templ is None:
                        return None                  
                    if MORPHNOMINAL and sName not in ['qcd', 'zll', 'wln', 'vvqq']:
                        return MorphHistW2(templ).get(shift=SF[year]['shift_SF']/smass('wcq') * smass(sName),
                                                      smear=SF[year]['smear_SF']
                                                      )
                    else:
                        return templ
                templ = smorph(templ)
                    
                sample = rl.TemplateSample(ch.name + '_' + sName, stype, templ)

                # Systematics
                #####################################################
                if not systs:  # Need at least one
                    sample.setParamEffect(sys_lumi, lumi_dict[year])
                else:
                    sample.setParamEffect(sys_lumi, lumi_dict[year])
                    sample.setParamEffect(sys_lumi_correlated, lumi_correlated_dict[year])
                    if year != '2016':
                        sample.setParamEffect(sys_lumi_1718, lumi_1718_dict[year])
                    sample.setParamEffect(sys_eleveto, 1.005)
                    sample.setParamEffect(sys_muveto, 1.005)
                    sample.setParamEffect(sys_tauveto, 1.005)

                    if sName in ["qcd"]:
                        continue

                    sys_names = [
                        'JES', 'JER', 'UES', 'jet_trigger', 'btagEffStat', 'btagWeight', 'pileup_weight',
                        'Z_d2kappa_EW', 'Z_d3kappa_EW', 'd1kappa_EW', 'd1K_NLO', 'd2K_NLO', 'd3K_NLO',
                        'L1Prefiring', 
                        'scalevar_7pt', 'scalevar_3pt',
                    ]
                    for sys_name in sys_names:
                        if (("NLO" in sys_name) or ("EW" in sys_name)) and not sName in ['zbb', 'zcc', 'zqq', 'wcq', 'wqq']:
                            continue
                        if ("Z_d" in sys_name) and sName not in ['zbb', 'zcc', 'zqq']:
                            continue
                        if sys_name == 'scalevar_7pt' and sName not in ['hcc', 'hbb']:
                            continue
                        if sys_name == 'scalevar_3pt' and sName not in ['vbfhcc', 'vbfhbb', 'whcc', 'zhcc', 'whbb', 'zhbb']:
                            continue
                        if opts.fast == 0:  # Convert to lnN for faster fitting
                            _sys_ef = shape_to_num(f, region, sName, ptbin,
                                                    sys_name, mask, bound=None if 'scalevar' not in sys_name else 0.25)
                            if _sys_ef is None:
                                continue
                            sample.setParamEffect(sys_shape_dict[sys_name], _sys_ef)
                        else:
                            _up = get_templ(f,
                                            region,
                                            sName,
                                            ptbin,
                                            syst=sys_name + "Up")
                            _dn = get_templ(f,
                                            region,
                                            sName,
                                            ptbin,
                                            syst=sys_name + "Down")
                            _up = smorph(_up)
                            _dn = smorph(_dn)
                            if _up is None or _dn is None:
                                continue
                            if badtemp_ma(_up[0]) or badtemp_ma(_dn[0]):
                                print("Skipping {} for sample {}, shape would be empty".format(sys_name, sName))
                                continue
                            sample.setParamEffect(sys_shape_dict[sys_name], _up[0], _dn[0])
                    # One sided
                    for sys_name in ['HEM18']:
                        if opts.fast == 0:  # Convert to lnN for faster fitting
                            _sys_ef = shape_to_num(f, region, sName, ptbin,
                                                    sys_name, mask)
                            if _sys_ef is None:
                                continue
                            sample.setParamEffect(sys_shape_dict[sys_name], _sys_ef)
                        else:                            
                            _up = get_templ(f,
                                            region,
                                            sName,
                                            ptbin,
                                            syst=sys_name + "Up")
                            _up = smorph(_up)
                            if _up is None:
                                continue
                            if badtemp_ma(_up[0]):
                                print("Skipping {} for sample {}, shape would be empty".format(sys_name, sName))
                                continue
                            sample.setParamEffect(sys_shape_dict[sys_name], _up[0])

                    if opts.mcstat and sName not in ['qcd']:
                        if opts.fast < 2:  # Convert to lnN for faster fitting
                            sample.autoMCStats(lnN=True)
                        else:
                            sample.autoMCStats(epsilon=1e-4)

                    # if sName not in ['tqq', 'stqq']:
                    if sName not in ["qcd", 'zll', 'wln']:
                        sample.scale(SF[year]['V_SF'])
                        sample.setParamEffect(
                            sys_veff, 1.0 + SF[year]['V_SF_ERR'] / SF[year]['V_SF'])
                    if sName in ["zcc", "hcc"]:
                        # _sf, _sfunc = passfailSF(f, region, sName, ptbin, mask,
                        #                          SF[year]['CC_SF'],
                        #                          SF[year]['CC_SF_ERR'])
                        _sf, _sfup, _sfdn = ddxSF(ptbin, 'cc')
                        _sf, _sfunc = passfailSF(f, region, sName, ptbin, mask,
                                                 _sf,
                                                 _sfdn)
                        sample.scale(_sf)
                        sample.setParamEffect(sys_ddxeff, _sfunc)
                    if 'bb' in sName:
                        _sf, _sfunc = passfailSF(f, region, sName, ptbin, mask, 1, 0.3)
                        sample.scale(_sf)
                        sample.setParamEffect(sys_ddxeffbb, _sfunc)
                    if sName in ["wcq", "wqq"]:
                        _sf, _sfunc = passfailSF(f, region, sName, ptbin, mask,
                                                 SF[year]['W_SF'], SF[year]['W_SF_ERR'])
                        sample.scale(_sf)
                        sample.setParamEffect(sys_ddxeffw, _sfunc)
                    # if sName == 'hcc':
                    #     sample.setParamEffect(sys_ggHpt, 1.2)
                    # if sName == 'hcc':
                    #     sample.setParamEffect(sys_VBFpt, 1.2)

                # Scale and Smear
                mtempl = AffineMorphTemplate(templ)

                if scale_syst and sName not in ['qcd', 'zll', 'wln', 'vvqq']:
                    # realshift = smass(sName) * SF[year]['shift_SF'] * SF[year]['shift_SF_ERR']
                    realshift = SF[year]['shift_SF_ERR']/smass('wcq') * smass(sName)
                    _up = mtempl.get(shift=realshift)
                    _down = mtempl.get(shift=-realshift)
                    if badtemp_ma(_up[0]) or badtemp_ma(_down[0]):
                        print("Skipping sample {}, scale systematic would be empty".format(sName))
                        continue
                    sample.setParamEffect(sys_scale, _up, _down, scale=1)

                if smear_syst and sName not in ['qcd', 'zll', 'wln', 'vvqq']:
                    _up = mtempl.get(smear=1 + SF[year]['smear_SF_ERR'])
                    _down = mtempl.get(smear=1 - SF[year]['smear_SF_ERR'])
                    if badtemp_ma(_up[0]) or badtemp_ma(_down[0]):
                        print("Skipping sample {}, scale systematic would be empty".format(sName))
                        continue
                    sample.setParamEffect(sys_smear, _up, _down)

                ch.addSample(sample)

            if not pseudo:
                data_obs = get_templ(f, region, 'data_obs',
                                     ptbin)[:-1]  # Don't pass variances
                if ptbin == 0 and region == "pass": print("Reading real data")

            else:
                raise NotImplementedError('Not up to date.')
                yields = []
                if 'qcd' not in include_samples:
                    include_samples = include_samples + ['qcd']
                for samp in include_samples:
                    if samp == "qcd" and opts.mockQCD and region == "pass":
                        _temp_yields = get_templ(
                            f, "fail", samp, ptbin)[0] * qcdeff * np.linspace(
                                0.8, 1.2, len(get_templ(f, "fail", samp, ptbin)[0]))
                    else:
                        _temp_yields = get_templ(f, region, samp, ptbin)[0]
                    if samp not in ['qcd', 'tqq', 'stqq'] and systs:
                        _temp_yields *= SF[year]['V_SF']
                    yields.append(_temp_yields)
                yields = np.sum(np.array(yields), axis=0)
                if throwPoisson:
                    yields = np.random.poisson(yields)

                data_obs = (yields, msd.binning, msd.name)
            ch.setObservation(data_obs)

            # drop bins outside rho validity
            ch.mask = mask

    if fitTF:
        if opts.transform:
            _transform = np.exp
        else:
            _transform = None
        degs = tuple([int(s) for s in opts.degs.split(',')])
        _basis = opts.basis.split(",")[1]
        if _basis == 'Bernstein':
            _inits = np.ones(tuple(n + 1 for n in degs))
        elif _basis == 'Chebyshev':
            _inits = np.zeros(tuple(n + 1 for n in degs))
            _inits[0,0] = 1
        else:
            raise ValueError("Basis ``{}`` not understood.".format(_basis))
        tf_dataResidual = rl.BasisPoly("tf{}_dataResidual".format(year),
                                       degs, ['pt', 'rho'], basis=_basis, init_params=_inits,
                                       limits=(0, 10), coefficient_transform=_transform)
        tf_dataResidual_params = tf_dataResidual(ptscaled, rhoscaled)
        if MCTF:
            tf_params = qcdeff * tf_MCtempl_params_final * tf_dataResidual_params
        else:
            tf_params = qcdeff * tf_dataResidual_params

        for ptbin in range(npt):
            failCh = model['ptbin{}fail{}'.format(ptbin, year)]
            passCh = model['ptbin{}pass{}'.format(ptbin, year)]

            qcdparams = np.array([
                rl.IndependentParameter(
                    'qcdparam{}_ptbin{}_msdbin{}'.format(year, ptbin, i), 0)
                for i in range(msd.nbins)
            ])
            initial_qcd = failCh.getObservation().astype(float)
            # was integer, and numpy complained about subtracting float from it
            for sample in failCh:
                initial_qcd -= sample.getExpectation(nominal=True)
            if np.any(initial_qcd < 0.):
                # raise ValueError("initial_qcd negative for some bins..", initial_qcd)
                warnings.warn("initial_qcd negative for some bins..", UserWarning)
                print(initial_qcd)
                warnings.warn("Negative bins will be forced positive", UserWarning)
                initial_qcd[0 > initial_qcd] = abs(initial_qcd[0 > initial_qcd])
            sigmascale = 10  # to scale the deviation from initial
            scaledparams = initial_qcd * (
                1 + sigmascale / np.maximum(1., np.sqrt(initial_qcd)))**qcdparams
            fail_qcd = rl.ParametericSample('ptbin{}fail{}_qcd'.format(ptbin, year),
                                            rl.Sample.BACKGROUND, msd, scaledparams)
            failCh.addSample(fail_qcd)
            pass_qcd = rl.TransferFactorSample('ptbin{}pass{}_qcd'.format(ptbin, year),
                                               rl.Sample.BACKGROUND,
                                               tf_params[ptbin, :], fail_qcd)
            passCh.addSample(pass_qcd)

    if muonCR:
        for ptbin in range(npt):
            failCh = model['ptbin{}fail{}'.format(ptbin, year)]
            passCh = model['ptbin{}pass{}'.format(ptbin, year)]
            tqqpass = passCh['tqq']
            tqqfail = failCh['tqq']
            stqqpass = passCh['stqq']
            stqqfail = failCh['stqq']
            sumPass = tqqpass.getExpectation(nominal=True).sum()
            sumFail = tqqfail.getExpectation(nominal=True).sum()
            sumPass += stqqpass.getExpectation(nominal=True).sum()
            sumFail += stqqfail.getExpectation(nominal=True).sum()
            tqqPF =  sumPass / sumFail
            tqqpass.setParamEffect(tqqeffSF, 1 * tqqeffSF)
            tqqfail.setParamEffect(tqqeffSF, (1 - tqqeffSF) * tqqPF + 1)
            tqqpass.setParamEffect(tqqnormSF, 1 * tqqnormSF)
            tqqfail.setParamEffect(tqqnormSF, 1 * tqqnormSF)
            stqqpass.setParamEffect(tqqeffSF, 1 * tqqeffSF)
            stqqfail.setParamEffect(tqqeffSF, (1 - tqqeffSF) * tqqPF + 1)
            stqqpass.setParamEffect(tqqnormSF, 1 * tqqnormSF)
            stqqfail.setParamEffect(tqqnormSF, 1 * tqqnormSF)

    # Fill in muon CR
    collapse = True
    if muonCR:
        for region in ['pass', 'fail']:
            ch = rl.Channel("muonCR{}{}".format(region, year))
            model.addChannel(ch)
            include_samples = [
                "qcd",
                "tqq",
                "stqq",
                "vvqq",
                'zll',
                'wln',
            ]

            for sName in include_samples:

                templ = get_templ(f_mu, region, sName, ptbin, muon=True)
                if collapse:
                    templ = one_bin(templ)
                stype = rl.Sample.BACKGROUND
                sample = rl.TemplateSample(ch.name + '_' + sName, stype, templ)

                if not systs:  # Need at least one
                    sample.setParamEffect(sys_lumi, 1.023)
                else:
                    sample.setParamEffect(sys_lumi, lumi_dict[year])
                    sample.setParamEffect(sys_lumi_correlated, lumi_correlated_dict[year])
                    if year != '2016':
                        sample.setParamEffect(sys_lumi_1718, lumi_1718_dict[year])
                    sample.setParamEffect(sys_eleveto, 1.005)
                    sample.setParamEffect(sys_tauveto, 1.005)

                    if sName in ["qcd"]:
                        continue

                    sys_names = [
                        'JES', 'JER', 'UES', 'mu_trigger', 'mu_isoweight', 'mu_idweight', 'btagEffStat', 'btagWeight', 'pileup_weight',
                        'Z_d2kappa_EW', 'Z_d3kappa_EW', 'd1kappa_EW', 'd1K_NLO', 'd2K_NLO', 'd3K_NLO',
                        'L1Prefiring',
                    ]
                    for sys_name in sys_names:
                        if (("NLO" in sys_name) or ("EW" in sys_name)) and not sName in ['zbb', 'zcc', 'zqq', 'wcq', 'wqq']:
                            continue
                        if ("Z_d" in sys_name) and sName not in ['zbb', 'zcc', 'zqq']:
                            continue
                        if opts.fast == 0:  # Convert to lnN for faster fitting
                            _sys_ef = shape_to_num(f_mu, region, sName, ptbin,
                                                    sys_name, mask, muon=True)
                            if _sys_ef is None:
                                continue
                            sample.setParamEffect(sys_shape_dict[sys_name], _sys_ef)
                        else:
                            _up = get_templ(f_mu,
                                            region,
                                            sName,
                                            ptbin,
                                            muon=True,
                                            syst=sys_name + "Up")
                            _dn = get_templ(f_mu,
                                            region,
                                            sName,
                                            ptbin,
                                            muon=True,
                                            syst=sys_name + "Down")
                            if _up is None or _dn is None:
                                continue
                            if badtemp_ma(_up[0]) or badtemp_ma(_dn[0]):
                                print("Skipping {} for sample {}, shape would be empty".format(sys_name, sName))
                                continue
                            if collapse:
                                _up = one_bin(_up)
                                _dn = one_bin(_dn)
                            sample.setParamEffect(sys_shape_dict[sys_name], _up[0], _dn[0])

                    if opts.mcstat and sName not in ['qcd'] and not collapse:
                        if opts.fast < 2:  # Convert to lnN for faster fitting
                            sample.autoMCStats(lnN=True)
                        else:
                            sample.autoMCStats(epsilon=1e-4)

                    if sName not in ["qcd", 'zll', 'wln']:
                        sample.scale(SF[year]['V_SF'])
                        sample.setParamEffect(
                            sys_veff, 1.0 + SF[year]['V_SF_ERR'] / SF[year]['V_SF'])

                ch.addSample(sample)

            if not pseudo:
                data_obs = get_templ(f_mu, region, 'data_obs', ptbin, muon=True)[:-1]
                if ptbin == 0 and region == "pass":
                    print("Reading real data")

            else:
                raise NotImplementedError('Not up to date.')
                yields = []
                for samp in include_samples:
                    _temp_yields = get_templ(f_mu, region, samp, ptbin, muon=True)[0]
                    yields.append(_temp_yields)
                yields = np.sum(np.array(yields), axis=0)
                if throwPoisson:
                    yields = np.random.poisson(yields)
                data_obs = (yields, msd.binning, msd.name)

            if collapse:
                data_obs = one_bin(data_obs)
            _nbinsmu = len(data_obs[0])

            ch.setObservation(data_obs)

        tqqpass = model['muonCRpass{}_tqq'.format(year)]
        tqqfail = model['muonCRfail{}_tqq'.format(year)]
        stqqpass = model['muonCRpass{}_stqq'.format(year)]
        stqqfail = model['muonCRfail{}_stqq'.format(year)]
        sumPass = tqqpass.getExpectation(nominal=True).sum()
        sumFail = tqqfail.getExpectation(nominal=True).sum()
        sumPass += stqqpass.getExpectation(nominal=True).sum()
        sumFail += stqqfail.getExpectation(nominal=True).sum()
        tqqPF = sumPass / sumFail
        tqqpass.setParamEffect(tqqeffSF, 1 * tqqeffSF)
        tqqfail.setParamEffect(tqqeffSF, (1 - tqqeffSF) * tqqPF + 1)
        tqqpass.setParamEffect(tqqnormSF, 1 * tqqnormSF)
        tqqfail.setParamEffect(tqqnormSF, 1 * tqqnormSF)
        stqqpass.setParamEffect(tqqeffSF, 1 * tqqeffSF)
        stqqfail.setParamEffect(tqqeffSF, (1 - tqqeffSF) * tqqPF + 1)
        stqqpass.setParamEffect(tqqnormSF, 1 * tqqnormSF)
        stqqfail.setParamEffect(tqqnormSF, 1 * tqqnormSF)

    with open("{}.pkl".format(model_name), "wb") as fout:
        pickle.dump(model, fout)

    model.renderCombine(model_name)

    conf_dict = vars(opts)
    # add info for F-test
    conf_dict['NBINS'] = np.sum(validbins)
    conf_dict['NBINSMU'] = _nbinsmu if muonCR else 0

    import json
    # Serialize data into file:
    json.dump(conf_dict,
              open("{}/config.json".format(model_name), 'w'),
              sort_keys=True,
              indent=4,
              separators=(',', ': '))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser.add_argument("--throwPoisson",
                        type=str2bool,
                        default='False',
                        choices={True, False},
                        help="If plotting data, redraw from poisson distribution")
    
    parser.add_argument("--basis",
                        type=str,
                        default='Bernstein,Bernstein',
                        help="Comma separated bases for TF fits (MC,Residual)."
                             "Choose from {'Bernstein' | 'Chebyshev'}"
                        )

    parser.add_argument("--fitTF",
                        type=str2bool,
                        default='True',
                        choices={True, False},
                        help="Fit TF for QCD")

    parser.add_argument("--MCTF",
                        type=str2bool,
                        default='True',
                        choices={True, False},
                        help="Fit QCD in MC first")

    parser.add_argument("--muCR", "--muonCR",
                        type=str2bool,
                        default='True',
                        choices={True, False},
                        help="Include muonCR to constrain ttbar")

    parser.add_argument("--scale",
                        type=str2bool,
                        default='True',
                        choices={True, False},
                        help="Include scale systematic")

    parser.add_argument("--smear",
                        type=str2bool,
                        default='True',
                        choices={True, False},
                        help="Include smear systematics")

    parser.add_argument("--systs",
                        type=str2bool,
                        default='True',
                        choices={True, False},
                        help="Include all systematics (separate from scale/smear)")

    parser.add_argument("--fast",
        type=int,
        default=0,
        choices=[0, 1, 2],
        help="0: Convert all shapes and mcstat to lnN (except scale/smear)"
             "1: Convert only mcstat shapes"
             "2: Don't convert anything"
        )

    parser.add_argument("--justZ",
                        type=str2bool,
                        default='False',
                        choices={True, False},
                        help="Only run Z sample with QCD")

    parser.add_argument("--justHZ",
                        type=str2bool,
                        default='False',
                        choices={True, False},
                        help="Only run H and Z sample with QCD")

    parser.add_argument("--year", type=int, default=2017, help="Year")

    parser.add_argument("--matched",
                        type=str2bool,
                        default='False',
                        choices={True, False},
                        help=("Use matched/unmatched templates"
                              "(w/o there is some W/Z/H contamination from QCD)"))

    parser.add_argument("--mcstat",
                        type=str2bool,
                        default='True',
                        choices={True, False},
                        help="Include mcstat unc")

    parser.add_argument("-t", "--templates", "--t", type=str, dest='templates', default=None, required=True, help="Primary templates")

    parser.add_argument("--mutemplates", "--mut",  "--tm", type=str, default=None, help="Muon templates")

    parser.add_argument('-o', "--model", type=str, default=None, help="Model directory")

    pseudo = parser.add_mutually_exclusive_group(required=True)
    pseudo.add_argument('--data', action='store_false', dest='pseudo')
    pseudo.add_argument('--MC', action='store_true', dest='pseudo')

    parser.add_argument('--unblind', action='store_true', dest='unblind')
    parser.add_argument('--higgs',
                        action='store_true',
                        dest='runhiggs',
                        help="Set Higgs as signal instead of z")
    parser.add_argument('--both',
                        action='store_true',
                        dest='runboth',
                        help="Both Z and H signals")
    parser.add_argument('--bb', action='store_true', dest='runbb')

    parser.add_argument('--transform', action='store_true', dest='transform')

    parser.add_argument(
        '--mockQCD',
        action='store_true',
        dest='mockQCD',
        help="Replace true pass QCD with scaled true fail QCD in pseudo data")

    parser.add_argument("--degs",
                        type=str,
                        default='1,2',
                        help="Polynomial degrees in the shape 'pt,rho' e.g. '2,2'")

    parser.add_argument("--degsMC",
                        type=str,
                        default='1,2',
                        help="Polynomial degrees in the shape 'pt,rho' e.g. '2,2'")

    args = parser.parse_args()
    print("Running with options:")
    print("    ", args)
    if args.mutemplates is None:
        args.mutemplates = args.templates.replace("templates_", "templatesmuCR_")

    dummy_rhalphabet(pseudo=args.pseudo,
                     throwPoisson=args.throwPoisson,
                     MCTF=args.MCTF,
                     scale_syst=args.scale,
                     smear_syst=args.smear,
                     systs=args.systs,
                     justZ=args.justZ,
                     blind=(not args.unblind),
                     runhiggs=args.runhiggs,
                     runboth=args.runboth,
                     fitTF=args.fitTF,
                     year=str(args.year),
                     muonCR=args.muCR,
                     opts=args)
