import os
import argparse
import random
import string

parser = argparse.ArgumentParser(description='')
parser.add_argument('base')
parser.add_argument('alt')
parser.add_argument('--gen', action='store_true', help='')
parser.add_argument('--fit', action='store_true', help='')
parser.add_argument('--collate', action='store_true', help='')
# Condor (get/fit)
parser.add_argument('--condor', action='store_true', help='')
# Optional
parser.add_argument('--tf', action='store_true', help='Toys Freq, freeze systs')
parser.add_argument('--set', choices=['2017CC', '2017tt'], help='', default=None)

args = parser.parse_args()

fit_call = "combineTool.py -M FitDiagnostics -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1  "
gen_call = "combineTool.py -M GenerateOnly  -d model_combined.root --cminDefaultMinimizerStrategy 0 --saveToys "
cond_str = """ --job-mode condor --sub-opts='+JobFlavour = "espresso"' --task-name task_{}"""
tf_srt = " --toysFrequentist --freezeParameters allConstrainedNuisances "

def rnd_suffix(N):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(N))

# steps - generate, fit1, fit 2, collate
# seed = "1"
seed = "0:20:1"
ntoys = 25 
parse_seed = range(*[int(a) for a in seed.split(":")])

# os.chdir('Tau17')
cwd = os.getcwd()
base_dir = 'TauRun2'
alt_dir = 'TauRun2Alt'
base_dir = args.base
alt_dir = args.alt

SPdict = {
    '2017CC': 'tf2017_MCtempl_deco0=-1.4921e+00,tf2017_MCtempl_deco1=-4.3915e-01,tf2017_MCtempl_deco2=4.6416e-01,tf2017_MCtempl_deco3=2.2583e-01,tf2017_MCtempl_deco4=-4.5086e-01,tf2017_MCtempl_deco5=8.3076e-02,tf2017_dataResidual_pt_par0_rho_par0=1.0227e+00',
    '2017tt': 'tf2017_MCtempl_deco0=1.6164e+00,tf2017_MCtempl_deco1=-4.8994e-02,tf2017_MCtempl_deco2=-1.5402e-01,tf2017_MCtempl_deco3=1.2641e-01,tf2017_MCtempl_deco4=4.5620e-01,tf2017_MCtempl_deco5=-2.0346e-01,tf2017_dataResidual_pt_par0_rho_par0=3.7624e-01',
}

for ss in [0, 1, 5] + range(10, 105, 5):
    tname = "{set}{ss}".format(set="" if args.set is None else "set", ss=ss)
    if args.gen:
        os.chdir(alt_dir)
        cmd = gen_call + "--expectSignal {ss} -n {tname} -t {ntoys} -s {seed}".format(ss=ss, ntoys=ntoys, seed=seed, tname=tname)
        if args.set is not None:
            cmd += ' --setParameters {} '.format(SPdict[args.set])
        if args.tf:
            cmd += tf_srt
        if args.condor:
            cmd += cond_str.format(rnd_suffix(6))
        os.system(cmd)
        os.chdir(cwd)

    if args.fit:
        os.chdir(alt_dir)
        for _s in parse_seed:
            cmd = fit_call + " --toysFile higgsCombine{tname}.GenerateOnly.mH120.{_s}.root".format(ss=ss, _s=_s, tname=tname)
            cmd += " --expectSignal {ss} -n {tname} -t {ntoys} -s {_s}".format(ss=ss, ntoys=ntoys, _s=_s, tname=tname)
            if args.tf:
                cmd += tf_srt
            if args.condor:
                cmd += cond_str.format(rnd_suffix(6))
            os.system(cmd)
        os.chdir(cwd)

        os.chdir(base_dir)
        for _s in parse_seed:
            cmd = fit_call + " --toysFile ../{alt_dir}/higgsCombine{tname}.GenerateOnly.mH120.{_s}.root".format(ss=ss, _s=_s, alt_dir=alt_dir, tname=tname)
            cmd += " --expectSignal {ss} -n {tname} -t {ntoys} -s {_s}".format(ss=ss, ntoys=ntoys, _s=_s, tname=tname)
            if args.tf:
                cmd += tf_srt
            if args.condor:
                cmd += cond_str.format(rnd_suffix(6))
            os.system(cmd)
        os.chdir(cwd)

    if args.collate:
        os.chdir(alt_dir)
        if args.set is not None:
            out_name = "bias_set"
        elif args.tf is not None:
            out_name = "bias_toysF"
        else:
            out_name = "bias" 
        cmd = 'hadd -f {out_name}_alt_{ss}_{_s}.root higgsCombine{ss}.FitDiagnostics.mH120.*.root'.format(ss=ss, _s=seed.replace(":", "-"), out_name=out_name)
        os.system(cmd)
        os.chdir(cwd)

        os.chdir(cwd)
        os.chdir(base_dir)
        cmd = 'hadd -f {out_name}_base{alt}_{ss}_{_s}.root higgsCombine{ss}.FitDiagnostics.mH120.*.root'.format(ss=ss, _s=seed.replace(":", "-"), alt=alt_dir, out_name=out_name)
        os.system(cmd)
        os.chdir(cwd)


