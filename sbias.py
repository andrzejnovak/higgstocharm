import os
import sys
import argparse
import random
import string

syslista = [
 'allConstrainedNuisances',
 '-',
 r"rgx{CMS_gghcc_s.*}",
 r'rgx{CMS_eff_bb_.*}',
 r'rgx{CMS_eff_cc_.*}',
 r'rgx{CMS_eff_w_.*}',
 r'rgx{CMS_gghcc_btagEffStat_.*}',
 r'rgx{CMS_gghcc_btagWeight_.*}',
 r'rgx{CMS_gghcc_e_veto_.*}',
 r'rgx{CMS_gghcc_ggHpt}',
 r'rgx{CMS_gghcc_m_veto_.*}',
 r'rgx{CMS_gghcc_scale_.*}',
 r'rgx{CMS_gghcc_smear_.*}',
 r'rgx{CMS_gghcc_trigger_.*}',
 r'rgx{CMS_gghcc_veff_.*}',
 r'rgx{CMS_gghcc_wznormEW}',
 r'rgx{CMS_gghcc_znormEW}',
 r'rgx{CMS_gghcc_znormQ}',
 r'rgx{CMS_lumi}',
 r'rgx{CMS_res_j_.*}',
 r'rgx{CMS_scale_j_.*}',
 r'rgx{CMS_ues_j_.*}',
 r'rgx{tf.*_MCtempl_deco0}',
 r'rgx{tf.*_MCtempl_deco1}',
 r'rgx{tf.*_MCtempl_deco2}',
 r'rgx{tf.*_MCtempl_deco3}',
 r'rgx{tf.*_MCtempl_deco4}',
 r'rgx{tf.*_MCtempl_deco5}',
 r'rgx{tf.*_dataResidual_pt_par0_rho_par0}',
 r'rgx{tf.*_dataResidual_pt_par1_rho_par0}',
 r'rgx{tf.*_dataResidual_pt_par0_rho_par1}',
 r'rgx{tf.*_dataResidual_pt_par1_rho_par1}',
 r'rgx{tqqeffSF_.*}',
 r'rgx{tqqnormSF_.*}'
 ]

syslistn = [
 'allConstrainedNuisances',
 'None',
 "scale_smear",
 'CMS_eff_bb_',
 'CMS_eff_cc_',
 'CMS_eff_w_',
 'CMS_gghcc_btagEffStat_',
 'CMS_gghcc_btagWeight_',
 'CMS_gghcc_e_veto_',
 'CMS_gghcc_ggHpt',
 'CMS_gghcc_m_veto_',
 'CMS_gghcc_scale_',
 'CMS_gghcc_smear_',
 'CMS_gghcc_trigger_',
 'CMS_gghcc_veff_',
 'CMS_gghcc_wznormEW',
 'CMS_gghcc_znormEW',
 'CMS_gghcc_znormQ',
 'CMS_lumi',
 'CMS_res_j_',
 'CMS_scale_j_',
 'CMS_ues_j_',
 'tf_MCtempl_deco0',
 'tf_MCtempl_deco1',
 'tf_MCtempl_deco2',
 'tf_MCtempl_deco3',
 'tf_MCtempl_deco4',
 'tf_MCtempl_deco5',
 'tf_dataResidual_pt_par0_rho_par0',
 'tf_dataResidual_pt_par1_rho_par0',
 'tf_dataResidual_pt_par0_rho_par1',
 'tf_dataResidual_pt_par1_rho_par1',
 'tqqeffSF_',
 'tqqnormSF_',
]

syslistn = [
 'allConstrainedNuisances',
 'None',
#  "scale_smear",
]
syslista = [
 'allConstrainedNuisances',
 '-',
#  r"rgx{CMS_gghcc_s.*}",
]


parser = argparse.ArgumentParser(description='')
parser.add_argument('dnames', nargs='*', default=os.getcwd())
parser.add_argument('--gen', action='store_true', help='')
parser.add_argument('--fit1', action='store_true', help='')
parser.add_argument('--fit2', action='store_true', help='')
parser.add_argument('--collate', action='store_true', help='')
parser.add_argument('--condor', action='store_true', help='')
parser.add_argument('--mc', action='store_true', help='')
args = parser.parse_args()

base_dir = args.dnames[0]
try:
    alt_dir = args.dnames[1]
except:
    alt_dir = None

fit_call = "combineTool.py -M FitDiagnostics -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 --setParameters z=1"
gen_call = "combineTool.py -M GenerateOnly  -d model_combined.root --cminDefaultMinimizerStrategy 0 --saveToys --setParameters z=1"
if not args.mc:
    fit_call += " --toysFrequentist "
    gen_call += " --toysFrequentist "
cond_str = """ --job-mode condor --sub-opts='+JobFlavour = "espresso"' --task-name task_{}"""

def rnd_suffix(N):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(N))

# seed = "0:19:1"
# ntoys = 25
seed = "0:19:1"
ntoys = 10
_pseed = [int(a) for a in seed.split(":")]
_pseed[1] += 1
parse_seed = range(*_pseed)

cwd = os.getcwd()
base_dir = os.path.join(cwd, base_dir)
alt_dir = os.path.join(cwd, alt_dir) if alt_dir is not None else None

ss = 100
sysa = "-"
sys = "-"
# for sys in syslist:
for sysa, sys in zip(syslista, syslistn):
    print("XXXXXX")
    print(sys)
    print("XXXXXX")
    if args.gen:
        os.chdir(base_dir)
        cmd = gen_call + "--expectSignal {ss} -n {sys} -t {ntoys} -s {seed}".format(ss=ss, ntoys=ntoys, seed=seed, sys=sys)
        cmd += " --freezeParameters {} ".format(sysa)
        if args.condor:
            cmd += cond_str.format(rnd_suffix(6))
        os.system(cmd)
        os.chdir(cwd)

    if args.fit1:
        os.chdir(base_dir)
        for _s in parse_seed:
            cmd = fit_call + " --toysFile higgsCombine{sys}.GenerateOnly.mH120.{_s}.root".format(ss=ss, _s=_s, sys=sys)
            cmd += " --expectSignal {ss} -n {sys} -t {ntoys} -s {_s}".format(ss=ss, ntoys=ntoys, _s=_s, sys=sys)
            cmd += " --freezeParameters {} ".format(sysa)
            if args.condor:
                cmd += cond_str.format(rnd_suffix(6))
            os.system(cmd)
        os.chdir(cwd)


    if args.collate:
        os.chdir(cwd)
        os.chdir(base_dir)
        oname = 'nfdz'
        if args.mc:
            oname += "MC"
        cmd = 'hadd -f '+oname+'_base_{ss}_{sys}.root higgsCombine{sys}.FitDiagnostics.mH120.*.root'.format(ss=ss, sys=sys)
        os.system(cmd)
        os.chdir(cwd)


