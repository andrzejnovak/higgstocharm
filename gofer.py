import os
import sys
import time
import argparse

BASE = os.getcwd()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make plots for ftest/gof files.')
    parser.add_argument('dname', nargs='?', default=os.getcwd())
    parser.add_argument('--make', action='store_true', help='')
    parser.add_argument('--blind', action='store_true', help='')
    parser.add_argument('--mc', action='store_true', help='')

    parser.add_argument('--build', action='store_true', help='')
    parser.add_argument('--run', action='store_true', help='')
    parser.add_argument('--condor', action='store_true', help='')
    parser.add_argument('--collect', action='store_true', help='')
    parser.add_argument('--plot', action='store_true', help='')
    parser.add_argument('--year', type=str, default=None, help="Year to display on plots.")
    parser.add_argument('-p', '--param', type=str, choices=['bern', 'cheby', 'exp'], default='bern',
                        help="Parametrization")
    args = parser.parse_args()

    dname = args.dname
    print("X", dname)
    start = time.time()

    if args.param == 'cheby':
        basis = ' --basis Bernstein,Chebyshev'
    elif args.param == 'exp':
        basis = ' --basis Bernstein,Bernstein --transform '
    else:
        basis = " "

    if args.make:
        # Made workspaces
        for pt in range(0, 4):
            for rho in range(0, 4):
                if args.mc:
                    os.system(
                        "python new_Hxx.py --MC --year {year} --templates tau/templates_ref{yearshort}_CC.root -o {dname}{p}{r} --degs {p},{r} "
                        " --justZ True --muCR False --MCTF False --muCR False --systs False --scale False --smear False --unblind "
                        .format(year=args.year, yearshort=args.year[-2:], dname=dname, p=str(pt), r=str(rho))
                        + basis
                    ) 
                else:
                    os.system(
                        "python new_Hxx.py --data --year {year} --templates reft/templates_ref{yearshort}_CC.root -o {dname}{p}{r} --degs {p},{r} "
                        "--mutemplates reft/templatesmuCR_ref{yearshort}_CC.root  --muCR True"
                        .format(year=args.year, yearshort=args.year[-2:], dname=dname, p=str(pt), r=str(rho))
                        + (" --unblind " if not args.blind else "")
                        + basis
                    )

    N = 3
    if args.build:
        N = 4
    for pt in range(0, N):
        for rho in range(0, N):
            path = "{}/{}{}{}".format(BASE, dname, pt, rho)
            print("XXXXXXXXXXXXXXXX")
            print("Working in:")
            print("    ", path)
            print("XXXXXXXXXXXXXXXX")
            os.chdir(path)

            if args.build:
                os.system("bash build.sh")

            if args.run:
                # os.system(
                #     "combine -M FitDiagnostics --expectSignal 1 -d model_combined.root  --cminDefaultMinimizerStrategy 0  --robustFit=1 "
                #     " -t -1 --toysFrequentist"
                #     " --redefineSignalPOIs z"
                #     )
                
                ### F-tests\
                # continue
                _seed = "1:10:1"
                _toys = '50'
                os.system(
                    "python ../workflows.py -w model_combined.root --method FTest --altmodel ../{}{}{}/model_combined.root "
                    " -t {} --seed {} ".format(dname, pt + 1, rho, _toys, _seed)
                    + (" --condor" if args.condor else "")
                    + (" --collect" if args.collect else "")
                    + " --freezeParameters r,z"
                    # " --debug"
                    )
                
                os.system(
                    "python ../workflows.py -w model_combined.root --method FTest --altmodel ../{}{}{}/model_combined.root "
                    " -t {} --seed {} ".format(dname, pt, rho + 1, _toys, _seed)
                    + (" --condor" if args.condor else "")
                    + (" --collect" if args.collect else "")
                    + " --freezeParameters r,z"
                    # " --debug"
                    )

                os.system(
                    "python ../workflows.py -w model_combined.root --method FTest --altmodel ../{}{}{}/model_combined.root "
                    " -t {} --seed {} ".format(dname, pt + 1, rho + 1, _toys, _seed)
                    + (" --condor" if args.condor else "")
                    + (" --collect" if args.collect else "")
                    + " --freezeParameters r,z"
                    # " --debug"
                    )

                ### Impacts
                # os.system("combineTool.py -M Impacts -d model_combined.root -m 125 --doInitialFit --robustFit 1 --setParameterRanges r=-1,5 --cminDefaultMinimizerStrategy 0 --X-rtd FITTER_DYN_STEP --expectSignal 1 -t -1 --toysFrequentist")
                # os.system("combineTool.py -M Impacts -d model_combined.root -m 125 --doFits --robustFit 1 --allPars --setParameterRanges r=-1,5 --expectSignal 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic" +
                #            """ --job-mode condor --sub-opts='+JobFlavour = "workday"' -t -1 --toysFrequentist """ +
                #            " --task-name 'tf%s%s' --exclude 'rgx{qcdparams*}'" % (pt, rho) )
                # os.system(
                #     "combineTool.py -M Impacts -d model_combined.root -m 125 --allPars -o impacts_toys.json"
                # )
                # os.system("plotImpacts.py -i impacts_toys.json -o impacts_toys")

    # # Collect plots
    if args.plot:
        os.chdir(BASE)
        os.mkdir('plots_{}'.format(dname))
        for pt in range(0, 3):
            for rho in range(0, 3):
                print("X", pt, rho)
                # os.system("python fplots.py {}{}{}  --year {} --out {}".format(dname, str(pt), str(rho), args.year, 'plots_{}'.format(dname)))
                os.system("python plot_ftest.py {}{}{}  --year {} -d {}".format(dname, str(pt), str(rho), args.year, 'plots_{}'.format(dname)))

    print("Done")
    end = time.time()
    print("TIME:", time.strftime("%H:%M:%S", time.gmtime(end - start)))

