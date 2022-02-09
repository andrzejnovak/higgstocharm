# from distutils import command
import os
import subprocess
from multiprocessing import Process
import argparse
import time

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', default="FTests/mconly16")
    parser.add_argument('-o', default="Ftests/Fouts")
    parser.add_argument('--outplots', default="Ftests/plots")

    parser_mc = parser.add_mutually_exclusive_group(required=True)
    parser_mc.add_argument('--data', action='store_false', dest='mc')
    parser_mc.add_argument('--mc', action='store_true', dest='mc')
    parser.add_argument('--year', type=str, choices=["2016", "2017", "2018"], required=True, help="Year to display.")
    parser.add_argument('--blind', action='store_true', help='')
    parser.add_argument('--param', type=str, choices=['bern', 'cheby', 'exp'], default='bern')

    parser.add_argument('--make', action='store_true', help='')
    parser.add_argument('--build', action='store_true', help='')

    parser.add_argument('--run', action="store_true")
    parser.add_argument('--toys', '-t', default="10")
    parser.add_argument('--seed', '-s', default="1")
    parser.add_argument('--base', action="store_true")
    parser.add_argument('--gen', action="store_true")
    parser.add_argument('--fits', action="store_true")
    parser.add_argument('--all', action="store_true")
    parser.add_argument('--condor', action="store_true")

    parser.add_argument('--plot', action="store_true")
    parser.add_argument('-p', action="store_true", help="parallel")
    parser.add_argument('--debug', action="store_true", help="Just print")

    args = parser.parse_args()

    start = time.time()
    commands = []

    rng_pt = 4
    rng_rho = 4

    if args.param == 'cheby':
        basis = ' --basis Bernstein,Chebyshev'
    elif args.param == 'exp':
        basis = ' --basis Bernstein,Bernstein --transform '
    else:
        basis = " "
    
    if args.make:
        # Made workspaces
        for pt in range(0, rng_pt + 1):
            for rho in range(0, rng_rho + 1):
                if args.mc:
                    cmd = (
                        "python new_Hxx.py --MC --year {year} --templates temps/templates_unblind{yearshort}_CC.root -o {dname}{p}{r} --degs {p},{r} "
                        " --justZ True --muCR False --MCTF False --muCR False --systs False --scale False --smear False --unblind "
                        .format(year=args.year, yearshort=args.year[-2:], dname=args.d, p=str(pt), r=str(rho))
                        + basis
                    ) 
                else:
                    cmd = (
                        "python new_Hxx.py --data --year {year} --templates temps/templates_unblind{yearshort}_CC.root -o {dname}{p}{r} --degs {p},{r} "
                        "--mutemplates temps/templatesmuCR_unblind{yearshort}_CC.root  --muCR True --degsMC 0,2 "
                        .format(year=args.year, yearshort=args.year[-2:], dname=args.d, p=str(pt), r=str(rho))
                        + (" --unblind " if not args.blind else "")
                        + basis
                    )
                commands.append(cmd)

    if args.build:
        for pt in range(0, rng_pt + 1):
            for rho in range(0, rng_rho + 1):
                _refdir = os.path.realpath(os.getcwd())
                path = "{}{}{}".format(args.d, str(pt), str(rho))
                print("XXXXXXXXXXXXXXXX")
                print("Working in:")
                print("    ", path)
                print("XXXXXXXXXXXXXXXX")
                os.chdir(path)

                os.system("bash build.sh")
                os.chdir(_refdir)


    if args.all:
        args.base = True
        args.gen = True
        args.fits = True

    
    if args.run:
        run_cfg = (" --base " if args.base else "") + (" --gen " if args.gen else "") + (" --fits " if args.fits else "")
        if args.condor:
            run_cfg += " --condor "
        toy_cfg = " -t {} -s {} ".format(args.toys, args.seed)
        for i in range(4):
            for j in range(4):
                base_cmd = "python custom.py -d {} {} ".format(args.o, "--mc" if args.mc else "--data")
                ws_base = " -w {base}{i}{j}/model_combined.root "
                ws_alt = " -a {base}{i}{j}/model_combined.root "

                cmd = base_cmd + ws_base.format(base=args.d, i=i, j=j) + ws_alt.format(base=args.d, i=i + 1, j=j) + toy_cfg + run_cfg
                commands.append(cmd)
                cmd = base_cmd + ws_base.format(base=args.d, i=i, j=j) + ws_alt.format(base=args.d, i=i, j=j + 1) + toy_cfg + run_cfg
                commands.append(cmd)
                cmd = base_cmd + ws_base.format(base=args.d, i=i, j=j) + ws_alt.format(base=args.d, i=i + 1, j=j + 1) + toy_cfg + run_cfg
                commands.append(cmd)


    if args.plot:
        base_cmd = "python new_plot_ftests.py -o {} --year {} {} ".format(args.outplots, args.year, "--mc" if args.mc else "--data")
        if args.mc:
            base_cmd += " --qplots "
        for i in range(4):
            for j in range(4):
                cmd = base_cmd + " --degs {},{} ".format(i, j) + " --degsalt {},{} ".format(i, j + 1)
                cmd += " -d {}/ftest_{}{}_{}{}".format(args.o, i, j, i, j + 1)
                commands.append(cmd)
                cmd = base_cmd + " --degs {},{} ".format(i, j) + " --degsalt {},{} ".format(i + 1, j)
                cmd += " -d {}/ftest_{}{}_{}{}".format(args.o, i, j, i + 1, j)
                commands.append(cmd)
                cmd = base_cmd + " --degs {},{} ".format(i, j) + " --degsalt {},{} ".format(i + 1, j + 1)
                cmd += " -d {}/ftest_{}{}_{}{}".format(args.o, i, j, i + 1, j + 1)
                commands.append(cmd)

    if args.debug:
        for cmd in commands:
            print(cmd)
        import sys
        sys.exit()
    processes = []
    for cmd in commands:
        if args.p:
            processes.append(subprocess.Popen(cmd, shell=True))
        else:
            processes.append(subprocess.Popen(cmd, shell=True).wait())

    while sum([p.wait() is not None for p in processes]) < len(processes):
        try:
            time.sleep(1)
            print([p.poll() is not None for p in processes])
            print([p.wait() for p in processes])
        except KeyboardInterrupt:
            term = [p.terminate() for p in processes]

    print("TIME:", time.strftime("%H:%M:%S", time.gmtime(time.time()-start)))

# python submit_ftests.py -d FTests/mconly16 -o Ftests/Fouts


