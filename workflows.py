#!/usr/bin/env python
from __future__ import print_function
import os, argparse, ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

import sys
import re
# sys.path.append('/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/python')
# import cms_style
# cms_style.extra_text="Preliminary Simulation"
# cms_style.cms_style()
cms_logo = False
global silence
silence = False


def exec_bash(command='echo "hello world"', debug=False):
    global silence
    if (not silence):
        print(command)
    if (not debug):
        os.system(command)
    return """%s\n""" % command


def get_parameters(args, query, workspace):
    w = ROOT.TFile(args.workspace, "READ").Get('w')
    allVars = w.allVars().contentsString().split(',')
    parameters = []
    for var in allVars:
        if '_In' in var:
            continue
        if query in var:
            parameters.append(var)
    return parameters


rMIN = -1.5
rMAX = 1.5


class CombineWorkflows:
    def __init__(self, parser=None):
        # print('Nothing to do here. This class is just a wrapper for some worflow methods using combine')
        self.methods = [
            func for func in dir(CombineWorkflows)
            if (callable(getattr(CombineWorkflows, func)) and '__' not in func)
        ]
        self.methods.remove('write_wrapper')
        # self.methods.remove('method')
        self.externToys = False
        self.justplots = False
        self.skipplots = True
        self._poi = ""
        self.workspace = ""
        self.altmodel = None
        self.workers = 1
        self.condor = False
        self._freezeParameters = ""
        self.lumi = 41.8
        self.name = ""
        self.seed = 123456
        self.toys = 50
        self.algo = "saturated"
        self.extraOptions = ""
        self.job_index = 0
        self.rhalphdir = os.getcwd(
        )  #"/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/rhalph"
        self.toysOptions = "--toysFrequentist"
        self.combineCMSSW = self.rhalphdir + '/CMSSW_10_2_13'
        self.modeldir = ""
        self.workspace = ""

        self._method = ""
        self.parser = parser,

        def dummyMethod(debug=True):
            raise BaseException(
                "You have not selected a CombineWorkflow method! Choose from: " +
                ", ".join(self.methods))

        self.combineString = dummyMethod


    @property
    def workspace(self):
        return self._workspace

    @workspace.setter
    def workspace(self, w):
        if (isinstance(w, str)):
            self._workspace = os.path.abspath(w)
            self.modeldir = self.workspace.replace(self.workspace.split('/')[-1], '')
        elif (isinstance(w, list)):
            self._workspace = np.array([os.path.abspath(iw) for iw in w], dtype=object)
            self.modeldir = np.array([os.path.dirname(iw) for iw in self.workspace],
                                     dtype=object)

    @property
    def altmodel(self):
        return self._altmodel

    @altmodel.setter
    def altmodel(self, w):
        if (isinstance(w, str)):
            self._altmodel = os.path.abspath(w)
        elif (w is None):
            self._altmodel = None

    @property
    def condor(self):
        return self._condor

    @condor.setter
    def condor(self, b):
        if (isinstance(b, bool)):
            self._altmodel = b
        else:
            raise TypeError("`condor` has to be a bool")

    @property
    def freezeParameters(self):
        if (self._freezeParameters == ''):
            return ""
        return '--freezeParameters ' + ','.join(self._freezeParameters)

    @freezeParameters.setter
    def freezeParameters(self, pars):
        if (isinstance(pars, str)):
            self._freezeParameters = pars.split(",")
        elif (isinstance(pars, list)):
            self._freezeParameters = pars
        else:
            raise TypeError(
                "freezeParameters must be either string (',' as delimiter) or list!\nYou provided a ",
                type(pars))

    @property
    def POI(self):
        return self._poi

    @POI.setter
    def POI(self, pois):
        if (isinstance(pois, str)):
            if (pois.lower() == "r"):
                self._poi = ""
            elif ('--redefineSignalPOIs ' in pois):
                self._poi = pois
            else:
                self._poi = '--redefineSignalPOIs ' + pois
        elif (isinstance(pois, list)):
            self._poi = '--redefineSignalPOIs ' + ','.join(pois)
            return
        else:
            raise TypeError(
                "POI must be either string (',' as delimiter) or list!\nYou provided a ",
                type(pois))

    @property
    def method(self):
        return self._method

    @method.setter
    def method(self, method_name):
        self._method = method_name
        self.combineString = getattr(self, method_name)

    def write_wrapper(self):
        global silence
        silence = True
        pathCMSSW = os.path.realpath(self.combineCMSSW)
        # pathCMSSW = os.path.realpath(self.rhalphdir + '/CMSSW_10_2_13')

        if (not os.path.isdir(self.modeldir)):
            os.makedirs(self.modeldir)
        wrapper_name = self.modeldir + '/wrapper.sh'
        with open(wrapper_name, 'w') as wrapper:
            wrapper.write("#!/bin/bash\n")
            wrapper.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
            wrapper.write("cd " + pathCMSSW + "/src/\n")
            wrapper.write("eval `scramv1 runtime -sh`\n")
            wrapper.write(self.combineString(debug=True))
            wrapper.close()
        os.system('chmod u+x ' + wrapper_name)
        silence = False

    def diagnostics(self, debug=True):
        command_string = """#FitDiagnostics Workflow\n"""
        command_string += exec_bash("cd " + self.modeldir + "\n", debug)
        command_string += exec_bash("source build.sh\n", debug)
        command_string += exec_bash(
            "combine -M FitDiagnostics {WORKSPACE} {POI} --saveShapes {EXTRA}".format(
                WORKSPACE=self.workspace, POI=self.POI, EXTRA=self.extraOptions), debug)
        command_string += exec_bash(
            "PostFitShapesFromWorkspace -w {WORKSPACE} -o {MODELDIR}fit_shapes.root --postfit --sampling -f {MODELDIR}fitDiagnostics.root:fit_s"
            .format(WORKSPACE=self.workspace, MODELDIR=self.modeldir), debug)
        return command_string

    def GOF(self, debug=True, merge=False):
        command_string = """#GOF test\n"""
        command_string += exec_bash("cd " + self.modeldir, debug)
        command_string += exec_bash("source build.sh", debug)

        if (not self.justplots and not merge):
            command_string += exec_bash(
                "combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --algo={ALGO} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Baseline\""
                .format(WORKSPACE=self.workspace,
                        FREEZEPARAMS=self.freezeParameters,
                        EXTRA=self.extraOptions,
                        NAME=self.name,
                        POI=self.POI,
                        ALGO=self.algo), debug)
            if (self.externToys):
                command_string += exec_bash(
                    "combine -M GenerateOnly -d {WORKSPACE} {POI} -m 0 -t {NTOYS} --toysFrequentist --saveToys -n {NAME} --seed {SEED}"
                    .format(WORKSPACE=self.workspace,
                            NAME=self.name,
                            POI=self.POI,
                            SEED=self.seed,
                            NTOYS=self.toys), debug)
                command_string += exec_bash(
                    "combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --algo={ALGO}  {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} --toysFrequentist  --toysFile higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --seed {SEED}"
                    .format(WORKSPACE=self.workspace,
                            FREEZEPARAMS=self.freezeParameters,
                            EXTRA=self.extraOptions,
                            NAME=self.name,
                            POI=self.POI,
                            SEED=self.seed,
                            NTOYS=self.toys,
                            ALGO=self.algo), debug)
            else:
                command_string += exec_bash(
                    "combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --algo={ALGO}  {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} {TOYSOPT} --seed {SEED}"
                    .format(WORKSPACE=self.workspace,
                            FREEZEPARAMS=self.freezeParameters,
                            EXTRA=self.extraOptions,
                            NAME=self.name,
                            POI=self.POI,
                            SEED=self.seed,
                            NTOYS=self.toys,
                            TOYSOPT=self.toysOptions,
                            ALGO=self.algo), debug)
        command_string += exec_bash(
            'python {RHALPHDIR}/CombinePlotter.py --method plot_gof_result --parameter "higgsCombine{NAME}Baseline.GoodnessOfFit.mH0.root;higgsCombine.{NAME}GoodnessOfFit.mH0.{SEED}.root;{ALGO};{LUMI}"'
            .format(RHALPHDIR=self.rhalphdir,
                    NAME=self.name,
                    SEED=self.seed,
                    ALGO=self.algo,
                    LUMI=self.lumi), debug)
        return command_string

    def FTestBatch(self, debug=True):
        command_string = "#FTest\n"
        qcd_fit = "qcdmodel" in self.modeldir
        import glob, json
        configs = [
            os.path.abspath(json_path)
            for json_path in glob.glob(self.modeldir + ("../" if qcd_fit else "") +
                                       "../*.json")
        ]
        alt_model_dirs = [
            os.path.dirname(config) + '/' + json.load(open(config))["ModelName"] +
            ("/qcdmodel/" if qcd_fit else "/") for config in configs
        ]
        alt_model_dirs.remove(self.modeldir)

        command_string += exec_bash("cd " + self.modeldir, debug)
        command_string += exec_bash("source build.sh", debug)

        GOF_extra = self.extraOptions + (" --fixedSignalStrength 1 "
                                         if "r" in self._freezeParameters else "")

        # # using snapshot
        # combine -M MultiDimFit -d workspace.root --saveWorkspace -m 0 -n '.snapshot' --setParameters r=1 --toysFrequentist --cminDefaultMinimizerStrategy 2 --cminDefaultMinim# izerTolerance 0.01 --seed 42 --freezeParameters r

        command_string += exec_bash(
            "combine -M MultiDimFit   -d {WORKSPACE} --saveWorkspace -m 0 {POI} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Snapshot\" --seed {SEED}"
            .format(WORKSPACE=self.workspace,
                    FREEZEPARAMS=self.freezeParameters,
                    EXTRA=self.extraOptions,
                    NAME=self.name,
                    POI=self.POI,
                    SEED=self.seed), debug)
        # # snapshot + observed
        # combine -d snapshot.root -m 0 --snapshotName MultiDimFit --bypassFrequentistFit -M GoodnessOfFit --algo saturated --setParameters r=1 --freezeParameters r --cminDefaultMinimizerStrategy 2 --cminDefaultMinimizerTolerance 0.01 --seed 42 --fixedSignalStrength 1
        command_string += exec_bash(
            "combine -M GoodnessOfFit -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root --snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Baseline\" --algo={ALGO}"
            .format(WORKSPACE=self.workspace,
                    NAME=self.name,
                    FREEZEPARAMS=self.freezeParameters,
                    POI=self.POI,
                    SEED=self.seed,
                    ALGO=self.algo,
                    EXTRA=GOF_extra), debug)
        # # snapshot + generate-only
        # combine -M GenerateOnly -d snapshot.root -m 0 --snapshotName MultiDimFit --bypassFrequentistFit -t 1 --toysFrequentist --seed 42 --cminDefaultMinimizerStrategy 2 --cminDefaultMinimizerTolerance 0.01 --setParameters r=1 --freezeParameters r --saveToys  -n '.snapshot2gen'
        command_string += exec_bash(
            "combine -M GenerateOnly  -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root --snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} {TOYSOPTIONS}  --saveToys  "
            .format(WORKSPACE=self.workspace,
                    NAME=self.name,
                    FREEZEPARAMS=self.freezeParameters,
                    POI=self.POI,
                    SEED=self.seed,
                    NTOYS=self.toys,
                    TOYSOPTIONS=self.toysOptions,
                    EXTRA=self.extraOptions), debug)
        # combine -d snapshot.root -m 0 --snapshotName MultiDimFit --bypassFrequentistFit -M GoodnessOfFit --algo saturated -t 1 --toysFrequentist --setParameters r=1 --cminDefaultMinimizerStrategy 2 --cminDefaultMinimizerTolerance 0.01 --toysFile toys.root --seed 42 -n '.snapshot2gen' --freezeParameters r
        command_string += exec_bash(
            "combine -M GoodnessOfFit -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root --snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} {TOYSOPTIONS}  --toysFile higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --algo={ALGO}     "
            .format(WORKSPACE=self.workspace,
                    NAME=self.name,
                    FREEZEPARAMS=self.freezeParameters,
                    POI=self.POI,
                    SEED=self.seed,
                    NTOYS=self.toys,
                    TOYSOPTIONS=self.toysOptions,
                    ALGO=self.algo,
                    EXTRA=GOF_extra), debug)

        #standalone
        # command_string += exec_bash("combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Baseline\" --algo={ALGO}".format(WORKSPACE=self.workspace,NAME=self.name,FREEZEPARAMS=self.freezeParameters,POI=self.POI,SEED=self.seed,ALGO=self.algo,EXTRA=GOF_extra),debug)
        # command_string += exec_bash("combine -M GenerateOnly  -d {WORKSPACE} -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} {TOYSOPTIONS}  --saveToys  ".format(WORKSPACE=self.workspace,NAME=self.name,FREEZEPARAMS=self.freezeParameters,POI=self.POI,SEED=self.seed,NTOYS=self.toys,TOYSOPTIONS=self.toysOptions,EXTRA=self.extraOptions),debug)
        # command_string += exec_bash("combine -M GoodnessOfFit -d {WORKSPACE} -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} {TOYSOPTIONS}  --toysFile higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --algo={ALGO}     ".format(WORKSPACE=self.workspace,NAME=self.name,FREEZEPARAMS=self.freezeParameters,POI=self.POI,SEED=self.seed,NTOYS=self.toys,TOYSOPTIONS=self.toysOptions,ALGO=self.algo,EXTRA=GOF_extra),debug)

        command_string += exec_bash("", debug)

        for model_dir in alt_model_dirs:
            command_string += exec_bash("", debug)
            command_string += exec_bash("cd " + model_dir, debug)
            command_string += exec_bash("source build.sh", debug)

            #using snapshot
            command_string += exec_bash(
                "combine -M MultiDimFit   -d *_combined.root --saveWorkspace -m 0 {POI} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Snapshot\" --seed {SEED}"
                .format(WORKSPACE=self.workspace,
                        FREEZEPARAMS=self.freezeParameters,
                        EXTRA=self.extraOptions,
                        NAME=self.name,
                        POI=self.POI,
                        SEED=self.seed), debug)
            command_string += exec_bash(
                "combine -M GoodnessOfFit -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root --snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Baseline\" --algo={ALGO}"
                .format(WORKSPACE=self.workspace,
                        NAME=self.name,
                        FREEZEPARAMS=self.freezeParameters,
                        POI=self.POI,
                        SEED=self.seed,
                        ALGO=self.algo,
                        EXTRA=GOF_extra), debug)
            command_string += exec_bash(
                "combine -M GoodnessOfFit -d higgsCombine{NAME}Snapshot.MultiDimFit.mH0.{SEED}.root --snapshotName MultiDimFit --bypassFrequentistFit -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} {TOYSOPTIONS}  --toysFile {BASEMODELDIR}/higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --algo={ALGO}     "
                .format(WORKSPACE=self.workspace,
                        NAME=self.name,
                        FREEZEPARAMS=self.freezeParameters,
                        POI=self.POI,
                        SEED=self.seed,
                        NTOYS=self.toys,
                        TOYSOPTIONS=self.toysOptions,
                        ALGO=self.algo,
                        EXTRA=GOF_extra,
                        BASEMODELDIR=self.modeldir), debug)
            #standalone
            # command_string += exec_bash("combine -M GoodnessOfFit -d *_combined.root -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Baseline\" --algo={ALGO}".format(WORKSPACE=self.workspace,NAME=self.name,FREEZEPARAMS=self.freezeParameters,POI=self.POI,SEED=self.seed,ALGO=self.algo,EXTRA=GOF_extra),debug)
            # command_string += exec_bash("combine -M GoodnessOfFit -d *_combined.root -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n \"{NAME}\" -t {NTOYS} {TOYSOPTIONS}  --toysFile {BASEMODELDIR}/higgsCombine{NAME}.GenerateOnly.mH0.{SEED}.root --algo={ALGO}".format(WORKSPACE=self.workspace,NAME=self.name,FREEZEPARAMS=self.freezeParameters,POI=self.POI,SEED=self.seed,NTOYS=self.toys,TOYSOPTIONS=self.toysOptions,ALGO=self.algo,EXTRA=GOF_extra,BASEMODELDIR=self.modeldir),debug)

        command_string += exec_bash("cd " + self.modeldir + "\n", debug)

        return command_string

    def FTest(self, debug=True):
        parsed_args = self.parser[0]
        command_string = "#FTest\n"
        print("MODELDIR", self.modeldir)
        if (not debug):
            os.chdir(self.modeldir + '/bkgtests')
        command_string += exec_bash("mkdir " + self.modeldir + '/bkgtests', debug)
        if not parsed_args.collect:
            command_string += exec_bash("rm " + self.modeldir + '/bkgtests/*', debug)

        GOF_extra = self.extraOptions + (" --fixedSignalStrength 1 "
                                         if "r" in self._freezeParameters else "")

        CONDOR_str = """ --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name "base{}" """ if parsed_args.condor else ""
        if not parsed_args.collect:
            # Ref Fit (base)
            command_string += exec_bash(
                "combineTool.py -M MultiDimFit -d {WORKSPACE} --saveWorkspace "
                " -m 0 {POI} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Snapshot\" "
                .format(WORKSPACE=self.workspace,
                        FREEZEPARAMS=self.freezeParameters,
                        EXTRA=self.extraOptions,
                        NAME=self.name,
                        POI=self.POI,
                        SEED=self.seed), debug)
            command_string += exec_bash("mv higgsCombineSnapshot.MultiDimFit.mH0.root mdfit.root", debug)

            # Ref GoF (base)
            command_string += exec_bash(
                "combineTool.py -M GoodnessOfFit -d mdfit.root"
                " --snapshotName MultiDimFit --bypassFrequentistFit"
                " -m 0 {POI} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Baseline\" --algo={ALGO}"
                .format(WORKSPACE=self.workspace,
                        NAME=self.name,
                        FREEZEPARAMS=self.freezeParameters,
                        POI=self.POI,
                        SEED=self.seed,
                        ALGO=self.algo,
                        EXTRA=GOF_extra), debug)
            command_string += exec_bash("mv higgsCombineBaseline.GoodnessOfFit.mH0.root refbase.root", debug)

            # Generate toys (base) - done only once
            command_string += exec_bash(
                "combineTool.py -M GenerateOnly  -d mdfit.root"
                " --snapshotName MultiDimFit --bypassFrequentistFit "
                " -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n \"{NAME}Toys\" -t {NTOYS} {TOYSOPTIONS}  --saveToys"
                .format(WORKSPACE=self.workspace,
                        NAME=self.name,
                        FREEZEPARAMS=self.freezeParameters,
                        POI=self.POI,
                        SEED=self.seed,
                        NTOYS=self.toys,
                        TOYSOPTIONS=self.toysOptions,
                        EXTRA=self.extraOptions), debug)

        toyfiles = sorted([a for a in os.listdir('.') if re.search("higgsCombineToys\.GenerateOnly\.mH0\..+\.root", a)])
        toysids = [int(a.split("mH0.")[-1].split(".")[0]) for a in toyfiles]

        # GoFs (base)
        if not parsed_args.collect:
            for sid, toyfile in zip(toysids, toyfiles):
                command_string += exec_bash(
                    "combineTool.py -M GoodnessOfFit -d mdfit.root"
                    " --snapshotName MultiDimFit --bypassFrequentistFit"
                    " -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n \"{NAME}GoFs\" -t {NTOYS} {TOYSOPTIONS}"
                    " --toysFile {TOYFILE} --algo={ALGO}"
                    .format(WORKSPACE=self.workspace,
                            NAME=self.name,
                            FREEZEPARAMS=self.freezeParameters,
                            POI=self.POI,
                            SEED=sid,
                            TOYFILE=toyfile,
                            NTOYS=self.toys,
                            TOYSOPTIONS=self.toysOptions,
                            ALGO=self.algo,
                            EXTRA=GOF_extra)
                    + CONDOR_str.format(sid),
                    debug)

        if not parsed_args.condor or parsed_args.collect:
            command_string += exec_bash("hadd -f gofsbase.root higgsCombineGoFs.GoodnessOfFit*root", debug)
            command_string += exec_bash("rm higgsCombineGoFs.GoodnessOfFit*root", debug)

        command_string += exec_bash("", debug)

        if not parsed_args.collect:
            # Ref fit (alt)
            command_string += exec_bash(
                "combineTool.py -M MultiDimFit -d {WORKSPACE} --saveWorkspace"
                " -m 0 {POI} {FREEZEPARAMS} {EXTRA} -n \"{NAME}SnapshotAltModel\" "
                .format(WORKSPACE=self.altmodel,
                        FREEZEPARAMS=self.freezeParameters,
                        EXTRA=self.extraOptions,
                        NAME=self.name,
                        POI=self.POI,
                        SEED=self.seed), debug)
            command_string += exec_bash("mv higgsCombineSnapshotAltModel.MultiDimFit.mH0.root mdfitalt.root")

            # Ref GoF (alt)
            command_string += exec_bash(
                "combineTool.py -M GoodnessOfFit -d mdfitalt.root"
                " --snapshotName MultiDimFit --bypassFrequentistFit"
                " -m 0 {POI} {FREEZEPARAMS} {EXTRA} -n \"{NAME}BaselineAltModel\" --algo={ALGO}"
                .format(NAME=self.name,
                        FREEZEPARAMS=self.freezeParameters,
                        POI=self.POI,
                        SEED=self.seed,
                        ALGO=self.algo,
                        EXTRA=GOF_extra), debug)
            command_string += exec_bash("mv higgsCombineBaselineAltModel.GoodnessOfFit.mH0.root refalt.root")

            # GoFs (alt)
            for sid, toyfile in zip(toysids, toyfiles):
                command_string += exec_bash(
                    "combineTool.py -M GoodnessOfFit -d mdfitalt.root"
                    " --snapshotName MultiDimFit --bypassFrequentistFit"
                    " -m 0 {POI} --seed {SEED} {FREEZEPARAMS} {EXTRA} -n \"{NAME}GoFsAlt\" -t {NTOYS} {TOYSOPTIONS}"
                    " --toysFile {TOYFILE} --algo={ALGO}"
                    .format(WORKSPACE=self.workspace,
                            NAME=self.name,
                            FREEZEPARAMS=self.freezeParameters,
                            POI=self.POI,
                            SEED=sid,
                            TOYFILE=toyfile,
                            NTOYS=self.toys,
                            TOYSOPTIONS=self.toysOptions,
                            ALGO=self.algo,
                            EXTRA=GOF_extra)
                    + CONDOR_str.format(sid),
                    debug)

        if not parsed_args.condor or parsed_args.collect:                  
            command_string += exec_bash("hadd -f gofsalt.root higgsCombineGoFsAlt.GoodnessOfFit*root", debug)
            command_string += exec_bash("rm higgsCombineGoFsAlt.GoodnessOfFit*root", debug)

        command_string += exec_bash("cd " + self.modeldir + "\n", debug)

        return command_string


if (__name__ == '__main__'):
    parser = argparse.ArgumentParser()

    parser.add_argument('--justplots', action="store_true")
    parser.add_argument('--skipplots', action="store_true")
    parser.add_argument('--debug', action="store_true")
    parser.add_argument('--condor', action="store_true")
    parser.add_argument('--collect', action="store_true")
    parser.add_argument('--method',
                        type=str,
                        choices=globals()['CombineWorkflows']().methods,
                        required=True)
    parser.add_argument('--POI', default="r")
    parser.add_argument('--workspace', '-w', default='WJetsOneScale_combined.root')
    parser.add_argument('--altmodel', '-a', default=None)
    parser.add_argument('--workers', default=5)
    parser.add_argument('--freezeParameters', default=None)
    parser.add_argument('--lumi', default=41.8)
    parser.add_argument('-n', '--name', default="")
    parser.add_argument('--seed', default="1234567")
    parser.add_argument('-t', '--toys', default=50)
    parser.add_argument('--algo', default="saturated")
    parser.add_argument('--extra',
                        default="",
                        help='pass extra arguments/options to combine commands')
    parser.add_argument('--job_index', default=0, type=int)
    parser.add_argument('--externToys', action="store_true")
    parser.add_argument(
        '--rhalphdir',
        type=str,
        default=
        "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/rhalph"
    )

    args = parser.parse_args()

    if (not os.path.isfile(args.workspace)):
        raise IOError('Could not find workspace file')

    args.modeldir = os.path.abspath('/'.join(args.workspace.split('/')[:-1]) +
                                    '/') if '/' in args.workspace else ''
    if (args.job_index > 0):
        args.seed = str(int(args.seed) + args.job_index)
        print('jobindex!=0. resetting seed to initial_seed+job_index:', args.seed)

    # print('workspace',args.workspace)
    # print('model_dir',args.modeldir)
    # print('using method',args.method)
    cw = CombineWorkflows(parser=args)
    # setting up CombineWorkflow (this is written with python2 in mind. So the property decorators defined above need to be "recreated" here)
    cw.method = args.method
    cw.POI = "" if args.POI == "r" else ("--redefineSignalPOIs" + args.POI)
    cw.workspace = os.path.abspath(args.workspace)
    cw.altmodel = os.path.abspath(args.altmodel)
    cw.freezeParameters = "" if args.freezeParameters == '' else (
        '--freezeParameters ' + args.freezeParameters)
    cw.name = args.name
    cw.seed = args.seed
    cw.toys = args.toys
    cw.algo = args.algo
    cw.condor = args.condor
    cw.extra = args.extra
    cw.job_index = args.job_index
    cw.externToys = args.externToys
    cw.rhalphdir = args.rhalphdir
    cw.modeldir = args.modeldir

    method = getattr(cw, args.method)
    print(cw)
    command_string = method(args.debug)
    if (args.debug):
        print()
        print()
        print(command_string)
