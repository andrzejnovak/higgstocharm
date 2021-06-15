# H->cc on top of rhalphalib

![Ralph](https://upload.wikimedia.org/wikipedia/en/thumb/1/14/Ralph_Wiggum.png/220px-Ralph_Wiggum.png)

## Install
Following the [recipe](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#cc7-release-cmssw_10_2_x-recommended-version) from combine. and clone to CMSSW environment.
```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.0.1
scramv1 b clean; scramv1 b

cd $CMSSW_BASE/src/
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
scram b -j8
```
Then install rhalphalib

```
cmsenv
cd $CMSSW_BASE/src/
git clone git@github.com:andrzejnovak/rhalphalib.git
cd rhalphalib
git fetch
git checkout origin/newhcc
# Need to update some packages against the ones in CMSSW (might need a few more)
pip install uproot --user --upgrade
pip install matplotlib --user --upgrade
pip install mplhep --user
# Install rhalphalib 
pip install --user -e .
```

and finally clone and run higgstocharm

```
git clone https://github.com/andrzejnovak/higgstocharm.git
cd higgstocharm

# Must chose --data or --MC, other options get printed
# python new_Hxx.py --data --unblind --year 2017 --templates n2nano/templates_nskim17_CC.root -o Test17
python new_Hxx.py --data --unblind --year 2017 -t tau/templates_new17_CC.root -o Test17 --degs 0,0 --fast 1
```

## Fitting

### Building workspace commands
```
bash build.sh
text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/*hcc*:r[1,-500,500]' --PO 'map=.*/zcc:z[1,-5,5]' model_combined.txt
# text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/*hcc*:r[1,-500,500]' model_combined.txt
# text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/zcc:r[1,-5,5]' model_combined.txt

combine -M FitDiagnostics --expectSignal 1 -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t -1 --toysFrequentist 
combine -M FitDiagnostics --expectSignal 1 -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 --saveShapes --saveWithUncertainties -t -1 --toysFrequentist --setParameters z=1
combine -M Significance model_combined.root --expectSignal 1 --redefineSignalPOIs z -t -1 --toysFrequentist
combineTool.py -M AsymptoticLimits -m 125 -d model_combined.root --expectSignal 1 --setParameters z=1 --redefineSignalPOIs r -t -1 --toysFrequentist 
python ../../../HiggsAnalysis/CombinedLimit/test/diffNuisances.py fitDiagnostics.root 


python ../plot.py --data 
python ../plotTF.py

# python ../plot.py --MC --year 2017 -o plots_MC_t1
```

### Running Impacts
Fitting Z
```
# Baseline
combineTool.py -M Impacts -d model_combined.root -m 125 --doInitialFit --robustFit 1 --setParameterRanges r=-1,5 --cminDefaultMinimizerStrategy 0 --X-rtd FITTER_DYN_STEP --expectSignal 1 -t -1 --toysFrequentist 
# Condor
combineTool.py -M Impacts -d model_combined.root -m 125 --doFits --robustFit 1 --allPars --setParameterRanges r=-1,5  -t -1 --toysFrequentist --expectSignal 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHccZ --exclude 'rgx{qcdparams*}'
# Collect
combineTool.py -M Impacts -d model_combined.root -o impacts.json
plotImpacts.py -i impacts.json -o plots/impacts_out_Z --blind
```

Fitting H
```
# Baseline
combineTool.py -M Impacts -d model_combined.root -m 125 --doInitialFit --robustFit 1 --setParameterRanges r=-500,500 --cminDefaultMinimizerStrategy 0 --X-rtd FITTER_DYN_STEP --expectSignal 40 -t -1 --toysFrequentist 
# Condor
combineTool.py -M Impacts -d model_combined.root -m 125 --doFits --robustFit 1 --allPars --setParameterRanges r=-500,500  -t -1 --toysFrequentist --expectSignal 40 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHccH --exclude 'rgx{qcdparams*}'
# Collect
combineTool.py -M Impacts -m 125 -d model_combined.root -o impacts.json
plotImpacts.py -i impacts.json -o plots/impacts_out_H --blind
```


### Running bias tests
Ensure signal min/max are sufficiently large
```
export BIAS=bias0
combineTool.py -M FitDiagnostics --expectSignal 0 -n $BIAS -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 20 -s 1:50:1  --toysFrequentist --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$BIAS --random
export BIAS=bias1
combineTool.py -M FitDiagnostics --expectSignal 1 -n $BIAS -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 25 -s 1:50:1  --toysFrequentist --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$BIAS --random
export BIAS=bias5
combineTool.py -M FitDiagnostics --expectSignal 5 -n $BIAS -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 25 -s 1:50:1  --toysFrequentist --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$BIAS --random
export BIAS=bias10
combineTool.py -M FitDiagnostics --expectSignal 10 -n $BIAS -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 25 -s 1:50:1  --toysFrequentist --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$BIAS --random
export BIAS=bias30
combineTool.py -M FitDiagnostics --expectSignal 30 -n $BIAS -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 25 -s 1:50:1  --toysFrequentist --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$BIAS --random
export BIAS=bias50
combineTool.py -M FitDiagnostics --expectSignal 50 -n $BIAS -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 25 -s 1:50:1  --toysFrequentist --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$BIAS --random
export BIAS=bias70
combineTool.py -M FitDiagnostics --expectSignal 70 -n $BIAS -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 25 -s 1:50:1  --toysFrequentist --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$BIAS --random
export BIAS=bias100
combineTool.py -M FitDiagnostics --expectSignal 100 -n $BIAS -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 25 -s 1:50:1  --toysFrequentist --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$BIAS --random
```
```
for BIAS in bias0 bias1 bias5 bias10 bias30 bias50 bias70 bias100
    do 
    hadd -f $BIAS.root *Combine$BIAS.*
    done

```

```
export BIAS=bias1k
combineTool.py -M FitDiagnostics --expectSignal 100 -n $BIAS -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 25 -s 1000:1020:1  --toysFrequentist --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$BIAS
export BIAS=bias5k
combineTool.py -M FitDiagnostics --expectSignal 100 -n $BIAS -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 25 -s 5000:5020:1  --toysFrequentist --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$BIAS
export BIAS=bias10k
combineTool.py -M FitDiagnostics --expectSignal 100 -n $BIAS -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 25 -s 10000:10020:1  --toysFrequentist --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$BIAS
export BIAS=bias100k
combineTool.py -M FitDiagnostics --expectSignal 100 -n $BIAS -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 25 -s 100000:100020:1  --toysFrequentist --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$BIAS
export BIAS=bias500k
combineTool.py -M FitDiagnostics --expectSignal 100 -n $BIAS -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 25 -s 500000:500020:1  --toysFrequentist --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$BIAS
```
```
for BIAS in bias1k bias5k bias10k bias100k bias500k
    do 
    hadd -f $BIAS.root *Combine$BIAS.*
    done

```


### Running likelihood scan
```
combineTool.py -M MultiDimFit -d model_combined.root --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 40 --setParameterRanges r=-20,20 -m 125 -t -1 --toysFrequentist

plot1DScan.py higgsCombine.Test.MultiDimFit.mH125.root -o plots/LScan_data_Zexp1 --y-max 30 --y-cut 30
```
