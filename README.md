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
python new_Hxx.py --data --unblind --year 2016 -t temps/templates_corr16_CC.root --mut temps/templatesmuCR_corr16_CC.root -o Correct16 --degs 1,0

python new_Hxx.py --data --unblind --year 2016 -t temps/templates_corr16_CC.root --mut temps/templatesmuCR_corr16_CC.root -o Unblind16mod --degs 1,0
```

## Fitting

### Building workspace commands
```
bash build.sh
# text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/*hcc*:r[1,-500,500]' --PO 'map=.*/zcc:z[1,-5,5]' model_combined.txt
# text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/*hcc*:r[1,-500,500]' model_combined.txt
# text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/zcc:r[1,-5,5]' model_combined.txt

combine -M FitDiagnostics -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1  --setParameters z=1,r=1 -n "" -t -1 --toysFrequentist 
combine -M FitDiagnostics --expectSignal 1 -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 --saveShapes --saveWithUncertainties -n "" -t -1 --toysFrequentist --setParameters z=1 
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
combineTool.py -M Impacts -d model_combined.root -m 125 --doInitialFit --robustFit 1 --setParameterRanges r=-1,5 --cminDefaultMinimizerStrategy 0 --X-rtd FITTER_DYN_STEP --expectSignal 1 -t -1 --toysFrequentist --redefineSignalPOIs z
# Condor
combineTool.py -M Impacts -d model_combined.root -m 125 --doFits --robustFit 1 --allPars --setParameterRanges r=-1,5  -t -1 --toysFrequentist --expectSignal 1 --redefineSignalPOIs z --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHccZ --exclude 'rgx{qcdparams*}'
# Collect
combineTool.py -M Impacts -d model_combined.root -m 125 --redefineSignalPOIs z -o impactsZ.json
plotImpacts.py -i impactsZ.json -o plots/impacts_out_Z
```

Fitting Z unblinding
```
# Baseline
combineTool.py -M Impacts -d model_combined.root -m 125 --doInitialFit --robustFit 1 --setParameterRanges r=-100,100 --cminDefaultMinimizerStrategy 0 --X-rtd FITTER_DYN_STEP --expectSignal 1 --redefineSignalPOIs z
# Condor
combineTool.py -M Impacts -d model_combined.root -m 125 --doFits --robustFit 1 --allPars --setParameterRanges r=-100,100  --redefineSignalPOIs z --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHccZ --exclude 'rgx{qcdparams*}'
# Collect
combineTool.py -M Impacts -d model_combined.root -m 125 --redefineSignalPOIs z -o impactsZunbl.json
plotImpacts.py -i impactsZunbl.json -o plots/impacts_out_Zunbl --blind
```

Fitting H
```
# Baseline
combineTool.py -M Impacts -d model_combined.root -m 125 --doInitialFit --robustFit 1 --setParameterRanges r=-500,500 --cminDefaultMinimizerStrategy 0 --X-rtd FITTER_DYN_STEP --expectSignal 40 -t -1 --toysFrequentist 
# Condor
combineTool.py -M Impacts -d model_combined.root -m 125 --doFits --robustFit 1 --allPars --setParameterRanges r=-500,500  -t -1 --toysFrequentist --expectSignal 40 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHccH --exclude 'rgx{qcdparams*}'
# Collect
combineTool.py -M Impacts -m 125 -d model_combined.root -o impactsH.json
plotImpacts.py -i impactsH.json -o plots/impacts_out_H --blind
```


### Running bias tests
Ensure signal min/max are sufficiently large

```
for bias in 0 1 `seq 5 5 100`
    do
    combineTool.py -M FitDiagnostics --expectSignal $bias -n bias$bias -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 20 -s 1:50:1 --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$bias
    done
```

```
for bias in 0 1 `seq 5 5 100`
    do 
    hadd -f bias$bias.root *Combinebias$bias.*
    done

```

```
for bias in 0 1 `seq 5 5 100`
    do
    combineTool.py -M FitDiagnostics --expectSignal $bias -n biastfFr$bias -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 20 -s 1:50:1 --toysFrequentist --freezeParameters allConstrainedNuisances --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$bias
    done
```

```
for bias in 0 1 `seq 5 5 100`
    do 
    hadd -f biastfFr$bias.root *CombinebiastfFr$bias.*
    done

```

```
for bias in 0 1 `seq 5 5 100`
    do
    combineTool.py -M FitDiagnostics --expectSignal $bias -n biasset$bias -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 -t 20 -s 1:50:1  --setParameters tf2017_MCtempl_deco0=1.62,tf2017_MCtempl_deco1=-5.34e-02,tf2017_MCtempl_deco2=-1.53e-01,tf2017_MCtempl_deco3=1.26e-01,tf2017_MCtempl_deco4=4.54e-01,tf2017_MCtempl_deco5=-1.92e-01,tf2017_dataResidual_pt_par0_rho_par0=1.0227e+00 --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$bias
    done
```

```
for bias in 0 1 `seq 5 5 100`
    do 
    hadd -f biasset$bias.root *Combinebiasset$bias.*
    done

```

tf2017_MCtempl_deco0=1.62,tf2017_MCtempl_deco1=-5.34e-02,tf2017_MCtempl_deco2=-1.53e-01,tf2017_MCtempl_deco3=1.26e-01,tf2017_MCtempl_deco4=4.54e-01,tf2017_MCtempl_deco5=-1.92e-01,tf2017_dataResidual_pt_par0_rho_par0=1.0227e+00


### Running likelihood scan
```
combineTool.py -M MultiDimFit -d model_combined.root --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 40 --setParameterRanges z=0.2,1.8 -m 125 --redefineSignalPOIs=z
plot1DScan.py higgsCombine.Test.MultiDimFit.mH125.root -o plots/LScan_data_Z_unbl --y-max 10 --y-cut 10 --POI=z

# Split year form full fit
text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/*hcc*:r[1,-500,500]' --PO 'map=.*2016/zcc:z16[1,0,2]' --PO 'map=.*2017/zcc:z17[1,0,2]' --PO 'map=.*2018/zcc:z18[1,0,2]' model_combined.txt

combineTool.py -M MultiDimFit -d model_combined.root --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 40 --setParameterRanges z16=0.2,1.8 -m 125 --redefineSignalPOIs=z16 -n .z16scan
plot1DScan.py higgsCombine.z16scan.MultiDimFit.mH125.root -o plots/LScan_data_Z_unbl16 --y-max 10 --y-cut 10 --POI=z16 

combineTool.py -M MultiDimFit -d model_combined.root --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 40 --setParameterRanges z16=0.2,1.8 -m 125 --redefineSignalPOIs=z16 -n .z16scan
plot1DScan.py higgsCombine.z16scan.MultiDimFit.mH125.root -o plots/LScan_data_Z_unbl16 --y-max 10 --y-cut 10 --POI=z16 

combineTool.py -M MultiDimFit -d model_combined.root --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 40 --setParameterRanges z16=0.2,1.8 -m 125 --redefineSignalPOIs=z16 -n .z16scan
plot1DScan.py higgsCombine.z16scan.MultiDimFit.mH125.root -o plots/LScan_data_Z_unbl16 --y-max 10 --y-cut 10 --POI=z16 

```

### Channel (year) compatibility 

```bash
combineTool.py -M MultiDimFit -d model_combined.root --algo singles --cminDefaultMinimizerStrategy 0 --setParameters r=1,z=1 --redefineSignalPOIs z -n .cmb

# Split year form full fit
text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/*hcc*:r[1,-500,500]' --PO 'map=.*2016/zcc:z16[1,0,2]' --PO 'map=.*2017/zcc:z17[1,0,2]' --PO 'map=.*2018/zcc:z18[1,0,2]' model_combined.txt

combineTool.py -M MultiDimFit -d model_combined.root --algo singles --cminDefaultMinimizerStrategy 0 --setParameters r=1,z16=1,z17=1,z18=1 --redefineSignalPOIs z16 -n .z16
combineTool.py -M MultiDimFit -d model_combined.root --algo singles --cminDefaultMinimizerStrategy 0 --setParameters r=1,z16=1,z17=1,z18=1 --redefineSignalPOIs z17 -n .z17
combineTool.py -M MultiDimFit -d model_combined.root --algo singles --cminDefaultMinimizerStrategy 0 --setParameters r=1,z16=1,z17=1,z18=1 --redefineSignalPOIs z18 -n .z18
```

```bash
combineTool.py -M PrintFit --json MultiDimFit_ccc.json -P z -i higgsCombine.cmb.MultiDimFit.mH120.root --algo singles

combineTool.py -M PrintFit --json MultiDimFit_ccc.json -P z16 -i higgsCombine.z16.MultiDimFit.mH120.root --algo singles
combineTool.py -M PrintFit --json MultiDimFit_ccc.json -P z17 -i higgsCombine.z17.MultiDimFit.mH120.root --algo singles
combineTool.py -M PrintFit --json MultiDimFit_ccc.json -P z18 -i higgsCombine.z18.MultiDimFit.mH120.root --algo singles

python ../plotCCC.py -i MultiDimFit_ccc.json -o plots/ccc_corr
```

```bash
bash build.sh
combine -M ChannelCompatibilityCheck -d model_combined.root --setParameters z=1 --redefineSignalPOIs=z  -g 2016 -g 2017 -g 2018
```
