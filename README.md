G# H->cc on top of rhalphalib

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

# Unblinding 
combine -M FitDiagnostics -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1  --setParameters z=1 -n ""  --freezeParameters z 
combineTool.py -M AsymptoticLimits -m 125 -d model_combined.root --setParameters z=1  --freezeParameters z --redefineSignalPOIs r


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
plotImpacts.py -i impactsZunbl.json -o plots/impacts_out_Zunbl 
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


Fitting H unblinding
```
# Baseline
combineTool.py -M Impacts -d model_combined.root -m 125 --doInitialFit --robustFit 1 --setParameterRanges r=-200,200 --cminDefaultMinimizerStrategy 0 --X-rtd FITTER_DYN_STEP --setParameters z=1 --freezeParameters=z --redefineSignalPOIs r
# Condor
combineTool.py -M Impacts -d model_combined.root -m 125 --doFits --robustFit 1 --allPars --setParameterRanges r=-200,200   --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHccH --exclude 'rgx{qcdparams*}' --setParameters z=1 --freezeParameters=z --redefineSignalPOIs r
# Collect
combineTool.py -M Impacts -m 125 -d model_combined.root -o impactsHunbl.json --redefineSignalPOIs r
plotImpacts.py -i impactsHunbl.json -o plots/impacts_out_Hunbl
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


### Running likelihood scan
#### Z
```
bash build.sh
combineTool.py -M MultiDimFit -d model_combined.root --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 40 --setParameterRanges z=0,2 -m 125 --redefineSignalPOIs=z
plot1DScan.py higgsCombine.Test.MultiDimFit.mH125.root -o plots/LScan_data_Z_unbl --y-max 10 --y-cut 10 --POI=z

# Split year form full fit
text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/*hcc*:r[1,-500,500]' --PO 'map=.*2016/zcc:z16[1,0,2]' --PO 'map=.*2017/zcc:z17[1,0,2]' --PO 'map=.*2018/zcc:z18[1,0,2]' model_combined.txt

combineTool.py -M MultiDimFit -d model_combined.root --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 40 --setParameterRanges z16=0,2 -m 125 --redefineSignalPOIs=z16 -n .z16scan
plot1DScan.py higgsCombine.z16scan.MultiDimFit.mH125.root -o plots/LScan_data_Z_unbl16 --y-max 10 --y-cut 10 --POI=z16 

combineTool.py -M MultiDimFit -d model_combined.root --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 40 --setParameterRanges z17=0,2 -m 125 --redefineSignalPOIs=z17 -n .z17scan
plot1DScan.py higgsCombine.z17scan.MultiDimFit.mH125.root -o plots/LScan_data_Z_unbl17 --y-max 10 --y-cut 10 --POI=z17

combineTool.py -M MultiDimFit -d model_combined.root --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 40 --setParameterRanges z18=0,2 -m 125 --redefineSignalPOIs=z18 -n .z18scan
plot1DScan.py higgsCombine.z18scan.MultiDimFit.mH125.root -o plots/LScan_data_Z_unbl18 --y-max 10 --y-cut 10 --POI=z18 

```

#### Higgs
```
bash build.sh
combineTool.py -M MultiDimFit -d model_combined.root --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 40 --setParameterRanges r=-100,100 -m 125 --redefineSignalPOIs=r -n Hscan --freezeParameters z --setParameters z=1
plot1DScan.py higgsCombineHscan.MultiDimFit.mH125.root -o plots/LScan_data_H_unbl --y-max 30 --y-cut 30 --POI=r

# Split year form full fit
text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/zcc:z[1,0,2]' --PO 'map=.*2016/*hcc*:r16[1,-200,400]' --PO 'map=.*2017/*hcc*:r17[1,-300,300]' --PO 'map=.*2018/*hcc*:r18[1,-300,300]' model_combined.txt

combineTool.py -M MultiDimFit -d model_combined.root --cminDefaultMinimizerStrategy 0  --robustFit 1 --algo grid --points 40 -m 125 --redefineSignalPOIs=r16 -n .h16scan --freezeParameters z --setParameters z=1
plot1DScan.py higgsCombine.h16scan.MultiDimFit.mH125.root -o plots/LScan_data_H_unbl16 --y-max 30 --y-cut 30 --POI=r16 

combineTool.py -M MultiDimFit -d model_combined.root --cminDefaultMinimizerStrategy 0  --robustFit 1 --algo grid --points 40 -m 125 --redefineSignalPOIs=r17 -n .h17scan --freezeParameters z --setParameters z=1
plot1DScan.py higgsCombine.h17scan.MultiDimFit.mH125.root -o plots/LScan_data_H_unbl17 --y-max 30 --y-cut 30 --POI=r17

combineTool.py -M MultiDimFit -d model_combined.root --cminDefaultMinimizerStrategy 0  --robustFit 1 --algo grid --points 40 -m 125 --redefineSignalPOIs=r18 -n .h18scan --freezeParameters z --setParameters z=1
plot1DScan.py higgsCombine.h18scan.MultiDimFit.mH125.root -o plots/LScan_data_H_unbl18 --y-max 30 --y-cut 30 --POI=r18 

```

### Channel (year) compatibility 
#### Z
```bash
bash build.sh
combineTool.py -M MultiDimFit -d model_combined.root --algo singles --cminDefaultMinimizerStrategy 0 --setParameters r=1,z=1 --redefineSignalPOIs z -n .cmb

# Split year form full fit
text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/*hcc*:r[1,-500,500]' --PO 'map=.*2016/zcc:z16[1,0,2]' --PO 'map=.*2017/zcc:z17[1,0,2]' --PO 'map=.*2018/zcc:z18[1,0,2]' model_combined.txt

combineTool.py -M MultiDimFit -d model_combined.root --algo singles --cminDefaultMinimizerStrategy 0 --setParameters r=1,z16=1,z17=1,z18=1 --redefineSignalPOIs z16 -n .z16
combineTool.py -M MultiDimFit -d model_combined.root --algo singles --cminDefaultMinimizerStrategy 0 --setParameters r=1,z16=1,z17=1,z18=1 --redefineSignalPOIs z17 -n .z17
combineTool.py -M MultiDimFit -d model_combined.root --algo singles --cminDefaultMinimizerStrategy 0 --setParameters r=1,z16=1,z17=1,z18=1 --redefineSignalPOIs z18 -n .z18

combineTool.py -M PrintFit --json MultiDimFit_ccc.json -P z -i higgsCombine.cmb.MultiDimFit.mH120.root --algo singles
combineTool.py -M PrintFit --json MultiDimFit_ccc.json -P z16 -i higgsCombine.z16.MultiDimFit.mH120.root --algo singles
combineTool.py -M PrintFit --json MultiDimFit_ccc.json -P z17 -i higgsCombine.z17.MultiDimFit.mH120.root --algo singles
combineTool.py -M PrintFit --json MultiDimFit_ccc.json -P z18 -i higgsCombine.z18.MultiDimFit.mH120.root --algo singles

python ../plotCCC.py -i MultiDimFit_ccc.json -o plots/ccc_corr -z
```

#### H
```bash
bash build.sh
combineTool.py -M MultiDimFit -d model_combined.root --algo singles --cminDefaultMinimizerStrategy 0 --setParameters r=1,z=1 --redefineSignalPOIs r -n .cmbHiggs --freezeParameters z

# Split year form full fit
text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/zcc:z[1,0,2]' --PO 'map=.*2016/*hcc*:r16[1,-200,300]' --PO 'map=.*2017/*hcc*:r17[1,-200,200]' --PO 'map=.*2018/*hcc*:r18[1,-200,200]' model_combined.txt

combineTool.py -M MultiDimFit -d model_combined.root --algo singles --cminDefaultMinimizerStrategy 0 --setParameters z=1,r16=1,r17=1,r18=1 --redefineSignalPOIs r16 -n .r16 --freezeParameters z 
combineTool.py -M MultiDimFit -d model_combined.root --algo singles --cminDefaultMinimizerStrategy 0 --setParameters z=1,r16=1,r17=1,r18=1 --redefineSignalPOIs r17 -n .r17 --freezeParameters z 
combineTool.py -M MultiDimFit -d model_combined.root --algo singles --cminDefaultMinimizerStrategy 0 --setParameters z=1,r16=1,r17=1,r18=1 --redefineSignalPOIs r18 -n .r18 --freezeParameters z 

combineTool.py -M PrintFit --json MultiDimFit_ccc_higgs.json -P r -i higgsCombine.cmbHiggs.MultiDimFit.mH120.root --algo singles
combineTool.py -M PrintFit --json MultiDimFit_ccc_higgs.json -P r16 -i higgsCombine.r16.MultiDimFit.mH120.root --algo singles
combineTool.py -M PrintFit --json MultiDimFit_ccc_higgs.json -P r17 -i higgsCombine.r17.MultiDimFit.mH120.root --algo singles
combineTool.py -M PrintFit --json MultiDimFit_ccc_higgs.json -P r18 -i higgsCombine.r18.MultiDimFit.mH120.root --algo singles

python ../plotCCC.py -i MultiDimFit_ccc_higgs.json -o plots/ccc_corr_Higgs
```

```bash
bash build.sh
combine -M ChannelCompatibilityCheck -d model_combined.root --setParameters z=1 --redefineSignalPOIs=z  -g 2016 -g 2017 -g 2018

combine -M ChannelCompatibilityCheck -d model_combined.root --setParameters z=1 --freezeParameters z --redefineSignalPOIs=r  -g 2016 -g 2017 -g 2018
```


### Unblind limi

```bash
bash build.sh
combineTool.py -M AsymptoticLimits -m 125 -d model_combined.root --setParameters z=1 --freezeParameters z --redefineSignalPOIs r -n .cmb

text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/zcc:z[1,0,2]' --PO 'map=.*2016/*hcc*:r16[1,-200,400]' --PO 'map=.*2017/*hcc*:r17[1,-200,200]' --PO 'map=.*2018/*hcc*:r18[1,-200,200]' model_combined.txt

combineTool.py -M AsymptoticLimits -m 125 -d model_combined.root --freezeParameters z,r17,r18 --setParameters z=1,r16=1,r17=1,r18=1 --redefineSignalPOIs r16 -n .r16 &
combineTool.py -M AsymptoticLimits -m 125 -d model_combined.root --freezeParameters z,r16,r18 --setParameters z=1,r16=1,r17=1,r18=1 --redefineSignalPOIs r17 -n .r17 &
combineTool.py -M AsymptoticLimits -m 125 -d model_combined.root --freezeParameters z,r16,r17 --setParameters z=1,r16=1,r17=1,r18=1 --redefineSignalPOIs r18 -n .r18 &
```


### Uncertainties split stat+syst
```
combineTool.py -M MultiDimFit -d model_combined.root  --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 10 --setParameterRanges r=-50,100 -m 125 -n .nominal.expSig --redefineSignalPOIs=r --freezeParameters z
combineTool.py -M MultiDimFit -d model_combined.root  --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo none --setParameterRanges r=-50,100  -m 125 -n .bestfit.expSig --saveWorkspace  --redefineSignalPOIs=r --freezeParameters z

combineTool.py -M MultiDimFit -d higgsCombine.bestfit.expSig.MultiDimFit.mH125.root --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 10 --setParameterRanges r=-50,100 -m 125 -n .statOnly.expSig --snapshotName MultiDimFit --freezeParameters allConstrainedNuisances,z --redefineSignalPOIs=r
plot1DScan.py higgsCombine.nominal.expSig.MultiDimFit.mH125.root  --others 'higgsCombine.statOnly.expSig.MultiDimFit.mH125.root:Freeze all:2' --breakdown syst,stat --output plots/breakdown_run2_higgs
```

```
combineTool.py -M MultiDimFit -d model_combined.root  --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 20  --setParameterRanges z=0.5,2 -m 125 -n .nominalZ.expSig --redefineSignalPOIs=z --freezeParameters r
combineTool.py -M MultiDimFit -d model_combined.root  --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo none --setParameterRanges z=0.5,2  -m 125 -n .bestfitZ.expSig --saveWorkspace  --redefineSignalPOIs=z --freezeParameters r

combineTool.py -M MultiDimFit -d higgsCombine.bestfitZ.expSig.MultiDimFit.mH125.root --cminDefaultMinimizerStrategy 0 --expectSignal 1 --robustFit 1 --algo grid --points 20 --setParameterRanges z=0.5,2 -m 125 -n .statOnlyZ.expSig --snapshotName MultiDimFit --freezeParameters allConstrainedNuisances,r  --redefineSignalPOIs=z
plot1DScan.py higgsCombine.nominalZ.expSig.MultiDimFit.mH125.root  --others 'higgsCombine.statOnlyZ.expSig.MultiDimFit.mH125.root:Freeze all:2' --breakdown syst,stat --POI z  --output plots/breakdown_run2_Z
```

#### Split with theory
Z
```
combine -M MultiDimFit model_combined.root --setParameters r=1,z=1 --cminDefaultMinimizerStrategy 0 --redefineSignalPOIs z --setParameterRanges z=0,2 --saveWorkspace -n .postfitZ
###
combine -M MultiDimFit higgsCombine.postfitZ.MultiDimFit.mH120.root -n .Ztotal --algo grid --snapshotName MultiDimFit --setParameterRanges z=0,2 --setParameters r=1,z=1  --redefineSignalPOIs z  &
combine -M MultiDimFit higgsCombine.postfitZ.MultiDimFit.mH120.root -n .Znoexpsyst --algo grid --snapshotName MultiDimFit --setParameterRanges z=0,2 --setParameters r=1,z=1  --redefineSignalPOIs z  --freezeParameters 'rgx{(?!.*_EW$)(?!.*_NLO$)(?!.*_th_scale.pt$).*}' &
combine -M MultiDimFit higgsCombine.postfitZ.MultiDimFit.mH120.root -n .Zfreezeall --algo grid --snapshotName MultiDimFit --setParameterRanges z=0,2 --setParameters r=1,z=1  --redefineSignalPOIs z  --freezeParameters allConstrainedNuisances &
### Replace X year
plot1DScan.py higgsCombine.Ztotal.MultiDimFit.mH120.root --POI z --main-label "Total Uncert." --others higgsCombine.Znoexpsyst.MultiDimFit.mH120.root:"Theory Uncert.":4 higgsCombine.Zfreezeall.MultiDimFit.mH120.root:"Statistical Only":2 --output plots/breakdown_Z_2017 --y-max 22 --y-cut 10 --breakdown "exp,thy,stat"
```
Higgs 2016 - shifted range
```
combine -M MultiDimFit model_combined.root --setParameters r=1,z=1 --cminDefaultMinimizerStrategy 0 --redefineSignalPOIs r --setParameterRanges r=-100,300 --saveWorkspace -n postfit  --freezeParameters=z
###
combine -M MultiDimFit higgsCombinepostfit.MultiDimFit.mH120.root -n total --algo grid --snapshotName MultiDimFit --setParameterRanges r=-100,300 --setParameters r=1,z=1  --freezeParameters=z &
combine -M MultiDimFit higgsCombinepostfit.MultiDimFit.mH120.root -n noexpsyst --algo grid --snapshotName MultiDimFit --setParameterRanges r=-100,300 --setParameters r=1,z=1 --freezeParameters 'rgx{(?!.*_EW$)(?!.*_NLO$)(?!.*_th_scale.pt$)(?!z$).*}' &
combine -M MultiDimFit higgsCombinepostfit.MultiDimFit.mH120.root -n freezeall --algo grid --snapshotName MultiDimFit --setParameterRanges r=-100,300 --setParameters r=1,z=1  --freezeParameters=allConstrainedNuisances,z &
### 
plot1DScan.py higgsCombinetotal.MultiDimFit.mH120.root --POI r --main-label "Total Uncert." --others higgsCombinenoexpsyst.MultiDimFit.mH120.root:"Theory Uncert.":4 higgsCombinefreezeall.MultiDimFit.mH120.root:"Statistical Only":2 --output plots/breakdown_H_2016 --y-max 22 --y-cut 10 --breakdown "exp,thy,stat"
```

Higgs
```
combine -M MultiDimFit model_combined.root --setParameters r=1,z=1 --cminDefaultMinimizerStrategy 0 --redefineSignalPOIs r --setParameterRanges r=-100,100 --saveWorkspace -n .postfit  --freezeParameters=z
###
combine -M MultiDimFit higgsCombine.postfit.MultiDimFit.mH120.root -n .total --algo grid --snapshotName MultiDimFit --setParameterRanges r=-100,100 --setParameters r=1,z=1  --freezeParameters=z &
combine -M MultiDimFit higgsCombine.postfit.MultiDimFit.mH120.root -n .noexpsyst --algo grid --snapshotName MultiDimFit --setParameterRanges r=-100,100 --setParameters r=1,z=1 --freezeParameters 'rgx{(?!.*_EW$)(?!.*_NLO$)(?!.*_th_scale.pt$)(?!z$).*}' &
combine -M MultiDimFit higgsCombine.postfit.MultiDimFit.mH120.root -n .freezeall --algo grid --snapshotName MultiDimFit --setParameterRanges r=-100,100 --setParameters r=1,z=1  --freezeParameters=allConstrainedNuisances,z &
### Replace X year
plot1DScan.py higgsCombine.total.MultiDimFit.mH120.root --POI r --main-label "Total Uncert." --others higgsCombine.noexpsyst.MultiDimFit.mH120.root:"Theory Uncert.":4 higgsCombine.freezeall.MultiDimFit.mH120.root:"Statistical Only":2 --output plots/breakdown_H_201X --y-max 22 --y-cut 10 --breakdown "exp,thy,stat"
```

### Running GoFs
```bash
# Like do_gosf.sh <step> <algo> ?<n-jobs>
bash ../do_gof.sh 1 0 
bash ../do_gof.sh 2 0 100
bash ../do_gof.sh 3 0 
bash ../do_gof.sh 4 0 
```

### Running F-Tests
```bash
python submit_ftests.py -t 50 -s 1:10:1 --year 2017 -d FTests/data17TF -o FTests/Fouts_data17TF/ --outplots Ftests/plots_data17/ --data --make -p
python submit_ftests.py -t 50 -s 1:10:1 --year 2017 -d FTests/data17TF -o FTests/Fouts_data17TF/ --outplots Ftests/plots_data17/ --data --build
python submit_ftests.py -t 50 -s 1:10:1 --year 2017 -d FTests/data17TF -o FTests/Fouts_data17TF/ --outplots Ftests/plots_data17/ --data --run --base -p
python submit_ftests.py -t 50 -s 1:10:1 --year 2017 -d FTests/data17TF -o FTests/Fouts_data17TF/ --outplots Ftests/plots_data17/ --data --run --gen --condor -p
python submit_ftests.py -t 50 -s 1:10:1 --year 2017 -d FTests/data17TF -o FTests/Fouts_data17TF/ --outplots Ftests/plots_data17/ --data --run --fits --condor -p
```

### Check 

```
text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/zcc:r_z[1,0,2]'  --PO 'map=.*/*h*:r[1,-100,100]'  model_combined.txt
```

## Make W fit

```
text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/zcc:r_z[1,0,2]' --PO 'map=.*/wcq:r_w[1,0,2]' --PO 'map=.*/*hcc*:r[1,-200,200]'  model_combined.txt

text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/zcc:r_z[1,0,2]' --PO 'map=.*/wcq:r_w[1,0,2]'  

PTS=20
combineTool.py -M MultiDimFit -m 125 model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 --redefineSignalPOIs r_z,r_w --setParameters r_w=1,r_z=1 --setParameterRanges r_w=0,3:r_z=0,3 --algo contour2d --points=$PTS --cl=0.68 -n .68
combineTool.py -M MultiDimFit -m 125 model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 --redefineSignalPOIs r_z,r_w --setParameters r_w=1,r_z=1 --setParameterRanges r_w=0,3:r_z=0,3 --algo contour2d --points=$PTS --cl=0.95 -n .95

combineTool.py -M MultiDimFit -m 125 model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1  --redefineSignalPOIs r_z,r_w --setParameters r_w=1,r_z=1 --setParameterRanges r_w=0,3:r_z=0,2 --algo grid --points 100 -n .grid


```

# zbb fit
```
text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/zcc:z_cc[1,0,2]' --PO 'map=.*/zbb:z_bb[1,0,2]' model_combined.txt

PTS=20
combineTool.py -M MultiDimFit -m 125 model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 --redefineSignalPOIs z_bb,z_cc --setParameters z_bb=1,z_cc=1 --setParameterRanges z_bb=0,3:z_cc=0,3 --algo contour2d --points=$PTS --cl=0.68 -n .68
combineTool.py -M MultiDimFit -m 125 model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 --redefineSignalPOIs z_bb,z_cc --setParameters z_bb=1,z_cc=1 --setParameterRanges z_bb=0,3:z_cc=0,3 --algo contour2d --points=$PTS --cl=0.95 -n .95

combineTool.py -M MultiDimFit -m 125 model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1  --redefineSignalPOIs z_bb,z_cc --setParameters z_bb=1,z_cc=1 --setParameterRanges z_bb=0,3:z_cc=0,2 --algo grid --points 100 -n .grid
```


# hcc/zcc 2D scan
```
text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO verbose --PO 'map=.*/zcc:r_z[1,0,2]' --PO 'map=.*/hcc:r_h[1,-40,40]' model_combined.txt

PTS=20
combineTool.py -M MultiDimFit -m 125 model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 --redefineSignalPOIs r_z,r_h --setParameters r_z=1,r_h=1 --setParameterRanges r_z=0,3:r_h=-40,40 --algo contour2d --points=$PTS --cl=0.68 -n .68
combineTool.py -M MultiDimFit -m 125 model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 --redefineSignalPOIs r_z,r_h --setParameters r_z=1,r_h=1 --setParameterRanges r_z=0,3:r_h=-40,40 --algo contour2d --points=$PTS --cl=0.95 -n .95

combineTool.py -M MultiDimFit -m 125 model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1  --redefineSignalPOIs r_z,r_h --setParameters r_z=1,r_h=1 --setParameterRanges r_z=0,3:r_h=-40,40 --algo grid --points 100 -n .grid
```



# Checks

### Vary DDX bb SFs
```
# 0.1
python new_Hxx.py --data --unblind --year 2016 -t temps/templates_preapproval16_CC.root --mut temps/templatesmuCR_preapproval16_CC.root -o TestSF16bb01 --degs 1,0  --degsMC 0,2  --bbsf 0.1 &
python new_Hxx.py --data --unblind --year 2017 -t temps/templates_preapproval17_CC.root --mut temps/templatesmuCR_preapproval17_CC.root -o TestSF17bb01 --degs 0,0  --degsMC 1,2  --bbsf 0.1 &
python new_Hxx.py --data --unblind --year 2018 -t temps/templates_preapproval18_CC.root --mut temps/templatesmuCR_preapproval18_CC.root -o TestSF18bb01 --degs 1,0  --degsMC 0,2  --bbsf 0.1 &
# 0.3
python new_Hxx.py --data --unblind --year 2016 -t temps/templates_preapproval16_CC.root --mut temps/templatesmuCR_preapproval16_CC.root -o TestSF16bb03 --degs 1,0  --degsMC 0,2  --bbsf 0.3 &
python new_Hxx.py --data --unblind --year 2017 -t temps/templates_preapproval17_CC.root --mut temps/templatesmuCR_preapproval17_CC.root -o TestSF17bb03 --degs 0,0  --degsMC 1,2  --bbsf 0.3 &
python new_Hxx.py --data --unblind --year 2018 -t temps/templates_preapproval18_CC.root --mut temps/templatesmuCR_preapproval18_CC.root -o TestSF18bb03 --degs 1,0  --degsMC 0,2  --bbsf 0.3 &
# 0.5
python new_Hxx.py --data --unblind --year 2016 -t temps/templates_preapproval16_CC.root --mut temps/templatesmuCR_preapproval16_CC.root -o TestSF16bb05 --degs 1,0  --degsMC 0,2  --bbsf 0.5 &
python new_Hxx.py --data --unblind --year 2017 -t temps/templates_preapproval17_CC.root --mut temps/templatesmuCR_preapproval17_CC.root -o TestSF17bb05 --degs 0,0  --degsMC 1,2  --bbsf 0.5 &
python new_Hxx.py --data --unblind --year 2018 -t temps/templates_preapproval18_CC.root --mut temps/templatesmuCR_preapproval18_CC.root -o TestSF18bb05 --degs 1,0  --degsMC 0,2  --bbsf 0.5 &

mkdir TestRun2bbSF01
cp TestSF16bb01/* TestRun2bbSF01
cp TestSF17bb01/* TestRun2bbSF01
cp TestSF18bb01/* TestRun2bbSF01

mkdir TestRun2bbSF03
cp TestSF16bb03/* TestRun2bbSF03
cp TestSF17bb03/* TestRun2bbSF03
cp TestSF18bb03/* TestRun2bbSF03

mkdir TestRun2bbSF05
cp TestSF16bb05/* TestRun2bbSF05/
cp TestSF17bb05/* TestRun2bbSF05/
cp TestSF18bb05/* TestRun2bbSF05/
```

### Scale bb

```
python new_Hxx.py --data --unblind --year 2016 -t temps/templates_preapproval16_CC.root --mut temps/templatesmuCR_preapproval16_CC.root -o TestScalebb16 --degs 1,0  --degsMC 0,2  --scalebb 4 &
python new_Hxx.py --data --unblind --year 2017 -t temps/templates_preapproval17_CC.root --mut temps/templatesmuCR_preapproval17_CC.root -o TestScalebb17 --degs 0,0  --degsMC 1,2  --scalebb 4 &
python new_Hxx.py --data --unblind --year 2018 -t temps/templates_preapproval18_CC.root --mut temps/templatesmuCR_preapproval18_CC.root -o TestScalebb18 --degs 1,0  --degsMC 0,2  --scalebb 4 &

mkdir TestScalebbRun2
cp TestScalebb16/* TestScalebbRun2/
cp TestScalebb17/* TestScalebbRun2/
cp TestScalebb18/* TestScalebbRun2/
```

### Extra unc bb

```
python new_Hxx.py --data --unblind --year 2016 -t temps/templates_preapproval16_CC.root --mut temps/templatesmuCR_preapproval16_CC.root -o TestExtrabb16 --degs 1,0  --degsMC 0,2  --extrauncbb True &
python new_Hxx.py --data --unblind --year 2017 -t temps/templates_preapproval17_CC.root --mut temps/templatesmuCR_preapproval17_CC.root -o TestExtrabb17 --degs 0,0  --degsMC 1,2   --extrauncbb True &
python new_Hxx.py --data --unblind --year 2018 -t temps/templates_preapproval18_CC.root --mut temps/templatesmuCR_preapproval18_CC.root -o TestExtrabb18 --degs 1,0  --degsMC 0,2   --extrauncbb True &

mkdir TestExtrabbRun2
cp TestExtrabb16/* TestExtrabbRun2/
cp TestExtrabb17/* TestExtrabbRun2/
cp TestExtrabb17/* TestExtrabbRun2/
```