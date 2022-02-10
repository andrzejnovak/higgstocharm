#!/bin/bash

if  [[ "$2" -eq 0 ]]; then
    ALGO=saturated
elif [[ "$2" -eq 1 ]]; then
    ALGO=KS
fi
echo "Using algo:" $ALGO

if [[ "$1" -eq 0 ]]; then
    combine -M FitDiagnostics -d model_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1  --setParameters z=1,r=1 -n "" --saveShapes --saveWithUncertainties #-t -1 --toysFrequentist --robustHesse 1
elif [[ "$1" -eq 1 ]]; then
    combine -M GoodnessOfFit model_combined.root --algo $ALGO  --expectSignal 1 --redefineSignalPOIs z -n DataGoF$ALGO
elif [[ "$1" -eq 2 ]]; then
    combineTool.py -M GoodnessOfFit model_combined.root --algo $ALGO  --expectSignal 1 -t 5 --toysFrequentist --redefineSignalPOIs z -n GoFs$ALGO --job-mode condor --sub-opts='+JobFlavour = "workday"' --task-name ggHcc$ALGO -s 1:"$3":1    
elif [[ "$1" -eq 3 ]]; then 
    hadd -f allgofs$ALGO.root higgsCombineGoFs$ALGO.GoodnessOfFit.mH120.*.root
elif [[ "$1" -eq 4 ]]; then 
    python ../plot_single_gof.py higgsCombineDataGoF$ALGO.GoodnessOfFit.mH120.root allgofs$ALGO.root --algo $ALGO --year "$3"
fi


# combineTool.py -M GoodnessOfFit --algorithm $ALGO -m 125 --there -d model_combined.root -n ".$ALGO.toys" --expectSignal 1 -t 500 --toysFrequentist 
# combineTool.py -M GoodnessOfFit --algorithm $ALGO -m 125 --there -d model_combined.root -n ".$ALGO" 
# combineTool.py -M CollectGoodnessOfFit --input higgsCombine.$ALGO.GoodnessOfFit.mH125.root higgsCombine.$ALGO.toys.GoodnessOfFit.mH125.*.root -o cmb_$ALGO.json   
# python ../../../CombineHarvester/CombineTools/scripts/plotGof.py --statistic $ALGO --mass 125.0 cmb_$ALGO.json --title-right="35.9 fb^{-1} (13 TeV)" --output='-$ALGO'
