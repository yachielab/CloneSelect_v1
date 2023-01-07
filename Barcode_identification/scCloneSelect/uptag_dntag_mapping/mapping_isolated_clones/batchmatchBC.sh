#!/bin/bash

#batchmatchBC.sh
#to run matchBC.sh for mutiple samples
#run in: library folder ~/DryLab/Projects/scCloneSelect/scCSExp17.3-3/CloneIsolation/barcodematch
#run by: bash batchmatchBC.sh

source ~/.bashrc

for i in {67..76};
do
    python matchBC.py ~/DryLab/Projects/scCloneSelect/scCSExp17.3-3/CloneIsolation/barcodeextraction/${i}-${i}/BClibsplit.csv ./scCSlibv3.csv --o matchedBC${i}.csv --distance 3
    echo $i
done