#!/bin/bash

for i in `seq 0 30`;
do
  #bsub -q psanaq -o"./output/fitOut_simulateReferenceDiff_phenoxyRadical_"${i}".txt" ./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/ -Fname phenoxyRadical_sMsPatternLineOut_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[555].dat
  bsub -q psanaq -o"./output/fitOut_20180627_1551_"${i}".txt" ./pairCorr.exe 20180627_1551 -fitTstep ${i}
  bsub -q psanaq -o"./output/fitOut_20180629_1630_"${i}".txt" ./pairCorr.exe 20180629_1630 -fitTstep ${i}
  bsub -q psanaq -o"./output/fitOut_20180630_1925_"${i}".txt" ./pairCorr.exe 20180630_1925 -fitTstep ${i}
  bsub -q psanaq -o"./output/fitOut_20180701_0746_"${i}".txt" ./pairCorr.exe 20180701_0746 -fitTstep ${i}
done
