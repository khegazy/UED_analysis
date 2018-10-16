#!/bin/bash 

./pairCorr.exe 20180627_1551 #-FillQ false
./pairCorr.exe 20180629_1630 #-FillQ false
./pairCorr.exe 20180630_1925 #-FillQ false
./pairCorr.exe 20180701_0746 #-FillQ false


#######################################################
#####  Difference PC with Simulated Final States  #####
#######################################################

#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/ -Fname phenoxyRadical_sMsPatternLineOut_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[555].dat -SubR true 

##############################################
#####  Reference Image / Theory Compare  #####
##############################################

#RUNNAMES=( "20180627_1551" "20180629_1630" "20180630_1925" "20180701_0746" )

#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/ -Fname nitrobenzene_sMsPatternLineOut_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[555].dat -FillQ false

#for run in "${RUNNAMES[@]}"
#do
#  ./pairCorr.exe ${run} -Idir ../staticDiffraction/results/ -Fname staticDiffraction_${run}[555].dat -Osuf -reference -FillQ theory
#done

### change molecule to phenoxyRadical
#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/ -Fname phenoxyRadical_sMsPatternLineOut_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[555].dat -FillQ false
