#!/bin/bash 

MAINDIR=/reg/neh/home/khegazy/analysis/nitroBenzene/UED
ROOTOUTDIR=/reg/ued/ana/scratch/nitroBenzene/rootFiles/
FILETORUN=preProcessing.exe

if [ -z "$1" ]; then
  echo "ERROR SUBMITTING JOBS!!!   Must give the type and and optionally the name of the run (i.e Run/Background and 20180701_0746)!"
  exit
fi

RUNTYPE=${1}

if [ -z "$2" ]; then
  RUNNAME="*"
else
  RUNNAME=${2}"*"
fi

if [ -z "$3" ]; then
  RESUBMIT=false
else
  RESUBMIT=true
fi

make clean; make
sleep 10
OUTPUTDIR=${ROOTOUTDIR}/logs/
for file in runLists/runList_${RUNTYPE}-${RUNNAME}
do

  strPos=$(($strPos + $subDLength + 11))

  CUTBEGINNING=${file#*runList_}
  OUTPUTFILENAME=${CUTBEGINNING%.*} #${1}-${2}
  ROOTFILENAME=${ROOTOUTDIR}${OUTPUTFILENAME}".root"

  if [[ "$RESUBMIT" = true ]] && [[ ! -f ${ROOTFILENAME} ]]; then
    continue
  fi

  echo "Submitting "${file}
  bsub -q psanaq -o${OUTPUTDIR}${OUTPUTFILENAME}".log" ./${FILETORUN} ${file}

done
