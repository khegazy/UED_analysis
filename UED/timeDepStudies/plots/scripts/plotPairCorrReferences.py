import numpy as np
import sys
sys.path.append('../../../plots/scripts')
from plotClass import plotCLASS
from plotParams import pltParams

params = pltParams()

plc = plotCLASS()


runs = ["20180629_1630", "20180627_1551", "20180630_1925", "20180701_0746"]
NradBins = "211"

files = ["../../results/sim-simulateReference-pairCorrOdd[" + NradBins + "].dat", 
         "/reg/ued/ana/scratch/nitroBenzene/simulations/references/nitrobenzene_pairCorr_Bins[500].dat",
         ""]

opts = {
  "xTitle"     : r"R [$\AA$]",
  "yLim"       : [-0.1, 1.1],
  "labels"     : ["Ifft simulated sMs", "Charge weighted distances", "Ifft reference data sMs"]
}


#######################################
#####  Plotting pair correlation  #####
#######################################

for i,run in enumerate(runs):
  files[2] = "../../results/data-" + run + "-reference-pairCorrOdd[" + NradBins + "].dat"


  plc.print1d(files,
        "../data-" + run + "-pairCorrCompare",
        xRange=params.Rrange,
        normalize="max",
        options=opts)

