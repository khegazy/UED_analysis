import numpy as np
import sys
sys.path.append('../../../plots/scripts')
from plotClass import plotCLASS
from plotParams import pltParams

params = pltParams()

plc = plotCLASS()


#runs = ["20180629_1630", "20180627_1551", "20180630_1925", "20180701_0746"]
runs = ["20180627_1551", "20180629_1630", "20180630_1925", "20180701_0746"]
timeSteps = [29, 18, 19, 19]
maxX = [2, 1.1, 1.1, 1.1]


#################################################
#####  Plotting time dependent diffraction  #####
#################################################
#colRange = [[-2e-3, 2e-3], [-2e-3, 2e-3], [-2e-3, 2e-3], [-2e-3, 2e-3]]
#colRange = [[-7e-3, 7e-3], [-7e-3, 7e-3], [-7e-3, 7e-3], [-7e-3, 7e-3]]
colRange = [-3e-3, 3e-3]
for i,run in enumerate(runs):
  opts = {
    "colorRange" : colRange,
    "xTitle"     : "Time [ps]",
    "yTitle"     : r"Q [$\AA^{-1}$]",
    }

  timeDelay = np.fromfile("../../../mergeScans/results/timeDelays["
        + str(timeSteps[i] + 1) + "].dat", np.double)
  
  plc.print2d("/reg/ued/ana/scratch/nitroBenzene/mergeScans/data-" + run + "-azmAvgDiff["
        + str(timeSteps[i]) + "," + str(params.NradAzmBins) + "].dat",
        "../data-" + run + "-azmAvgDiffFull",
        X=timeDelay,
        yRange=params.QrangeAzm,
        options=opts)

  opts["xSlice"] = [-0.3, maxX[i]]
  plc.print2d("/reg/ued/ana/scratch/nitroBenzene/mergeScans/data-" + run + "-azmAvgDiff["
        + str(timeSteps[i]) + "," + str(params.NradAzmBins) + "].dat",
        "../data-" + run + "-azmAvgDiff",
        X=timeDelay,
        yRange=params.QrangeAzm,
        options=opts)


######################################################
#####  Plotting time dependent pair correlation  #####
######################################################

colRange = [-7e-2, 7e-2]#, [-1e-2, 1e-2], [-1e-2, 1e-2], [-1e-2, 1e-2]]
for i,run in enumerate(runs):
  opts = {
    "colorRange" : colRange,
    "xTitle"     : "Time [ps]",
    "yTitle"     : r"R [$\AA$]",
    #"interpolate": [200, 100]
    }

  timeDelay = np.fromfile("../../../mergeScans/results/timeDelays["
        + str(timeSteps[i] + 1) + "].dat", np.double)
  
  plc.print2d("../../results/data-" + run + "-pairCorrOdd["
        + str(timeSteps[i]) + "," + str(params.NpairCorrBins) + "].dat",
        "../data-" + run + "-pairCorrFull",
        X=timeDelay,
        yRange=params.Rrange,
        options=opts)

  opts["xSlice"] = [-0.3, maxX[i]]
  plc.print2d("../../results/data-" + run + "-pairCorrOdd["
        + str(timeSteps[i]) + "," + str(params.NpairCorrBins) + "].dat",
        "../data-" + run + "-pairCorr",
        X=timeDelay,
        yRange=params.Rrange,
        options=opts)

files = ["../../results/sim-phenoxyRadical-pairCorrOdd[" + str(params.NpairCorrBins) + "].dat"]
plc.print1d(files,
    "../sim-phenoxyRadical-pairCorrDiff",
    xRange=params.Rrange)
