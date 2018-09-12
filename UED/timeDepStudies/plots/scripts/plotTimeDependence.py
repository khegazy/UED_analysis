import numpy as np
import sys
sys.path.append('../../../plots/scripts')
from plotClass import plotCLASS
from plotParams import pltParams

params = pltParams()

plc = plotCLASS()


runs = ["20180629_1630", "20180627_1551", "20180630_1925", "20180701_0746"]
timeSteps = [18, 29, 19, 19]
maxX = [1.1, 2, 1.1, 1.1]


#################################################
#####  Plotting time dependent diffraction  #####
#################################################

colRange = [[-2e-3, 2e-3], [-2e-3, 2e-3], [-2e-3, 2e-3], [-2e-3, 2e-3]]
for i,run in enumerate(runs):
  opts = {
    "colorRange" : colRange[i],
    "xTitle"     : "Time [ps]",
    "yTitle"     : r"Q [$\AA^{-1}$]",
    }

  timeDelay = np.fromfile("../../../mergeScans/results/timeDelays["
        + str(timeSteps[i] + 1) + "].dat", np.double)
  
  plc.print2d("/reg/ued/ana/scratch/nitroBenzene/mergeScans/data-" + run + "-sMsAzmAvgDiff["
        + str(timeSteps[i]) + "," + str(params.NradAzmBins) + "].dat",
        "../data-" + run + "-azmAvgDiffFull",
        X=timeDelay,
        yRange=params.QrangeAzm,
        options=opts)

  opts["xSlice"] = [-0.3, maxX[i]]
  plc.print2d("/reg/ued/ana/scratch/nitroBenzene/mergeScans/data-" + run + "-sMsAzmAvgDiff["
        + str(timeSteps[i]) + "," + str(params.NradAzmBins) + "].dat",
        "../data-" + run + "-azmAvgDiff",
        X=timeDelay,
        yRange=params.QrangeAzm,
        options=opts)


######################################################
#####  Plotting time dependent pair correlation  #####
######################################################

colRange = [[-1e-2, 1e-2], [-1e-2, 1e-2], [-1e-2, 1e-2], [-1e-2, 1e-2]]
for i,run in enumerate(runs):
  opts = {
    "colorRange" : colRange[i],
    "xTitle"     : "Time [ps]",
    "yTitle"     : r"R [$\AA$]",
    }

  timeDelay = np.fromfile("../../../mergeScans/results/timeDelays["
        + str(timeSteps[i] + 1) + "].dat", np.double)
  
  plc.print2d("../../results/data-" + run + "-pairCorr["
        + str(timeSteps[i]) + "," + str(params.NpairCorrBins) + "].dat",
        "../data-" + run + "-pairCorrFull",
        X=timeDelay,
        yRange=params.Rrange,
        options=opts)

  opts["xSlice"] = [-0.3, maxX[i]]
  plc.print2d("../../results/data-" + run + "-pairCorr["
        + str(timeSteps[i]) + "," + str(params.NpairCorrBins) + "].dat",
        "../data-" + run + "-pairCorr",
        X=timeDelay,
        yRange=params.Rrange,
        options=opts)


