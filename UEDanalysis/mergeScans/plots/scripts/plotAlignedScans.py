import numpy as np
import sys
sys.path.append('../../../plots/scripts')
from plotClass import plotCLASS
from plotParams import pltParams

params = pltParams()

plc = plotCLASS()


runs = ["20180629_1630", "20180627_1551", "20180630_1925", "20180701_0746"]
timeSteps = [18, 29, 19, 19]
colRange = [[-2e-3, 2e-3], [-2e-3, 2e-3], [-2e-3, 2e-3], [-2e-3, 2e-3]]
maxX = [2, 10, 3, 3]

for i,run in enumerate(runs):
  opts = {
    "colorRange" : colRange[i],
    "xTitle"     : "Time [ps]",
    "yTitle"     : r"Q [$\AA^{-1}$]",
    "xSlice" : [-0.3, maxX[i]]
    }
  timeDelay = np.fromfile("../../results/timeDelays["
        + str(timeSteps[i] + 1) + "].dat", np.double)
  
  plc.print2d("/reg/ued/ana/scratch/nitroBenzene/mergeScans/data-" + run + "-sMsAzmAvgDiff["
        + str(timeSteps[i]) + "," + str(params.NradAzmBins) + "].dat",
        "../data-" + run + "-sMsAzmAvgDiff",
        X=timeDelay,
        yRange=params.QrangeAzm,
        options=opts)



