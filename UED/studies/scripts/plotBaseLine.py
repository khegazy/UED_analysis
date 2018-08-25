import glob
import numpy as np
import sys
import os
import threading
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.append('../../../plots/scripts')
from plotClass import plotCLASS
from plotParams import pltParams


params = pltParams()
plc = plotCLASS()

runNames = ["20180627_1551", "20180629_1630", "20180630_1925", "20180701_0746"]
for run in runNames:

  readout, _ = plc.importImage("../data/readoutNoise_" + run + ".dat", False);
  diffAmplitude, _ = plc.importImage("../data/diffAmplitude_" + run + ".dat", False);
  imgNorm, _ = plc.importImage("../data/imgNorm_" + run + ".dat", False);
  pressure, _ = plc.importImage("../data/pressures_" + run + ".dat", False);
  print(pressure)
  pressureDers, _ = plc.importImage("../data/pressureDers_" + run + ".dat", False);

  fig = plt.figure(figsize=(10,12))
  
  plotGrid = (5,1)
  axPrD = plt.subplot2grid(plotGrid, (0,0))
  axPrD.plot(pressureDers)
  axPrs = plt.subplot2grid(plotGrid, (1,0))
  axPrs.plot(pressure)
  axRdO = plt.subplot2grid(plotGrid, (4,0))
  axRdO.plot(readout)
  axImN = plt.subplot2grid(plotGrid, (2,0))
  axImN.plot(imgNorm)
  axDfA = plt.subplot2grid(plotGrid, (3,0))
  axDfA.plot(diffAmplitude)

  plt.tight_layout()
  fig.savefig("../baseLine_" + run + ".png")
  plt.close()




