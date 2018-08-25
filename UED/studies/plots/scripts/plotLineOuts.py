import numpy as np
import sys
sys.path.append('../../../plots/scripts')
from plotClass import plotCLASS
from plotParams import pltParams

params = pltParams()

plc = plotCLASS()

for i in range(1,121):
  files = ["../data/timeSlice-1548450_" + str(i) + ".dat", 
           "../data/timeSlice-1549200_" + str(i) + ".dat"]

  plc.print1d(files, 
        "../timeSlices-" + str(i),
        xRange=params.Qrange)


