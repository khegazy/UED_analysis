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


stage = ["1530000", "1530500", "1542350", "1542450", "1542550", "1542650", "1542750"]
#for f in glob.glob("../data/polNormLO*1542350*"):
#  plc.print1d([f], "../" + f[8:-4], xRange=[0,2*3.14159])

opts = {
    "xTitle"   : "Pixel Value",
    "line"     : [[0,0], [0,0], "b", 4],
    "text"     : [0, 0, "Var: blay"]
}

#for f in glob.glob("../data/val*"):
#plc.printHist("../data/vals_1543050_425.dat", 40, "../vals_1543050_425")
#sys.exit(0)
#for f in glob.glob("../data/vals_*"):
folder = "/reg/ued/ana/scratch/nitroBenzene/radialPixelDist/"
for stg in stage:
  for f in glob.glob(folder + "data/vals_" + stg +"_*"):

    print(f)
    vInd = f.find("vals_")
    ipos = f.find("_", vInd+7) + 1
    fpos = f.find("_", ipos+1)
    print(f[ipos:fpos])
    rad = float(f[ipos:f.find("_", ipos+1)])
    ipos = fpos + 1
    fpos = f.find("_", ipos+1)
    print(f[ipos:fpos])
    meanL = float(f[ipos:fpos])
    ipos = fpos + 1
    fpos = f.find("_", ipos+1)
    print(f[ipos:fpos])
    var = float(f[ipos:fpos])
    ipos = fpos + 1
    print(f[ipos:-4])
    meanR = float(f[ipos:-4])
    if os.path.isfile(folder + f[8:-4]):
      continue

 
    opts["line"][0] = [meanL, meanL]
    opts["line"][1] = [0,100]
    fig,ax = plc.printHist(f, 50, folder + f[vInd:-4], 
        options=opts, returnPlt=True)

    l = Line2D([meanR, meanR], [0,100], color="r", linewidth=4)
    ax.add_line(l)

    xmin,xmax = ax.get_xlim()
    ymin,ymax = ax.get_ylim()
    ax.text(0.45*xmax, 0.9*ymax, "std: " + str(var))

    print('saving')
    fig.savefig(folder + "hists/" + f[vInd:-4] + ".png")
    plt.close()

#for f in glob.glob("../data/ang*"):
#  plc.printHist(f, 400, "../" + f[8:-4], options=opts)
#for f in glob.glob("../data/ind*"):
#  plc.printHist(f, 400, "../" + f[8:-4], options=opts)


"""
opts = {
    "colorMap"    : "Greys",
    "colorNorm"   : "log",
    "colorRange"    : [0.1, 1500]
}
for f in glob.glob("../data/outl*"):
  plc.print2d(f, "../" + f[8:-4], options=opts)
"""
"""

class threadPol(threading.Thread):
  def __init__(self,threadID, name, counter, files):
    threading.Thread.__init__(self)
    self.threadID = threadID
    self.name = name
    self.counter = counter
    self.files = files
    self.plc = plotCLASS()
    opts = {
      "xTitle"   : "Angle [rad]",
    }

  def run(self):
    for f in self.files:
      self.plc.print1d([f], "../" + f[8:-4], xRange=[0,2*3.14159])

class threadVals(threading.Thread):
  def __init__(self,threadID, name, counter, files):
    threading.Thread.__init__(self)
    self.threadID = threadID
    self.name = name
    self.counter = counter
    self.files = files
    self.plc = plotCLASS()
    opts = {
      "xTitle"   : "Pixel Value",
    }

  def run(self):
    for f in self.files:
      self.plc.printHist(f, 400, "../" + f[8:-4], options=opts)

Nthreads = 8
threads = []
polFiles = glob.glob("../data/pol*")
itr = float(len(polFiles))/(2*Nthreads - 1.)
for i in range(2*Nthreads):
  threads.append(threadPol(i, "Thread-" + str(i), i, 
      polFiles[int(i*itr):int((i+1)*itr)]))

i = 2*Nthreads
threads.append(threadPol(i, "Thread-" + str(i), i, 
    polFiles[int(i*itr):]))



valFiles = glob.glob("../data/val*")
sft = len(threads)
itr = float(len(valFiles))/(Nthreads - 1.)
for i in range(Nthreads):
  threads.append(threadVals(i + sft, "Thread-" + str(i + sft), i + sft, 
      valFiles[int(i*itr):int((i+1)*itr)]))

i = Nthreads
threads.append(threadVals(i + sft, "Thread-" + str(i + sft), i + sft, 
    valFiles[int(i*itr):]))

for th in threads:
  th.start()
"""






