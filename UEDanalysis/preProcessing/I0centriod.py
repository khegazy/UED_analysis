import numpy as np
import sys
import os
import glob
import argparse
from scipy.misc import imread
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
sys.path.append('../plots/scripts')
sys.path.append('/reg/neh/home/khegazy/baseTools/UEDanalysis/plots/scripts')
from plotClass import plotCLASS
from plotParams import pltParams

params = pltParams()

plc = plotCLASS()

outFolder = "/reg/ued/ana/scratch/nitroBenzene/I0/plots/"
opts = {
    "colorRange"  : [-200, 200],
    "colorMap"    : "seismic",
    #"colorNorm"   : "log"
    }

X,Y = np.meshgrid(np.arange(512), np.arange(512))
ellVals = [0.25, 0.5, 0.75, 0.85]
deltaVals = [50, 50, 40, 30]

def ellipsoid(x, **kwargs):
  return x[0]*(x[1] - kwargs["R"])**2 + x[2]*(x[3] - kwargs["C"])**2 - 1

def ellipsoidXYvals(x, Xi):
  dYi = np.sqrt((1 - x[2]*(x[3] - Xi)**2)/x[0])
  Xi = np.concatenate((Xi,np.flip(Xi, axis=0)), axis=0)
  Yi = np.concatenate((x[1] + dYi, np.flip(x[1] - 1*dYi, axis=0)), axis=0)
  return Xi, Yi


if __name__ == '__main__':
  # Import command line arguments
  argP = argparse.ArgumentParser()
  argP.add_argument('--fileName', type = str, default = None,
      help = "Name of the file to find the centroid")
  argP.add_argument('--minPixVal', type = int, default = 500,
      help = "Remove all pixels below in center of mass center")
  argP.add_argument('--approxR', type = int, default = None,
      help = "Approximate row value of the centroid")
  argP.add_argument('--approxC', type = int, default = None,
      help = "Approximate column value of the centroid")

  args = argP.parse_args()

  fileName  = args.fileName
  minPixVal = args.minPixVal
  approxR   = args.approxR
  approxC   = args.approxC

  if "Run" not in fileName:
    indI = fileName.find("Background") + 11
  else:
    indI = fileName.find("Run") + 4
  indF = fileName.find("/", indI+1)
  run = fileName[indI:indF]
  indI = fileName.find("scan") + 4
  indF = fileName.find("/", indI+1)
  scan = int(fileName[indI:indF])

  image = imread(fileName)
  plt.imshow(image, vmin=0, vmax=6000)
  plt.colorbar()
  smooth = gaussian_filter(image, sigma=4)

  # Center of mass center
  delta = 70
  smooth[smooth<minPixVal] = 0
  smooth[:,:approxC-delta] = 0
  smooth[:,approxC+delta:] = 0
  smooth[:approxR-delta,:] = 0
  smooth[approxR+delta:,:] = 0
  norm = np.sum(smooth)


  centerR = []
  centerC = []
  results = []
  beamMax = np.amax(smooth)
  for rat,dv in zip(ellVals, deltaVals):
    # Gather indices and values of pixels with similar values
    val = rat*beamMax 
    inds = np.where((smooth < val+dv) & (smooth > val-dv))
    kw = {
        "R" : inds[0],
        "C" : inds[1],
        "V" : smooth[inds[0],inds[1]]
        }
    #print(inds)
    #print(smooth[inds[0],inds[1]])
    #res = least_squares(ellipsoid, np.array([1e-2, 243, 1e-2, 157]), 
    res = least_squares(ellipsoid, np.array([1e-2, approxR, 1e-2, approxC]),
        kwargs=kw, verbose=0)
    if np.any(res.x < 0) or (res.cost == 0) or np.any(res.x < 1e-5):
      print("ERROR IN FITTING ELLIPSE ", scan, rat)
      plt.scatter(inds[1],inds[0])
      results.append(np.array([[np.nan]*4]))
      continue

    centerR.append(res.x[1])
    centerC.append(res.x[3])
    results.append(np.expand_dims(np.array(res.x), axis=0))

    dMax = np.sqrt(0.99/res.x[2])
    Xe = np.linspace(res.x[3]-dMax, res.x[3] + dMax, 500)
    Xe, Ye = ellipsoidXYvals(res.x, Xe)

    plt.scatter(inds[1],inds[0])
    plt.plot(Xe, Ye, "k-")
  
  ind = fileName.rfind("/") + 1
  name = fileName[ind:-4]
  plt.ylim(approxR-35, approxR+35)
  plt.xlim(approxC-35, approxC+35)
  plt.axes().set_aspect('equal')
  plt.savefig("/reg/ued/ana/scratch/nitroBenzene/I0/plots/ellFit/"
      + "run-" + run + "_scan-" 
      + str(scan) + "_" + name + ".png")
  plt.close()

  center = np.array([
      np.round(np.mean(centerR)),
      np.round(np.mean(centerC))]).astype(np.int32)
  ind = fileName.rfind("/")
  center.tofile(
      os.path.join("/reg/ued/ana/scratch/",
          "nitroBenzene/preProcessing/centers/temp",
          "centers_" + name + "[2].dat"))

  results = np.concatenate(results, axis=0).astype(np.double)
  results.tofile(
      os.path.join("/reg/ued/ana/scratch/",
          "nitroBenzene/preProcessing/centers/temp",
          "results_" + name\
              + "[{},{}].dat".format(results.shape[0], results.shape[1])))

  I0norm = np.sum(smooth[smooth > beamMax*0.2]).astype(np.double)
  I0norm.tofile(
      os.path.join("/reg/ued/ana/scratch/",
          "nitroBenzene/preProcessing/centers/temp",
          "norm_" + name + "[1].dat"))

