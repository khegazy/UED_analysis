import numpy as np
import sys
import os.path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from subprocess import call
#from fitFunctions import *
from fittingFuncts import *


if __name__ == '__main__':

  # Import command line arguments
  argP = argparse.ArgumentParser()
  argP.add_argument('--QperPix', type=float, default=0,
      help="Q density per pixel")

  args = argP.parse_args()

  NdataBins = 50
  scaleSimBins = 1
  NsimBins = NdataBins*scaleSimBins

  if args.QperPix:
    QperPix = args.QperPix
  else:
    QperPix = 0.0223

  Iebeam     = 5.
  elEnergy   = 3.7e6
  screenDist = 4.
  beamCut    = 4

  resultFolder = "/reg/neh/home5/khegazy/analysis/nitroBenzene/UED/qScale/results/"
  simFolder = "/reg/ued/ana/scratch/nitroBenzene/simulations/references/"
  simExeFile = "/reg/neh/home5/khegazy/analysis/nitroBenzene/simulation/diffractionPattern/simulateRefPatterns.exe"

  dataList = []
  #runs = ["20180701_0746", "20180629_1630"]#, "20180627_1551", "20180630_1925"]#, "20161102_LongScan1", "All"]
  runs = ["20180629_1630", "20180627_1551", "20180630_1925", "20180701_0746"]#, "20161102_LongScan1", "All"]
  #runs = ["20180627_1551", "20180630_1925"]#, "20161102_LongScan1", "All"]
  #runs = ["20180701_0746"]
  for rname in runs:
    fileName = "../mergeScans/results/referenceAzm-" + rname + "[555].dat"
    if os.path.isfile(fileName):
      dataList.append(np.fromfile(fileName, dtype = np.double))
    else:
      print("ERROR: Cannot find dataFile: " + fileName)
      sys.exit(0)


  ###  Get size of vector  ###
  size = dataList[0].shape[0]
  print size
  #maxQ = size*QperPix/0.86
  maxQ = size*QperPix


  #####  Import simulation  #####

  ###  Find zero crossings  ###
  atmDiffFile, molDiffFile = getSimNames(simFolder, size, 
      maxQ, Iebeam, screenDist, elEnergy)
 
  print("maxQ", maxQ)
  if not os.path.isfile(atmDiffFile):
    call([simExeFile,
          #"-mol", "nitrobenzene",
          "-maxQ", str(maxQ),
          #"-QperPix", str(QperPix),
          "-Odir", simFolder,
          "-Nbins", str(size)])

  atmDiff = np.fromfile(atmDiffFile, dtype=np.double)
  molDiff = np.fromfile(molDiffFile, dtype=np.double)
  zeroX = getZeroXQ(np.arange(size)*maxQ/size, molDiff) 


 
  ################################
  #####  Fitting background  #####
  ################################

  for idt,data in enumerate(dataList):
    Q = np.arange(size)*maxQ/size
    diffXvals = np.zeros(zeroX.shape[0])
    res = np.zeros(zeroX.shape[0])
    for i,q in enumerate(zeroX):
      ind = np.argmin(np.abs(Q - q))
      diffXvals[i] = atmDiff[ind]
      res[i] = data[ind]
      """
      if abs(Q[ind-1] - q) <= abs(Q[min(ind+1, Q.shape[0]-2)] - q):
        diffXvals[i] = atmDiff[ind-1] +\
                        (atmDiff[ind] - atmDiff[ind-1])/(Q[ind] - Q[ind-1])*(q - Q[ind-1])
        res[i]       = data[ind-1] +\
                        (data[ind] - data[ind-1])/(Q[ind] - Q[ind-1])*(q - Q[ind-1])
      else:
        diffXvals[i] = atmDiff[ind] +\
                        (atmDiff[ind+1] - atmDiff[ind])/(Q[ind+1] - Q[ind])*(q - Q[ind])
        res[i]       = data[ind] +\
                        (data[ind+1] - data[ind])/(Q[ind+1] - Q[ind])*(q - Q[ind])
      """
   
    print("zeroX",zeroX)
    print("res", res)
    fig, ax = plt.subplots()
    ax.plot(Q, molDiff)
    for q in zeroX:
      circ = plt.Circle((q, 0.0), 0.05, color='k')
      ax.add_artist(circ)
    ax.set_ylim([-0.25,0.3])

    fig.savefig("plots/zeroCrossings-" + runs[idt] + ".png") 
    print(diffXvals)

    


    nullst = zeroX
    null = res
    print("start fit")
    expfit = fit_scatt_bkg(null, nullst, Q)
    print("end fit")

    plt.figure()
    plt.semilogy(Q,data,linewidth=2,label='Experimental Diffraction Pattern')
    plt.semilogy(Q,expfit,linewidth=2,label='Exponential Fit')
    plt.semilogy(nullst,null,'o',color='red',markersize=10)
    #plt.ylim(0,600)
    plt.xlabel('Q [$\AA^{-1}$]')
    plt.ylabel('Diffraction Intensity / a.u.')
    plt.xlim(0,12)
    plt.legend(loc='best',fontsize=25)
    plt.tight_layout()
    plt.savefig("plots/expfit-" + runs[idt] + ".png")





    #sMolexp  = (patt-expfit)*Q/Phepatterns.I_at1Diso
    sMolexp  = (data-expfit)*Q/atmDiff
    sMolsim  = molDiff*Q/atmDiff
    sMolsim  *= np.min(sMolexp[size*0.25:size*0.5])/np.min(sMolsim[size*0.25:size*0.5]) 

    plt.figure()
    plt.plot(Q,sMolexp,linewidth=2,label='Experimental Molecular Diffraction')
    plt.plot(Q,sMolsim,linewidth=2,label='Simulated Molecular Diffraction')

    plt.plot(Q,np.zeros_like(Q),color='black')
    #plt.ylim(-1250,2000)
    plt.ylim(-30,60)
    if runs[idt] is "20161102_LongScan1":
      plt.ylim(-0.5e7,0.5e7)
    plt.xlabel('Q [$\AA^{-1}$]')
    plt.ylabel('Molecular Diffraction Intensity / a.u.')
    #plt.xlim(0,10.2)
    plt.xlim(0,maxQ)
    plt.legend(loc='best',fontsize=20)
    plt.tight_layout()
    plt.grid()
    plt.xticks(range(0, int(maxQ), 1))
    plt.savefig("plots/compareFit-" + runs[idt] + ".png")

