import sys
import numpy as np
import argparse
import copy
from scipy.optimize import curve_fit
sys.path.append("/reg/neh/home/khegazy/baseTools/UEDanalysis/plots/scripts")
from plotClass import plotCLASS
import matplotlib.pyplot as plt
import matplotlib.colors as colors


if __name__ == '__main__':

  # Import command line arguments
  argP = argparse.ArgumentParser()
  argP.add_argument('--X', type=str, default='NULL',
      help='Function to fit to Y data.')
  argP.add_argument('--Y', type=str, default='NULL',
      help='Function to fit to.')
  argP.add_argument('--W', type=str, default='NULL',
      help='Standard deviation of Y data.')
  argP.add_argument('--plot', type=str, default='NULL',
      help='Plot name.')
  argP.add_argument('--scanCoeffs', type=str, default='NULL',
      help='Grid search of initial points.')
  argP.add_argument('--globMin', type=str, default='NULL',
      help='Grid search for global minimum.')
  args = argP.parse_args()

  plc = plotCLASS()

  X, xSize = plc.importImage(args.X)
  Y, ySize = plc.importImage(args.Y)
  W, wSize = plc.importImage(args.W)

  print("XSIZE", xSize)

  coeffs, errs = None, None
  if (xSize[0] == 1) or (len(xSize) == 1):
    def func(x, a):
      return a*x

    coeffs, errs = curve_fit(func, X, Y, sigma=W)

  elif xSize[0] == 2:
    def func(x, a, b):
      return a*x[0,:] + b*x[1,:]

    if args.scanCoeffs != 'NULL':
      c1 = np.linspace(0,0.5,num=10)
      c2 = np.linspace(0,0.5,num=10)
      fitV = np.zeros((2,len(c1),len(c2)))
      for i in range(len(c1)):
        for j in range(len(c2)):
          coeffs, errs = curve_fit(func, X, Y, sigma=W, 
              bounds=(0,100), p0=np.array([c1[i],c2[j]]))
          fitV[0,i,j] = coeffs[0]
          fitV[1,i,j] = coeffs[1]

      X,Y = np.meshgrid(c1,c2)
      plt.close()
      plt.pcolor(X,Y,fitV[0,:,:])
      plt.colorbar()
      plt.savefig(args.plot + "_coeffSmear1.png")
      plt.close()
      plt.pcolor(X,Y,fitV[1,:,:])
      plt.colorbar()
      plt.savefig(args.scanCoeffs + "_coeffSmear2.png")
    elif args.globMin != 'NULL':
      c1 = np.linspace(0,2.,num=30)
      c2 = np.linspace(0,2.,num=30)
      loss = np.zeros((len(c1),len(c2)))
      minLoss = np.inf
      bestCoeffs  = None
      bestErrs    = None
      for i in range(1, len(c1)):
        for j in range(1, len(c2)):
          coeffs, errs = curve_fit(func, X, Y, sigma=W,
              bounds=(0,100), p0=np.array([c1[i],c2[j]]))
          loss[i,j] = np.sqrt(np.mean((Y 
              - np.sum(np.reshape(coeffs, (-1,1))*X, axis=0))/W)**2)
          if loss[i,j] < minLoss:
            minLoss = loss[i,j]
            bestErrs   = copy.deepcopy(errs)
            bestCoeffs = copy.deepcopy(coeffs)
            bestStartVals = [c1[i],c2[j]]


      loss[0,0] = 0
      coeffs  = bestCoeffs
      errs    = bestErrs
      Xm,Ym   = np.meshgrid(c1,c2)
      plt.close()
      plt.pcolor(Xm, Ym, loss, 
          #norm=colors.LogNorm(vmin=minLoss, vmax=np.amax(loss)),
          vmin=minLoss,
          vmax=np.amax(loss))

      plt.xlabel("Coefficient 1")
      plt.ylabel("Coefficient 2")
      plt.colorbar()
      plt.savefig(args.globMin + "_loss.png")

      with open("./results/bestFitInitialVals.txt", 'a') as f:
        f.write(str(minLoss) 
            + "\t" + str(bestStartVals[0]) 
            + "\t" + str(bestStartVals[1]) + "\n") 
    else:
      coeffs, errs = curve_fit(func, X, Y, sigma=W, bounds=(0,10), p0=np.ones(2)*0.5)
      loss = np.sqrt(np.mean((Y 
          - np.sum(np.reshape(coeffs, (-1,1))*X, axis=0))/W)**2)
      with open("./results/bestFitDefaultInitialVals.txt", 'a') as f:
        f.write(str(loss) + "\n")


  elif xSize[0] == 3:
    def func(x, a, b, c):
      return a*x[0,:] + b*x[1,:] + c*x[2,:]

    coeffs, errs = curve_fit(func, X, Y, bounds=(0,100))

  elif xSize[0] == 4:
    def func(x, a, b, c, d):
      return a*x[0,:] + b*x[1,:] + c*x[2,:] + d*x[3,:]

  elif xSize[0] == 5:
    def func(x, a, b, c, d, e):
      return a*x[0,:] + b*x[1,:] + c*x[2,:] + d*x[3,:] + e*x[4,:]

  else:
    print("Currently do not support fitting N number of sims, please add")
    sys.exit(0)


  coeffs, errs = curve_fit(func, X, Y, bounds=(0,100))




  errs = np.sqrt(np.diag(errs))
  #print("coeff ",coeffs)
  #print("errs ", errs)

  coeffsFileName = "./results/fitCoefficients_Bins["\
      + str(coeffs.shape[0]) + "].dat"
  with open(coeffsFileName, "wb") as f: 
    coeffs.astype(np.double)
    coeffs.tofile(f) 

  errsFileName = "./results/fitCoefficientErrors_Bins["\
      + str(coeffs.shape[0]) + "].dat"
  with open(errsFileName, "wb") as f: 
    errs.astype(np.double)
    errs.tofile(f) 

  if args.plot != 'NULL':
    plt.plot(X[0,:])
    plt.plot(X[1,:])
    plt.savefig("testingInpFit.png")
    plt.close()

    plt.plot(Y)
    fit = np.zeros_like(Y)
    for i in range(coeffs.shape[0]):
      plt.plot(X[i,:]*coeffs[i])
      fit += X[i,:]*coeffs[i]
    plt.plot(fit)

    plt.savefig(args.plot + ".png")



