import sys
import numpy as np
import argparse
from scipy.optimize import curve_fit
sys.path.append("/reg/neh/home/khegazy/baseTools/UEDanalysis/plots/scripts")
from plotClass import plotCLASS
import matplotlib.pyplot as plt


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
  args = argP.parse_args()

  plc = plotCLASS()

  X, xSize = plc.importImage(args.X)
  Y, ySize = plc.importImage(args.Y)
  W, wSize = plc.importImage(args.W)

  coeffs, errs = None, None
  if (xSize[0] == 1) or (len(xSize) == 1):
    def func(x, a):
      return a*x

    coeffs, errs = curve_fit(func, X, Y, sigma=W)

  elif xSize[0] == 2:
    def func(x, a, b):
      return a*x[0,:] + b*x[1,:]

    coeffs, errs = curve_fit(func, X, Y, bounds=(0,100))

  elif xSize[0] == 3:
    def func(x, a, b, c):
      return a*x[0,:] + b*x[1,:] + c*x[2,:]

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



