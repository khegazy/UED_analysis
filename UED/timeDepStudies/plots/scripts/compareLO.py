import numpy as np
from plotClass import plotCLASS
from plotParams import pltParams

params = pltParams()

plc = plotCLASS()

selectTimes = [1.8, 2.8, 4.5, 6]

times = np.fromfile("../data/timeDelays[72].dat", np.double)

timeInds = np.searchsorted(times, selectTimes)
timeInds[timeInds>=times.shape[0]-1] = -1
print(timeInds)

diffOpts = {
    "labels" : [str(tm) for tm in selectTimes],
    "xLabel" : r"Q [$\AA^{-1}$]",
    "ySlice" : [-5, 5]
    }

pcorOpts = {
    "labels" : [str(tm) for tm in selectTimes],
    "xLabel" : r"R [$\AA$]"
    }
phnoxyOpts = {
    "labels" : ["Phenoxy Simulation", "Data"],
    "ySlice" : [-6,6],
    }
phnylOpts = {
    "labels" : ["Phenyl Simulation", "Data"],
    "ySlice" : [-6,6],
    }
for std in ["0.750000"]: #smearSTD:

  files = ["../data/sim_phenoxyRadicalDiffractionDiff[" 
                + params.NradBins + "].dat",
           "../data/data_sMsFinalL0DiffSmear"
           + std + "[" + params.NradBins + "].dat"]
  plc.print1d(files, "../compare_phenoxysMsFinalDiff",
            xRange=params.Qrange,
            normalize=True,
            options=phnoxyOpts)

  files = ["../data/sim_phenylRadicalDiffractionDiff[" 
                + params.NradBins + "].dat",
           "../data/data_sMsFinalL0DiffSmear"
           + std + "[" + params.NradBins + "].dat"]
  plc.print1d(files, "../compare_phenylsMsFinalDiff",
            xRange=params.Qrange,
            normalize=True,
            options=phnylOpts)

  files = ["../data/sim_phenoxyRadicalPairCorrDiff[" 
                + params.NpairCorrBins1 + "].dat",
           "../data/data_pairCorrFinalDiffSmear"
                + std + "[" + params.NpairCorrBins + "].dat"]
  plc.print1d(files, "../compare_phenoxyPairCorrFinalDiff",
            xRange=params.Rrange,
            normalize=True,
            options=phnoxyOpts)

  files = ["../data/sim_phenylRadicalPairCorrDiff[" 
                + params.NpairCorrBins1 + "].dat",
           "../data/data_pairCorrFinalDiffSmear"
                + std + "[" + params.NpairCorrBins + "].dat"]
  plc.print1d(files, "../compare_phenylPairCorrFinalDiff",
            xRange=params.Rrange,
            normalize=True,
            options=phnylOpts)

  fileName = "../data/data_sMsL0DiffSmear" + std\
      + "[" + params.NtimeSteps + "," + params.NradBins + "].dat"
  plc.printLineOut(fileName, 0, timeInds, "../data_comparesMsLOSmear" + std +"", 
      xRange=params.Qrange, options=diffOpts)
  
  fileName = "../data/data_pairCorrDiffSmear" + std\
      + "[" + params.NtimeSteps + "," + params.NpairCorrBins + "].dat" 
  plc.printLineOut(fileName, 0, timeInds, "../data_comparePairCorrLOSmear" + std +"", 
      xRange=params.Rrange, options=pcorOpts)
