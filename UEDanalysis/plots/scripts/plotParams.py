import numpy as np

class pltParams:
  def __init__(self): 

    self.Nlegendres     = 5
    self.QperPix        = 0.0223

    self.NradAzmBins    = 555
    self.NradLegBins    = 50
    self.QrangeAzm      = [0, self.QperPix*self.NradAzmBins]
    self.QrangeLeg      = [0, self.QperPix*self.NradLegBins]

    self.NpairCorrBins  = 369 
    self.Rrange         = [0, 9.86]
    #self.Rrange         = [0, 5.6354]
    #self.Rrange         = [0, self.NpairCorrBins
    #                       *(2*np.pi/(self.QperPix*(2*self.NradAzmBins - 1)))]

    self.smearStr  = "0.750000"
    self.smearList = ["1.000000", "0.750000", "0.500000", "0.250000", "0.100000", "0.050000", "0.025000"]

    self.mergeResultFolder = "/reg/ued/ana/scratch/nitroBenzene/mergeScans/"
