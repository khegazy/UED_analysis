
class pltParams:
  def __init__(self): 

    self.QperPix            = 0.0223
    self.NradAzmBins        = 555
    self.NpairCorrAzmBins   = 398
    self.NradLegBins        = 50
    self.NpairCorrLegBins   = 398
    self.Nlegendres         = 5

    self.QrangeAzm = [0, self.QperPix*self.NpairCorrAzmBins]
    self.RrangeAzm = [0, 6.255]
    self.QrangeLeg = [0, self.QperPix*self.NradLegBins]
    self.RrangeLeg = [0, 6]

    self.smearStr  = "0.750000"
    self.smearList = ["1.000000", "0.750000", "0.500000", "0.250000", "0.100000", "0.050000", "0.025000"]

