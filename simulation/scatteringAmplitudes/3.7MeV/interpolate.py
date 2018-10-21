import numpy as np
from scipy.interpolate import interp1d


for atm in ["hydrogen","carbon", "oxygen","nitrogen","iodine"]:
  angStr = []
  sctStr = []

  print("opening: ",atm,atm + "_dcs.dat")
  with open(atm + "_dcs.dat", 'r') as inpFile:
    ind=0
    for line in inpFile:
      if ind < 31:
        ind += 1
        continue

      angStr.append(line[2:11])
      sctStr.append(line[39:50])

  angs = np.array(angStr).astype(np.float64)
  scts = np.sqrt(np.array(sctStr).astype(np.float64))

  f = interp1d(angs, scts, 'cubic')

  outAngs = np.arange(0, 0.7, 1e-4)
  outAngs = np.round(outAngs, 4)
  outScts = f(outAngs)

  np.savetxt("interpolationAngs.dat", outAngs)
  np.savetxt(atm+"_interpolation.dat", outScts)
  #outAngs.tofile("interpolationAngs.dat")
  #outScts.tofile(atm+"_interpolation.dat")
