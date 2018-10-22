import sys
sys.path.append("../plots/scripts/")
from plotClass import plotCLASS
import tensorflow as tf
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt




class fitPairCorr():

  def __init__(self):

    self.q_per_pixel        = 0.0223
    self.electron_energy    = 3.7e6  #eV


    self.atoms        = ["carbon", "nitrogen", "oxygen"]
    self.singleAtoms  = ["nitrogen"]
    self.NfitFxns     = 1000
    self.Nbins        = 555
    self.NbinsSkip    = 0
    self.NbondTypes   = np.math.factorial(len(self.atoms))\
                            - len(self.singleAtoms)
    self.scatAmps     = {}
    self.bondTypes    = []


    #####  Variables  #####
    self.qPerPix  = tf.get_variable(
                        "qPerPix",
                        initializer=self.q_per_pixel,
                        trainable=False)
    self.elEnergy = tf.get_variable(
                        "elEnergy",
                        initializer=self.electron_energy,
                        trainable=False)


    #####  Constants  #####
    self.scatAngs   = None
    self.scatAmps   = None
    self.atomicNorm = None
    self.C  = tf.constant(299792458, dtype=tf.float32)
    self.h  = tf.constant(4.135667e-15, dtype=tf.float32)
    self.PI = tf.constant(np.pi, dtype=tf.float32)
    self.eMass = tf.constant(0.510999e6, dtype=tf.float32)
    self.rVals = tf.constant(
                     np.reshape(
                       np.linspace(1, 9, self.NfitFxns),
                       (1,-1)),
                     dtype=tf.float32)
    self.qBins = tf.constant(
                     np.tile(
                       np.expand_dims(
                         np.expand_dims(
                           np.arange(self.NbinsSkip, self.Nbins),
                           1),
                         0),
                       (len(self.atoms), 1, 1)),
                     dtype=tf.float32)


    #####  Get Data  #####
    simDir = "/reg/ued/ana/scratch/nitroBenzene/simulations/"
    phnxy = "phenoxyRadical_sMsPatternLineOut_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[555].dat"
    refer = "nitrobenzene_sMsPatternLineOut_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[555].dat"
    atmc  = "nitrobenzene_atmDiffractionPatternLineOut_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[555].dat"

    self.data = np.fromfile(simDir+phnxy, dtype=np.double)\
                    - np.fromfile(simDir+refer, dtype=np.double)

    self.atomicScat = np.fromfile(simDir+atmc)
    self.atomicNorm = tf.constant(
                          1./np.expand_dims(self.atomicScat, 0),
                          dtype=tf.float32)


    _scatAmps = None
    inpAngs = []
    with open("/reg/neh/home/khegazy/simulations/scatteringAmplitudes/3.7MeV/"
        + "interpolationAngs.dat", 'r') as inpFile:
      for line in inpFile:
        inpAngs.append(line)

    for atm in self.atoms:
      inpAmps = []
      with open("/reg/neh/home/khegazy/simulations/scatteringAmplitudes/3.7MeV/"
          + atm + "_interpolation.dat", 'r') as inpFile:
        for line in inpFile:
          inpAmps.append(line)
      """
      with open("/reg/neh/home/khegazy/simulations/scatteringAmplitudes/3.7MeV/"
          + atm + "_dcs.dat", 'r') as inpFile:
        ind = 0
        for line in inpFile:
          if ind < 31:
            ind += 1
            continue

          inpAmps.append(line[39:50])
          if atm == self.atoms[0]:
            inpAngs.append(np.sqrt(line[2:11].asfloat()))
      """
 
      if _scatAmps is None:
        _scatAmps = np.reshape(
                      np.array(inpAmps).astype(np.float32),
                      (1,-1))
      else:
        _scatAmps = np.concatenate((_scatAmps,
                      np.reshape(
                        np.array(inpAmps).astype(np.float32),
                        (1,-1))),
                      axis=0)

    print("Ishapes", len(inpAngs), _scatAmps.shape)
    self.scatAngs = tf.constant(
                        np.tile(
                          np.expand_dims(
                            np.expand_dims(
                              np.array(inpAngs).astype(np.float64), 1),
                            0),
                          (_scatAmps.shape[0],1,1)),
                        dtype=tf.float32)

    self.scatAmps = tf.constant(np.expand_dims(_scatAmps, axis=2),
                        dtype=tf.float32)
        
    


          
  def makeFit(self):

    self.fitCoeffs  = tf.get_variable(
                        "fitCoeff", 
                        initializer=np.zeros((self.NbondTypes, self.NfitFxns),
                            dtype=np.float32))

    self.deBrogW    = self.h*self.C/tf.sqrt(
                        tf.square(self.elEnergy + self.eMass)
                        - tf.square(self.eMass))*1e10

    self.qInp       = 4*self.PI/self.deBrogW\
                        *tf.sin((self.scatAngs/2)*(self.PI/180))

    self.qEval      = self.qPerPix*self.qBins

    print("SHAPE", self.qInp.shape.as_list())
    print("SHAPE", self.scatAmps.shape.as_list())
    print("SHAPE", self.qEval.shape.as_list())
    self.interp     = tf.contrib.image.interpolate_spline(
                          self.qInp,
                          self.scatAmps,
                          self.qEval,
                          order=1,
                          name="interpolation")

    self.bondScatAmpsList = []
    self.bondTypes        = []
    for i in range(len(self.atoms)):
      for j in range(i, len(self.atoms)):
        if (i == j) and (self.atoms[i] in self.singleAtoms):
          continue
        self.bondTypes.append(self.atoms[i] + "-" + self.atoms[j])
        print(self.bondTypes[-1])
        self.bondScatAmpsList.append(
            tf.multiply(
              self.interp[i,:,0], 
              self.interp[j,:,0],
              name="bondSA_" + self.bondTypes[-1]))

    self.bondScatAmps   = tf.stack(self.bondScatAmpsList)
    print("bondScatAmps ", self.bondScatAmps.shape.as_list())
    self.bondNormAmps = tf.multiply(
                            self.bondScatAmps,
                            self.atomicNorm)
    print("bondNormAmps ", self.bondNormAmps.shape.as_list())


    self.Sargs = tf.einsum(
                    'ij,jk->ik',
                    self.qEval[0],
                    self.rVals,
                    name = "Sargs")
    print("Sargs ", self.Sargs.shape.as_list())

    self.sinusiods    = tf.sin(self.Sargs)
    self.sinusiodsDR  = tf.multiply(
                              self.sinusiods,
                              self.rVals)
    self.fitFxns      = tf.multiply(
                            tf.expand_dims(self.sinusiodsDR, 0),
                            tf.expand_dims(self.bondNormAmps,2))
    print("fitFxns ", self.fitFxns.shape.as_list())

    self.prediction = tf.einsum(
                          'ijk,ik->j',
                          self.fitFxns,
                          self.fitCoeffs)




if __name__ == "__main__":
 
  plc = plotCLASS()
  debug = True

  fitCls = fitPairCorr()
  fitCls.makeFit()

  
  with tf.Session() as sess:

    sess.run(tf.global_variables_initializer())
    sess.run(tf.local_variables_initializer())

    if debug:
      output = [
          fitCls.deBrogW,
          fitCls.atomicNorm,
          fitCls.scatAmps,
          fitCls.qInp,
          fitCls.qEval,
          fitCls.bondScatAmps,
          fitCls.bondNormAmps,
          fitCls.interp,
          fitCls.Sargs,
          fitCls.sinusiods,
          fitCls.sinusiodsDR]

      deBrog, atmNorm, scatAmps, qInp, qEval, bondScatAmps,\
      bondNormAmps, interp, Sargs, sinusiods, sinusiodsDR\
          = sess.run(output)

      opts = {"yLog" : True}
      print("shape qs", qInp.shape, qEval.shape, bondScatAmps.shape)
      print("DeBroglie Wavelength (angs): ", deBrog) # forgot to m->angs
      plc.print1d(atmNorm, "./plots/testFit_atomicNorm", options=opts, isFile=False)
      opts["labels"] = fitCls.atoms
      plc.print1d(scatAmps, "./plots/testFit_scatteringAmps", options=opts, isFile=False)
      plc.print1d(qInp[0,:,0], "./plots/testFit_qInp", isFile=False)
      plc.print1d(qEval[0,:,0], "./plots/testFit_qEval", isFile=False)
      opts["labels"] = fitCls.bondTypes
      plc.print1d(bondScatAmps, "./plots/testFit_bondScatAmps", 
          options=opts, isFile=False)
      opts["labels"] = fitCls.bondTypes
      plc.print1d(bondNormAmps, "./plots/testFit_bondNormAmps", isFile=False)
      opts["labels"] = fitCls.atoms
      plc.print1d(interp, "./plots/testFit_interp", options=opts, isFile=False)
      plc.print2d(Sargs, "./plots/testFit_Sargs", isFile=False)
      plc.print2d(sinusiods, "./plots/testFit_sinusiods", isFile=False)



