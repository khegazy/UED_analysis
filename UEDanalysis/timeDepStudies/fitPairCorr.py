import sys
sys.path.append("../plots/scripts/")
import os
from plotClass import plotCLASS
import tensorflow as tf
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt




class fitPairCorr():

  def __init__(self, 
      parameters,
      _fitData, 
      _variance,
      _sineFit,
      debug=False):

    ###############################
    #####  Import Parameters  #####
    ###############################

    self.sineFit    = _sineFit
    self.parameters = parameters

    if "maxTrain" in parameters:
      self.maxTrainIters = parameters["maxTrain"] 
    else:
      raise RuntimeError("Must provide parameter `maxTrain`.")
    if "saveEvery" in parameters:
      self.saveEvery = parameters["saveEvery"]
    else:
      raise RuntimeError("Must provide parameter `saveEvery`.")
    if "atoms" in parameters:
      self.atoms = parameters["atoms"]
    else:
      raise RuntimeError("Must provide parameter `atoms`.")
    if "Natoms" in parameters:
      self.Natoms = parameters["Natoms"]
    else:
      raise RuntimeError("Must provide parameter `Natoms`.")
    if "Rrange" in parameters:
      self.minR = parameters["Rrange"][0]
      self.maxR = parameters["Rrange"][1]
    else:
      raise RuntimeError("Must provide parameter `singleAtoms`.")
    if "Nbins" in parameters:
      self.NinpBins = parameters["Nbins"]
    else:
      raise RuntimeError("Must provide parameter `Nbins`.")
    if "NfitFxns" in parameters:
      self.NfitFxns = parameters["NfitFxns"]
    else:
      raise RuntimeError("Must provide parameter `NfitFxns`.")
    if "qPerPix" in parameters:
      self.q_per_pixel = parameters["qPerPix"] 
    else:
      raise RuntimeError("Must provide parameter `qPerPix`.")
    if "doNormalEqn" in parameters:
      self.doNormalEqn = parameters["doNormalEqn"]
    else:
      raise RuntimeError("Must provide parameter `doNormalEqn`.")
    if "elEnergy" in parameters: 
      self.electron_energy = parameters["elEnergy"]  #eV
    else:
      raise RuntimeError("Must provide parameter `elEnergy`.")
    self.NbinsSkip    = 0
    self.scatAmps     = {}
    self.bondTypes    = []

    for i in range(len(self.atoms)):
      if self.atoms[i] == "hydrogen":
        continue
      for j in range(i, len(self.atoms)):
        if self.atoms[j] == "hydrogen":
          continue
        if (i == j) and (self.Natoms[i] == 1):
          continue
        self.bondTypes.append(self.atoms[i] + "-" + self.atoms[j])
    self.NbondTypes = len(self.bondTypes)


    #####  Data  #####
    self.NANVAL     = 1.234567e-10
    if _fitData.shape[0] != self.NinpBins:
      raise RuntimeError("Length of data (%i) and self.NinpBins (%i) must match" 
          % (_fitData.shape[0], self.NinpBins))

    self.startBin = 0
    while _fitData[self.startBin] is self.NANVAL:
      self.startBin += 1
    if "startBin" in parameters:
      self.startBin = parameters["startBin"]

    self.Nbins = self.NinpBins - self.startBin

    print("start", self.startBin)
    self.data = _fitData[self.startBin:].astype(np.float32)
    self.var  = _variance[self.startBin:].astype(np.float32)


    #####  Get Sine Transform  #####
    STfileName = parameters["sineTransFile"]
    STmaxR  = float(STfileName[STfileName.find("maxR")+5:\
                STfileName.find("[")])
    STbins  = int(STfileName[STfileName.find("[")+1:\
                STfileName.find("]")])
    print(STmaxR)
    self.inputSineTransform = np.fromfile(STfileName, dtype=np.double)
    self.sineTransformIntrp = interp1d(
                                np.linspace(0, 1, STbins)*STmaxR, 
                                self.inputSineTransform)
    self.sineTransform      = tf.constant(
                                self.sineTransformIntrp(
                                  np.linspace(self.minR, self.maxR, self.NfitFxns)),
                                dtype=tf.float32)
 

    #####  Training  #####
    self.beta1 = 0.99
    self.beta2 = 0.999
    self.learning_rate = 1e-4


    #####  Variables  #####
    self.global_step  = tf.Variable(
                            0, 
                            name="global_step", 
                            trainable=False)
    self.qPerPix      = tf.get_variable(
                            "qPerPix",
                            initializer=self.q_per_pixel,
                            trainable=False)
    self.fitCoeffs    = None



    #####  Constants  #####
    self.scatAngs   = None
    self.scatAmps   = None
    self.C  = tf.constant(299792458, dtype=tf.float32)
    self.h  = tf.constant(4.135667e-15, dtype=tf.float32)
    self.PI = tf.constant(np.pi, dtype=tf.float32)
    self.eMass    = tf.constant(0.510999e6, dtype=tf.float32)
    self.elEnergy = tf.constant(self.electron_energy)
    self.atomMult = tf.constant(self.Natoms, dtype=tf.float32)
    self.rVals    = tf.constant(
                      np.reshape(
                        np.linspace(self.minR, self.maxR, self.NfitFxns),
                        (-1,1)),
                      dtype=tf.float32)
    self.qBins    = tf.constant(
                      np.tile(
                        np.expand_dims(
                          np.expand_dims(
                            np.arange(self.startBin, self.NinpBins),
                            1),
                          0),
                        (len(self.atoms), 1, 1)),
                      dtype=tf.float32)


    #####  Placeholders  #####
    self.fitData    = tf.placeholder(tf.float32, [self.Nbins], name="data")
    self.fitDataVar = tf.placeholder(tf.float32, [self.Nbins], name="var")


      
    if not self.doNormalEqn:

      if self.sineFit:
        self.fitCoeffs  = tf.get_variable(
                            "fitCoeff",
                            initializer=np.zeros((1, self.NfitFxns), 
                                        dtype=np.float32))
      else:
        self.fitCoeffs  = tf.get_variable(
                            "fitCoeff", 
                            initializer=np.zeros(
                              (self.NbondTypes, self.NfitFxns), 
                              dtype=np.float32))



    #####  Get Simulations  #####
    _scatAmps = None
    inpAngs = []
    with open("/reg/neh/home/khegazy/baseTools/simulation/scatteringAmplitudes/3.7MeV/"
        + "interpolationAngs.dat", 'r') as inpFile:
      for line in inpFile:
        inpAngs.append(line)

    for atm in self.atoms:
      inpAmps = []
      with open("/reg/neh/home/khegazy/baseTools/simulation/scatteringAmplitudes/3.7MeV/"
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
        
   

    #########################
    #####  Build Model  #####
    #########################

    ##### Build Graph  #####
    self.makeFit()

    ##### Optimization  #####
    self.add_loss()
    if not self.doNormalEqn:
      self.add_optimizer()


          
  def makeFit(self):

    self.deBrogW    = self.h*self.C/tf.sqrt(
                        tf.square(self.elEnergy + self.eMass)
                        - tf.square(self.eMass))*1e10

    self.qInp       = 4*self.PI/self.deBrogW\
                        *tf.sin((self.scatAngs/2)*(self.PI/180))

    self.qEval      = self.qPerPix*self.qBins

    self.interp     = tf.contrib.image.interpolate_spline(
                          self.qInp,
                          self.scatAmps,
                          self.qEval,
                          order=1,
                          name="interpolation")

    self.bondScatAmpsList = []
    self.atomicScat       = []
    bondItr               = 0
    for i in range(len(self.atoms)):
      if self.atoms[i] == "hydrogen":
        continue
      self.atomicScat.append(self.atomMult[i]*tf.square(self.interp[i,:,0]))
      for j in range(i, len(self.atoms)):
        if self.atoms[j] == "hydrogen":
          continue
        if (i == j) and (self.Natoms[i] == 1):
          continue
        self.bondScatAmpsList.append(
            tf.multiply(
              self.interp[i,:,0], 
              self.interp[j,:,0],
              name="bondSA_" + self.bondTypes[bondItr]))
        bondItr += 1


    self.bondScatAmps   = tf.stack(self.bondScatAmpsList)
    self.molAtomicScat  = tf.reduce_sum(tf.stack(self.atomicScat), 
                              0, 
                              keepdims=True)
    print("bondScatAmps ", self.bondScatAmps.shape.as_list())
    print("molAtomicScat ", self.molAtomicScat.shape.as_list())
    self.bondNormAmps = tf.divide(
                            self.bondScatAmps,
                            self.molAtomicScat)
    print("bondNormAmps ", self.bondNormAmps.shape.as_list())


    self.Sargs = tf.einsum(
                    'ij,kj->ik',
                    self.rVals,
                    self.qEval[0],
                    name = "Sargs")
    print("Sargs ", self.Sargs.shape.as_list())

    if self.sineFit:
      self.fitFxns = tf.expand_dims(tf.sin(self.Sargs), 0)

    else:
      self.sinusiods    = tf.sin(self.Sargs)
      self.sinusiodsDR  = tf.divide(
                            self.sinusiods,
                            self.rVals)
      self.fitFxns      = tf.multiply(
                            tf.expand_dims(self.sinusiodsDR, 0),
                            tf.expand_dims(self.bondNormAmps,1))
    print("fitFxns ", self.fitFxns.shape.as_list())

    if self.doNormalEqn:
      self.X    = tf.reshape(self.fitFxns, 
                    (self.NbondTypes*self.NfitFxns, self.Nbins))
      self.varM = tf.eye(self.Nbins) #tf.diag(1./(self.var + 1e-3))
      self.NE = tf.matmul(self.X, tf.einsum('ij,kj->ik', self.varM, self.X))
      self.NE_norm    = tf.matrix_inverse(tf.matmul(self.X,
                            tf.einsum('ij,kj->ik', self.varM, self.X)),
                          name="normEqnNorm")
      self.Linv = tf.matmul(self.NE_norm, self.NE)
      self.Rinv = tf.matmul(self.NE, self.NE_norm)
      print("X / varM / NE_norm  ", 
          (self.X.shape.as_list(),
            self.varM.shape.as_list(),
            self.NE_norm.shape.as_list()))
      self.flatCoeffs = tf.matmul(self.NE_norm,
                          tf.matmul(self.X,
                            tf.matmul(self.varM, 
                              tf.expand_dims(self.fitData,1))),
                          name="flatCoeffs")
      self.fitCoeffs  = tf.reshape(self.flatCoeffs, 
                          (self.NbondTypes, self.NfitFxns),
                          name="fitCoeffs")
      print("fitCoeffs ", self.fitCoeffs.shape.as_list())

    self.prediction = tf.einsum(
                       'ijk,ij->k',
                        self.fitFxns,
                        self.fitCoeffs,
                        name="prediction")

    print("prediction ", self.prediction.shape.as_list())



  def add_loss(self):

    self.chiSq = tf.reduce_mean(
                    tf.divide(
                      tf.square(self.prediction - self.fitData),
                      self.fitDataVar))
    self.sineT = tf.reduce_mean(
                    tf.divide(
                      tf.square(
                        self.sineTransform - tf.reduce_sum(self.fitCoeffs, 0)),
                      tf.abs(self.sineTransform) + 1e-10))
    self.lowR  = tf.reduce_mean(
                    tf.matmul(
                      tf.abs(self.fitCoeffs),
                      self.rVals))

    self.smooth =   tf.reduce_mean(
                      tf.abs(
                        self.fitCoeffs[:,:-1] - self.fitCoeffs[:,1:]))
    self.smooth +=  0.5*tf.reduce_mean(
                      tf.abs(
                        self.fitCoeffs[:,:-2] - self.fitCoeffs[:,2:]))
    self.smooth *= 100*self.chiSq
    self.nonZero  = tf.reduce_mean(
                      tf.minimum(tf.zeros_like(self.fitCoeffs), 
                        self.fitCoeffs))
    self.nonZero  *= -1000
 
    self.noise    = tf.reduce_mean(
                      tf.square(
                        tf.where(
                          tf.less(self.fitCoeffs, 
                            0.01*tf.ones_like(self.fitCoeffs)), 
                          1 + tf.abs(self.fitCoeffs), 
                          tf.zeros_like(self.fitCoeffs))))
    """
    self.noise    = tf.reduce_mean(
                      tf.abs(
                        1./(0.0051 - tf.where(
                          tf.less(self.fitCoeffs, 
                            0.005*tf.ones_like(self.fitCoeffs)), 
                          tf.abs(self.fitCoeffs), 
                          tf.zeros_like(self.fitCoeffs)))))
    self.noise    /= 10000
    """

    #self.loss = self.chiSq + 100*self.sineT + self.lowR*tf.square(self.chiSq)
    self.loss = self.chiSq + self.smooth# + self.nonZero #+ self.noise


  def add_optimizer(self):

    with tf.variable_scope("optimizer"):
      self.solver = tf.train.AdamOptimizer(
                        learning_rate = self.learning_rate,
                        beta1         = self.beta1,
                        beta2         = self.beta2)

      self.update = self.solver.minimize(self.loss, global_step=self.global_step)


  def run_fit_step(self, sess):

    feed_dict = {
        self.fitData    : self.data,
        self.fitDataVar : self.var}

    output = [
        self.global_step,
        self.loss,
        self.chiSq,
        self.sineT,
        self.lowR,
        self.smooth,
        self.nonZero,
        self.noise,
        self.update]

    global_step, currentLoss,\
    currentChiSq, currentSineT,\
    currentLowR, currentSmooth,\
    currentNonZero, currentNoise, _ = sess.run(output, feed_dict)

    return global_step, currentLoss, currentChiSq, 100*currentSineT, currentLowR, currentSmooth, currentNonZero, currentNoise


  def fit(self, sess):

    self.history = []

    global_step = 0
    while (global_step < self.maxTrainIters) or (self.maxTrainIters is 0):
      global_step, currentLoss,\
      currentChiSq, currentSineT,\
      currentLowR, currentSmooth,\
      currentNonZero, currentNoise = self.run_fit_step(sess)
      self.history.append(currentLoss)
      print("Loss %i: \t%.2E \t%.2E \t%.2E \t%.2E \t%.2E \t%.2E" % 
          (global_step, currentLoss, currentChiSq, currentSineT, currentSmooth, currentNonZero, currentNoise))

      if global_step % self.saveEvery is 0:
        self.saver.save(sess,
            fileName=self.checkpointPath + "/fitting",
            global_step=global_step)

        pl.dump(self.history, 
            open(self.checkpointPath 
              + "/fitting-history-" 
              + str(global_step)+".pl", "wb"))

        if global_step is not 0:
          os.remove(self.checkpointPath 
              + "/fitting-history-"
              + str(global_step - self.saveEvery) + ".pl")



  def get_fitCoeff(self, sess):

    feed_dict = {
        self.fitData    : self.data,
        self.fitDataVar : self.var}

    return sess.run(self.fitCoeffs, feed_dict)


  def get_fit(self, sess):

    feed_dict = {
        self.fitData    : self.data,
        self.fitDataVar : self.var}

    if "startBin" in self.parameters:
      fit = np.zeros(self.parameters["Nbins"])
      fit[self.parameters["startBin"]:] = sess.run(self.prediction, feed_dict)
      return fit
    else:
      return sess.run(self.prediction, feed_dict)

  
  def evaluate(self, sess, output):
    
    feed_dict = {
        self.fitData    : self.data,
        self.fitDataVar : self.var}

    return sess.run(output, feed_dict)


  def debugFxn(self, sess):
    
    feed_dict = {
        self.fitData    : self.data,
        self.fitDataVar : self.var}

    output = [self.Linv, self.Rinv]

    return sess.run(output, feed_dict)


