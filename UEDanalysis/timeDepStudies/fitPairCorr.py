import sys
sys.path.append("../plots/scripts/")
import os
from plotClass import plotCLASS
import tensorflow as tf
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt




class fitPairCorr():

  def __init__(self, parameters):

    ###############################
    #####  Import Parameters  #####
    ###############################

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
    if "singleAtoms" in parameters:
      self.singleAtoms  = parameters["singleAtoms"]
    else:
      raise RuntimeError("Must provide parameter `singleAtoms`.")
    if "Nbins" in parameters:
      self.Nbins = parameters["Nbins"]
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
    if "elEnergy" in parameters: 
      self.electron_energy = parameters["elEnergy"]  #eV
    else:
      raise RuntimeError("Must provide parameter `elEnergy`.")
    self.NbinsSkip    = 0
    self.NbondTypes   = np.math.factorial(len(self.atoms))\
                            - len(self.singleAtoms)
    self.scatAmps     = {}
    self.bondTypes    = []


    #####  Training  #####
    self.beta1 = 0.99
    self.beta2 = 0.999
    self.learning_rate = 1e-4


    #####  Variables  #####
    self.mask         = None
    self.global_step  = tf.Variable(
                            0, 
                            name="global_step", 
                            trainable=False)
    self.qPerPix      = tf.get_variable(
                            "qPerPix",
                            initializer=self.q_per_pixel,
                            trainable=False)
    self.elEnergy     = tf.get_variable(
                            "elEnergy",
                            initializer=self.electron_energy,
                            trainable=False)
    self.fitCoeffs    = tf.get_variable(
                          "fitCoeff", 
                          initializer=np.zeros((self.NbondTypes, self.NfitFxns),
                              dtype=np.float32))



    #####  Constants  #####
    self.scatAngs   = None
    self.scatAmps   = None
    self.atomicNorm = None
    self.NANVAL     = 1.234567e-10
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


    #####  Placeholders  #####
    self.fitData    = tf.placeholder(tf.float32, [self.Nbins], name="data")
    self.fitDataVar = tf.placeholder(tf.float32, [self.Nbins], name="var")
    self.fitMask    = tf.placeholder(tf.float32, [self.Nbins], name="mask")


    #####  Get Simulations  #####
    simDir = "/reg/ued/ana/scratch/nitroBenzene/simulations/"
    atmc  = "nitrobenzene_atmDiffractionPatternLineOut_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[555].dat"


    self.atomicScat = np.fromfile(simDir+atmc)
    self.atomicNorm = tf.constant(
                          1./np.expand_dims(self.atomicScat, 0),
                          dtype=tf.float32)


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
    self.add_optimizer()


  def add_data(self, data, variance):
    if data.shape[0] != self.Nbins:
      raise RuntimeError("Length of data (%i) and self.Nbins (%i) must match" 
          % (data.shape[0], self.Nbins))

    self.data = data
    self.var  = variance
    self.mask = np.ones(self.Nbins)

    ind = 0
    while self.data[ind] is self.NANVAL:
      self.mask[ind] = 0;
      ind += 1

          
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
    self.sinusiodsDR  = tf.divide(
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


  def add_loss(self):

    self.loss = tf.reduce_sum(
                    tf.divide(
                      tf.multiply(
                        self.fitMask,
                        tf.square(self.prediction - self.fitData)),
                      self.fitDataVar)
                    )
    self.loss /= tf.reduce_sum(self.fitMask)


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
        self.fitDataVar : self.var,
        self.fitMask    : self.mask}

    output = [
        self.global_step,
        self.loss,
        self.update]

    global_step, currentLoss, _ = sess.run(output, feed_dict)

    return global_step, currentLoss


  def fit(self, sess):

    self.history = []

    global_step = 0
    while (global_step < self.maxTrainIters) or (self.maxTrainIters is 0):
      global_step, currentLoss = self.run_fit_step(sess)
      self.history.append(currentLoss)
      print("Loss %i: \t%f" % (global_step, currentLoss))

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
        self.fitDataVar : self.var,
        self.fitMask    : self.mask}

    return sess.run(self.fitCoeffs, feed_dict)
