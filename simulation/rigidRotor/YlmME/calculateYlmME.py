import sys
import os
import glob
import numpy as np
import py3nj
import argparse
from multiprocessing import Pool
from enum import Enum


parser = argparse.ArgumentParser()
parser.add_argument('L', type=int)
parser.add_argument('M', type=int)
parser.add_argument('maxJ', type=int)

args = parser.parse_args()

L = args.L
M = args.M
maxJ = args.maxJ
ovr_prefix = "YlmME_"
ind_prefix = "YlmME_indices_"

count = 0
overlap = []
indices = []

ovr_fileName =\
    ovr_prefix + "L-" + str(L)\
    + "_M-" + str(M)\
    + "_maxJ-" + str(maxJ) + ".dat"
ind_fileName =\
    ind_prefix + "L-" + str(L)\
    + "_M-" + str(M)\
    + "_maxJ-" + str(maxJ)

#if os.path.exists(ovr_fileName):
#  sys.exit(0)
for l1 in np.arange(maxJ+1):
  for m1_ in np.arange(-1*l1,l1+1):
    m1 = -1*m1_
    m2 = -1*(m1 + M)
    if np.amax(np.abs(np.array([m2, l1-L]))) > np.amin([maxJ, l1+L]):
      continue
    for l2 in np.arange(np.amax(np.abs(np.array([m2, l1-L]))), np.amin([maxJ, l1+L])+1):
      if np.abs(m2) > l2:
        continue

      integral = ((-1)**m1)*np.sqrt((2*l1+1)*(2*L+1)*(2*l2+1)/(4*np.pi))\
          *py3nj.wigner3j(2*l1,2*L,2*l2,0,0,0)\
          *py3nj.wigner3j(2*l1,2*L,2*l2,2*m1,2*M,2*m2)

      if integral:
        overlap.append(integral)
        indices.append(np.array([[l1, m1, l2, m2]]))

ovr = np.real(np.array(overlap)).astype(np.double)
ind = np.concatenate(indices, axis=0).astype(np.int32)

#for o,i in zip(ovr,ind):
#  print(o,i)

ovr.tofile(ovr_fileName)
ind_fileName += "_bins[" + str(ind.shape[0]) + "," + str(ind.shape[1]) + "].dat"
ind.tofile(ind_fileName)
print(ind[:5,:])

