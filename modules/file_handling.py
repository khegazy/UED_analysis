import os, sys, glob
import numpy as np

def get_shape_binary(file_name):
    if folder is not None:
        file_prefix = os.path.join(folder, file_prefix)
    files = glob.glob(file_prefix + "*")
    
    if len(files) == 0:
        print("ERROR: Cannot find file corresponding to " + file_prefix)
        sys.exit(1)
    elif len(files) > 1:
        print("ERROR: Found too many files corresponding to " + file_prefix)
        sys.exit(1)
    
    shape = files[0][files[0].find("[")+1:files[0].find("]")]
    shape = shape.split(",")
    for i in range(len(shape)):
        shape[i] = int(shape[i])
    
    return np.array(shape).astype(int)

