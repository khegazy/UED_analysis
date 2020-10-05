import os, sys, glob
import numpy as np


def get_shape(fName):
    ipos = fName.find("Shape[") + 6
    fpos = fName.find("]")
    nums = fName[ipos:fpos].split(",")
    
    shape = []
    for n in nums:
        shape.append(int(n))
    return np.array(shape)


def get_images(
        molecule, run, scans=None, scan_ind=None, 
        base_dir="/reg/ued/ana/scratch/", skip=None, form='ts'):
    
    if scans is not None:
        if type(scans) is not list:
            scans = np.array([scans])
    base_dir = os.path.join(base_dir, molecule, "preProcessing", "Run-"+run)
    
    if scans is None:
        # Get all scans
        scans = []
        file_search = os.path.join(base_dir, "Scan*procImg_Shape*")
        for fl in glob.glob(file_search):
            ipos = fl.find("Scan-")+5
            scan = int(fl[ipos:fl.find("_", ipos)])
            if scan not in scans:
                scans.append(scan)
        scans = np.sort(np.array(scans))
    if scan_ind is not None:
        scans = np.array([scans[scan_ind]])

    if skip is not None:
        for sk in skip:
            scans = scans[scans!=sk]
    print("Scans: ", scans)

        
    # Get processed images
    print("Getting Images")
    images = []
    for scn in scans:
        file_search = os.path.join(base_dir, "Scan-{}_procImg_Shape*".format(scn))
        for fl in glob.glob(file_search):
            images.append(np.reshape(np.fromfile(fl, np.double), np.insert(get_shape(fl), 0, 1)))
    images = np.concatenate(images)

    # Get processed image nanMap
    print("Getting nanMaps")
    image_nanMaps = []
    for scn in scans:
        file_search = os.path.join(base_dir, "Scan-{}_procImg_nanM*".format(scn))
        for fl in glob.glob(file_search):
            image_nanMaps.append(np.reshape(np.fromfile(fl, np.int32), np.insert(get_shape(fl), 0, 1)))
    image_nanMaps = np.concatenate(image_nanMaps)
 
    # Get stage positions
    print("Getting stage positions")
    stage_pos = []
    for scn in scans:
        file_search = os.path.join(base_dir, "Scan-{}_stagePos_*".format(scn))
        for fl in glob.glob(file_search):
            stage_pos.append(np.expand_dims(np.fromfile(fl, np.int32), 0))
    stage_pos = np.concatenate(stage_pos)
    
    # Get image number
    print("Getting image number")
    image_num= []
    for scn in scans:
        file_search = os.path.join(base_dir, "Scan-{}_imgNum_*".format(scn))
        for fl in glob.glob(file_search):
            image_num.append(np.expand_dims(np.fromfile(fl, np.int32), 0))
    image_num = np.concatenate(image_num)
        
    # Get image norms
    print("Getting image normalization")
    image_norms = []
    for scn in scans:
        file_search = os.path.join(base_dir, "Scan-{}_imgNorm_*".format(scn))
        for fl in glob.glob(file_search):
            image_norms.append(np.expand_dims(np.fromfile(fl, np.float32), 0))
    image_norms = np.concatenate(image_norms)

    # Get image centers
    print("Getting image centers")
    image_centers = []
    for scn in scans:
        file_search = os.path.join(base_dir, "Scan-{}_centerR_*".format(scn))
        for fl in glob.glob(file_search):
            image_centers.append(np.expand_dims(np.concatenate([\
                np.expand_dims(np.fromfile(fl, np.int32), -1),
                np.expand_dims(np.fromfile(fl.replace("centerR", "centerC"), np.int32), -1)], axis=-1), 0))
    image_centers = np.concatenate(image_centers)

    # Get image center std ratios
    print("Getting image centers")
    image_center_stdRs = []
    for scn in scans:
        file_search = os.path.join(base_dir, "Scan-{}_centerRstd*".format(scn))
        for fl in glob.glob(file_search):
            image_center_stdRs.append(np.expand_dims(np.concatenate([\
                np.expand_dims(np.fromfile(fl, np.float32), -1),
                np.expand_dims(np.fromfile(fl.replace("centerR", "centerC"), np.float32), -1)], axis=-1), 0))
    image_center_stdRs = np.concatenate(image_center_stdRs)
    
    # Get image reference indicator
    print("Getting image reference indicator")
    image_isRef = []
    for scn in scans:
        file_search = os.path.join(base_dir, "Scan-{}_imgIsRef*".format(scn))
        for fl in glob.glob(file_search):
            image_isRef.append(np.expand_dims(np.fromfile(fl, np.int32), 0))
    image_isRef = np.concatenate(image_isRef)

    # Reformat dimensions
    if form == 'ts':
        images = np.transpose(images, (1,0,2,3))
        image_nanMaps = np.transpose(image_nanMaps, (1,0,2,3))
        stage_pos = np.transpose(stage_pos, (1,0))
        image_num = np.transpose(image_num, (1,0))
        image_norms = np.transpose(image_norms, (1,0))
        image_centers = np.transpose(image_centers, (1,0,2))
        image_center_stdRs = np.transpose(image_center_stdRs, (1,0,2))
        image_isRef = np.transpose(image_isRef, (1,0))
    elif form == 'st':
      pass
    else:
      print("Must specify a format 'st' or 'ts' where s=scans and t=time.")
      sys.exit(0)


    return scans, images, image_nanMaps, stage_pos, image_num, image_norms,\
            image_centers, image_center_stdRs, image_isRef
