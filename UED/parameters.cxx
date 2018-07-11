#include "parameters.h"

parameterClass::parameterClass(std::string runName) {

  // Molecule
  molName = "nitrobenzene";

  // Image parameters
  Nlegendres = 1;
  NradBins = 50;
  imgSize = 895;
  hotPixel = 1750;
  hasRef = false;
  refStagePos = -1;

  hasLaserBkg = false;
  // Pair correlation parameters
  NautCpadding = 1000;
  holeRat = 0.15;
  rMaxRat = 0.75;
  padDecayRat = 0.5;


  subtractT0 = true;
  timeSmearStd = 0.025;

  // Simulation parameters
  compareSims     = false;
  Iebeam          = 5;
  elEnergy        = 3.7e6;
  screenDist      = 4;
  simReferenceDir = "/reg/ued/ana/scratch/nitroBenzene/simulations/references/";
  //std::string simReferenceDir = "../simulation/diffractionPattern/output/references/";

  preProcOutputDir = "./";//"/reg/ued/ana/scratch/nitroBenzene/rootFilesTest/";
  alignScansOutputDir = "/reg/ued/ana/scratch/nitroBenzene/alignScans/";
  backgroundImage = "NULL";
  backgroundFolder = "/reg/ued/ana/scratch/nitroBenzene/background/";
  //nanVal = 400;//1.23456789e12;
  pltCent = false;
  verbose = false;
  pltVerbose = true;

  scaleStagePos = 1e4;

  if (runName.compare("20161102_LongScan1") == 0) {

    NradBins = 30;
    imgMatType = "uint32";

    QperPix = 0.0233;
    maxQ = QperPix*imgSize/2;

    hasRef = true;
    refStagePos = 38.0;

    timeZero = 0.5;
    legStdCut = 1.5;

    /////  Center Finding Parameters  /////
    // Rough center finding
    sigma = 8;
    blockCentR = 560;
    blockCentC = 485;
    minRad = 90;
    maxRad = 340;
    meanInd = 335;
    stdScale = 0.025;

    // Fine center finding
    NavgCenters = 15;
    symVLeg = 1;
    minRadBin = 60;
    centShellW = 130; //70;
    holeR = 547;
    holeC = 492;
    holeRad = 60;
 
    cntrScale = 10;
    cntrMinScale = 1;
    cntrPowellTol = 0.2;
    cntrFracTol1d = 0.01;

    // Laser Background Removal Parameters
    hasLaserBkg = true;
    decayConst = -1.0/400.0;
    coreValThresh = 6.0e-3; //90; //65; //5e5;
    coreRad = 4; 
    clusterRad = 1; 
    minClusterSize = 550;
    minPixelSize = 300;
    minDensity = 0.2;//0.2;
    borderValThresh = 5.8e-3; //75; //1e5;
    borderRad = 3;
    padRad = 10;
 
    nanMap.resize(imgSize);
    for (int ir=0; ir<imgSize; ir++) {
      nanMap[ir].resize(imgSize, 0);
      if ((ir > 350) && (ir < 500)) {
        for (int ic=(int)imgSize*0.8; ic<imgSize; ic++) {
          nanMap[ir][ic] = NANVAL;
        }
      }
    }

    /*
    for (int ir=imgSize/2-40; ir<imgSize/2+40; ir++) {
      for (int ic=imgSize/2-120; ic<imgSize/2; ic++) {
        nanMap[ir][ic] = nanVal;
      }
    }
    for (int ir=imgSize/2-200; ir<imgSize/2+50; ir++) {
      for (int ic=imgSize-225; ic<imgSize-115; ic++) {
        nanMap[ir][ic] = nanVal;
      }
    }
    for (int ir=imgSize/2-60; ir<imgSize/2+35; ir++) {
      for (int ic=imgSize-75; ic<imgSize; ic++) {
        nanMap[ir][ic] = nanVal;
      }
    }
    */


  }
  else if ((runName.compare("20180701_0746") == 0)
            || (runName.compare("20180701_0738") == 0)) {

    QperPix = 0.0223;
    imgMatType = "uint16";

    //backgroundImage = "backgroundImg-20180701_0738.dat";

    hasRef = true;
    refStagePos = 154.0;

    /////  Center Finding Parameters  /////
    // Rough center finding
    sigma = 8;
    blockCentR = 560;
    blockCentC = 505;
    minRad = 70;
    maxRad = 325;
    meanInd = 350;
    stdScale = 0.075;

    // Fine center finding
    NavgCenters = 15;
    symVLeg = 1;
    minRadBin = 125;
    centShellW = 170;
    holeR = 590;
    holeC = 513;
    holeRad = 45;

    cntrScale = 10;
    cntrMinScale = 1;
    cntrPowellTol = 1;
    cntrFracTol1d = 0.01;
 
    timeZero = 8.8;
    legStdCut = 3.0;
  }
  else if ((runName.compare("20180629_1630") == 0)
            || (runName.compare("20180701_0738") == 0)) {

    QperPix = 0.0223;
    imgMatType = "uint16";

    //backgroundImage = "backgroundImg-20180701_0738.dat";

    hasRef = true;
    refStagePos = 154.0;

    /////  Center Finding Parameters  /////
    // Rough center finding
    sigma = 8;
    blockCentR = 560;
    blockCentC = 505;
    minRad = 70;
    maxRad = 325;
    meanInd = 350;
    stdScale = 0.075;

    // Fine center finding
    NavgCenters = 15;
    symVLeg = 1;
    minRadBin = 125;
    centShellW = 170;
    holeR = 590;
    holeC = 513;
    holeRad = 45;

    cntrScale = 10;
    cntrMinScale = 1;
    cntrPowellTol = 1;
    cntrFracTol1d = 0.01;
 
    timeZero = 8.8;
    legStdCut = 3.0;
  }
 
  else if (runName.compare("doRunLists") == 0) {
  }
  else if (runName.compare("other") == 0) {
  }
  else {
    std::cerr << "ERROR: Cannot set values for run " + runName + "!!!\n";
    exit(0);
  }

  refStagePos *= scaleStagePos;
  maxQ = QperPix*imgSize/2;
  maxR = (rMaxRat*NradBins)*(2*PI/(2*maxQ));

  std::cout << "maxQ " << maxQ<<std::endl;
}


