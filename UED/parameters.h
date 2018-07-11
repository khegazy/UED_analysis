#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include "/reg/neh/home/khegazy/baseScripts/constants.h"


class parameterClass {
  public:

    std::string molName;

    int Nlegendres;
    int NradBins;
    // Diffraction parameters
    double QperPix;
    double maxQ;
    double timeZero;
    double timeSmearStd;
    double legStdCut;
    int imgSize;
    double hotPixel;
    std::string imgMatType;
    //double nanVal;
    bool hasLaserBkg;

    std::string preProcOutputDir;
    std::string alignScansOutputDir;

    bool hasRef;
    double refStagePos;
    bool subtractT0;

    // Pair correlation parameters
    int NautCpadding;
    double holeRat;
    double rMaxRat;
    double padDecayRat;
    double maxR;

    // Simulation parameters
    bool compareSims;
    double Iebeam;
    double elEnergy;
    double screenDist;
    std::string simReferenceDir;

    // Center finding
    int NavgCenters;
    int symVLeg;
    int minRadBin;
    int centShellW;
    int holeR;
    int holeC;
    int holeRad;

    double sigma;
    int blockCentR;
    int blockCentC;
    int minRad;
    int maxRad;
    int meanInd;
    double stdScale;

    double cntrScale = 10;
    double cntrMinScale = 1;
    double cntrPowellTol = 0.2;
    double cntrFracTol1d = 0.01;

    // Background removal
    std::string backgroundFolder;
    std::string backgroundImage;

    // Laser Background Removal Parameters
    double  decayConst;
    double  coreValThresh; //90; //65; //5e5;
    int     coreRad; 
    int     clusterRad;
    int     minClusterSize;
    int     minPixelSize;
    double  minDensity;//0.2;
    double  borderValThresh; //75; //1e5;
    int     borderRad;
    int     padRad;
    std::vector< std::vector<double> > nanMap;

 
    double scaleStagePos;


    // Debugging
    bool verbose;
    bool pltCent;
    bool pltVerbose;

    parameterClass(std::string);
};
