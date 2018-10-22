#ifndef PARAMETERCLASS_H
#define PARAMETERCLASS_H


#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <map>
#include "/reg/neh/home/khegazy/baseTools/tools/constants.h"


class parameterClass {
  public:

    enum radicalEnum 
        {nitrobenzene, 
         phenoxyRadical, 
         phenylRadical, 
         Nradicals};
    radicalEnum molecule;
    std::string molName;

    std::vector<int> badScans;
    std::map< int, std::vector<int> > badImages;
    int Nlegendres;
    int NradLegBins;
    int NmaxRadBins;
    int NradAzmBins;
    // Diffraction parameters
    double QperPix;
    double maxQleg;
    double maxQazm;
    double timeZero;
    double timeSmearStd;
    double legStdCut;
    double smearSTD;
    int imgSize;
    int buffer;
    int imgShutterTime;
    std::string imgMatType;
    bool hasLaserBkg;
    bool laserClusterRemoval;

    std::string preProcOutputDir;
    std::string mergeScansOutputDir;
    std::string scanSearchOutputDir;
    std::string radialPixelDist;
    std::string indexPath;


    bool hasRef;
    double refStagePos;
    bool subtractReference;
    float imgNormRadMin;
    float imgNormRadMax;

    // Pair correlation parameters
    int NautCpadding;
    int NbinsSkip;
    double filterVar;
    double holeRat;
    double rMaxLegRat;
    double rMaxAzmRat;
    double padDecayRat;
    double maxR;

    // Simulation parameters
    bool compareSims;
    bool simPairCorr;
    bool simPltVerbose;
    int NradSimBins;
    double Iebeam;
    double elEnergy;
    double screenDist;
    std::string xyzDir;
    std::string simOutputDir;
    std::vector<std::string> xyzFiles;
    std::vector<std::string> radicalNames;

    // Background removal
    std::string backgroundFolder;
    std::string backgroundImage;

    // Hot pixel removal
    double hotPixel; 
    int    shellWidth;
    int    Npoly;
    double stdIncludeLeft;
    double distSTDratioLeft; 
    double stdCutLeft;       
    int    meanBinSize;      
    double stdIncludeRight;  
    double distSTDratioRight;
    double stdChangeRatio;   
    double stdCutRight;      
    bool   outlierVerbose;
    double outlierMapSTDcut;
    double outlierCoreValThresh; //90; //65; //5e5;
    int    outlierCoreRad; 
    int    outlierClusterRad;
    int    outlierMinClusterSize;
    int    outlierMinPixelSize;
    double outlierMinDensity;//0.2;
    double outlierShapeVarCut;
    double outlierShapeEdgeCut;
    double outlierBorderValThresh; //75; //1e5;
    int    outlierBorderDistLimit;
    int    outlierBorderRad;
    int    outlierPadRad;
    int    outlierrMaxScale;
    int    outlierrMinScale;
    int    outliercMaxScale;
    int    outliercMinScale;

    std::string indicesPath;

    // Center finding
    int NavgCenters;
    int centerFxnType;
    int centerMinRadBin;
    int centerShellWidth;
    int holeR;
    int holeC;
    int holeRad;

    double sigma;
    int blockCentR;
    int blockCentC;
    int minRad;
    int maxRad;
    int meanInd;
    double COMstdScale;

    double cntrScale;
    double cntrMinScale;
    double cntrPowellTol;
    double cntrFracTol1d;

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


    // Image Filtering Values
    int R1Bin;
    int R8Bin;
    int order;
    int suppressBins;
    float padMaxHeight;
    float WnLow;
    float WnHigh;
    std::string filterType;

    // Merging Scans
    bool  Qnormalize;
    bool  normalizeImgs;
    float mergeSTDscale;
    float mergeImageSTDScale;
    float legImageNoiseCut;
    float azmImageNoiseCut;
    float labSTDcut;
    float labParamSmear;

    // PV
    int pvSampleTimes;
    float pressureSmear;
    std::string pressureFileName;

    // Misc
    double scaleStagePos;

    // Debugging
    bool verbose;
    bool pltCent;
    bool pltVerbose;

    parameterClass(std::string);
};

#endif
