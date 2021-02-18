#ifndef PARAMETERCLASS_H
#define PARAMETERCLASS_H


#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <stdlib.h>
#include "constants.h"


class parameterClass {
  public:

    std::string experiment;
    std::string run;
    bool        hasI0;

    enum radicalEnum 
        {initialState, 
         finalState1, 
         finalState2, 
         finalState3,
         finalState4,
         finalState5};
    radicalEnum molecule;
    std::string molName;
    std::string home, scratch;

    std::vector<int> badScans;
    std::map< int, std::vector<int> > badImages;
    std::map< int, std::vector< std::pair<float, float> > > badRegions;
    int Nlegendres;
    int NradLegBins;
    int NmaxRadBins;
    int NradAzmBins;
    int NfinalPoints;

    // Diffraction parameters
    double QperPix;
    double maxQleg;
    double maxQazm;
    double timeZero;
    double legStdCut;
    int imgSize;
    int imgEdgeBuffer;
    int imgShutterTime;
    std::string imgMatType;
    bool hasLaserBkg;
    bool laserClusterRemoval;

    std::string preProcOutputDir;
    std::string preProcI0OutputDir;
    std::string mergeScansOutputDir;
    std::string scanSearchOutputDir;
    std::string radialPixelDist;
    std::string indexPath;


    bool hasRef;
    double refStagePosCut;
    bool subtractReference;
    float imgNormRadMin;
    float imgNormRadMax;
    std::vector<int> refSubtractStagePos;

    // Pair correlation parameters
    int NautCpadding;
    int NbinsSkip;
    int maxRbins;
    double filterVar;
    double holeRat;
    double rMaxLegRat;
    double rMaxAzmRat;
    double padDecayRat;
    double maxR;
    double lowQfillSimScale;
    double pCorrQcut;
    std::string fillLowQfile;
    bool fillLowQtheory;
    bool fillLowQzeros;
    bool fillLowQsine;
    bool fillLowQfitTheory;
    bool useFilledHole;

    bool pCorrGaussFilter;
    bool pCorrButterFilter;
    std::string pCorrFilterType;
    int   pCorrFilterOrder;
    float pCorrWnHigh;

    // Analysis parameters
    std::vector< std::vector<double> > signalRranges;
    std::vector< std::vector<double> > signalQranges;

    // Simulation parameters
    bool compareSims;
    bool simPairCorr;
    bool getBonds;
    bool fsFitOffset;
    bool fsFilterSMS;
    bool simPltVerbose;
    int NradSimBins;
    int fitQbegin;
    int fitQend;
    int fsQfitBegin;
    int fsQfitEnd;
    int fsRfitBegin;
    int fsRfitEnd;
    double fsFilterVar;
    double Iebeam;
    double elEnergy;
    double screenDist;

    bool    simHotFinalState;
    double  hotFSrefVar;
    double  hotFStdepVar;

    std::string xyzDir;
    std::string simOutputDir;
    std::vector<std::string> finalStates;
    std::vector<std::string> intermediateStates;
    std::vector<std::string> xyzFiles;
    std::vector<std::string> radicalNames;

    // Background removal
    int readoutAzmBinStart;
    int readoutAzmBinEnd;
    int readoutLegBinStart;
    int readoutLegBinEnd;
    double readoutStart;
    double readoutEnd;
    double bkgSTDcut;
    std::string backgroundFolder;
    std::string backgroundImage;
    std::string refCorrection;

    // Hot pixel removal
    double XrayHighCut;
    double XrayLowCut;
    double XraySTDcut; 
    int    XrayWindow; 
    bool   xRayHitDist;

    double hotPixel; 

    int    NshellOutlierLoops;
    int    shellWidth;
    int    Npoly;
    double stdChangeRatioLeft;
    double stdAccRatioLeft;
    double stdOutlierCutLeft;
    double stdCutLeft;
    double fracShellSTDcutLeft;
    double stdChangeRatioRight;
    double stdAccRatioRight;
    double stdOutlierCutRight;      
    double stdCutRight;      
    double fracShellSTDcutRight;
    double meanBinSize;
    double stdChangeRatio;
    double outlierSTDcut;
    bool   plotRadPixDist;

    double stdIncludeLeft;
    double stdIncludeRight;  
    double distSTDratioLeft; 
    double distSTDratioRight;

    double outlierSimpleSTDcut;

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
    bool scanAvgCenter;
    bool I0centers;
    bool computeCenters;
    int centerFxnType;
    int centerMinRadBin;
    int centerShellWidth;
    int holeR;
    int holeC;
    int holeRad;

    bool doCOM_center;
    double sigma;
    int centR_estimate;
    int centC_estimate;
    int minRad;
    int maxRad;
    int meanInd;
    std::vector<int> meanInds;
    double COMstdScale;

    double cntrScale;
    double cntrMinScale;
    double cntrPowellTol;
    double cntrFracTol1d;
    double centerSTDcut;

    int I0minPixVal;
    int I0approxR;
    int I0approxC;

    std::vector<double> I0ellRats;
    std::string centerDir;

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
    std::vector< std::vector<int> > nanMap;


    // Image Filtering Values
    int R1Bin;
    int R8Bin;
    int order;
    int suppressBins;
    float padMaxHeight;
    float WnLow;
    float WnHigh;
    bool  pltFilterVerbose;
    std::string filterType;

    // Remove low order polynomial noise
    int NlowOrderPoly;
    bool lowPolySubtractStudy;

    // Merging Scans
    bool  Qnormalize;
    bool  mergeNormalizeImgs;
    bool  useBootstrapSEM;
    bool  computeBootstrapSEM;
    bool  testMergeNbootStrap;
    int   timeFiltOrder;
    int   smearTimeBinWindow;
    int   mergeNbootstrap;
    float mergeSTDscale;
    float mergeImageSTDScale;
    float legImageNoiseCut;
    float azmImageNoiseCut;
    float labSTDcut;
    float labParamSmear;
    float timeWnHigh;
    float timeWnLow;
    float timeSmearSTD;
    float scanImgAzmSTDcut;
    float scanImgAzmRefSTDcut;
    std::string timeFilterType;

    bool  mergeGaussSmoothRef;
    float mergeGSmoothSTD;

    bool saveMergeIntermediates;
    std::string saveMergeInterFolder;


    // PV
    bool  getPVs;
    int   pvSampleTimes;
    float pressureSmear;
    float throttle;
    std::string pvFolder;
    std::map< std::string, std::string > pvMap;

    double gasShiftCut;

    // Power Scans
    double range1Qbegin;
    double range1Qend;
    double range2Qbegin;
    double range2Qend;

    // Time Zero
    std::vector< std::vector<double> > tZeroQranges;
    std::vector<int> tZeroRatio;

    // Misc
    std::vector< std::vector<double> > bkgStudyRanges;

    // Debugging
    bool verbose;
    bool pltCent;
    bool pltVerbose;

    parameterClass(std::string);
};

#endif
