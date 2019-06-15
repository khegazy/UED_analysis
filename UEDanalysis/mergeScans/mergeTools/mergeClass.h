#ifndef MERGECLASS_H
#define MERGECLASS_H


#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <ctime>

//Home Grown
#include "/reg/neh/home/khegazy/baseTools/tools/tools.h"
#include "/reg/neh/home/khegazy/baseTools/tools/imageProcessing.h"
#include "/reg/neh/home/khegazy/baseTools/tools/plotClass.h"
#include "/reg/neh/home/khegazy/baseTools/tools/saveClass.h"
#include "/reg/neh/home/khegazy/baseTools/tools/saving.h"
#include "/reg/neh/home/khegazy/baseTools/tools/parameters.h"


class mergeClass : public parameterClass {

  private :

    PLOTclass* plt;
  
    bool _compareReference;
    std::vector<double> compReference;

    int curScan;
    int NlegBins;
    int NrefBins;
    std::string runInd;
    std::string curDate;
    std::string curRun;
    std::string simFileNameSuffix;

    // Simulation results

    std::vector<double> Qazm;
    std::vector<double> Qleg;

    struct LTparamStruct {
      int scan;
      int64_t stagePos;
      double  imgNorm;
    };

    struct referenceStruct {
      std::vector<double> azmRef;
      std::vector<double> legRef;
      double imgNorm;
      double scale;
    };


    void getNrefBins();
    void initializeVariables();
    void calculateSMS();


  public :

    std::vector<double> atmLegDiff, atmAzmDiff;
    std::vector<double> molLegDiff, molAzmDiff;
    // Constructors/Destructor
    mergeClass(std::string runName);
    void compareReference(std::string);
    ~mergeClass();

    // Variables
    std::string runName;
    bool smearedTime;
    bool didSMSnormalize;
    bool didPairCorrSTD;
    bool didSubtractT0;


    // Output files and location
    string fileName = "mergedScans.txt";

    // Statistics
    bool SEMisBootstrap;
    double scanImgAzmRefMean, scanImgAzmRefSTD;
    double scanImgNormAzmRefMean, scanImgNormAzmRefSTD;
    std::vector<double> runLegRefMean, runAzmRefMean, runsMsRefMean;
    std::vector<double> runLegRefSTD, runAzmRefSTD, runsMsRefSTD;
    std::vector<double> runLegRefSEM, runAzmRefSEM, runsMsRefSEM;
    std::vector<double> scanImgAzmMeans, scanImgAzmSTDs;
    std::vector<double> scanImgNormAzmMeans, scanImgNormAzmSTDs;
    std::vector< std::vector<double> >  imgAzmRefMeans, imgAzmMeans;
    std::vector< std::vector<double> >  runLegMeans, runAzmMeans, 
                                        runsMsMeans, runPCorrMeans;
    std::vector< std::vector<double> >  runLegSEM, runAzmSEM, 
                                        runsMsSEM, runPCorrSEM;
    std::vector< std::vector<double> >  runLegSTD, runAzmSTD, 
                                        runsMsSTD, runPCorrSTD;

    // Data containers
    std::vector<double> timeDelays;
    std::vector<double> sMsLegNorm;
    std::vector<double> sMsAzmNorm;
    std::vector<double> azmReference;
    std::vector< std::vector<double> > azimuthalAvg, azimuthalsMs;
    std::vector< std::vector<double> > smearedAzmAvg, smearedAzmsMs;
    std::vector< std::vector<double> > legReference;
    std::vector< std::vector< std::vector<double> > > legendres, legendresMs;
    std::vector< std::vector< std::vector<double> > > smearedImg;
    std::map< int64_t, int > stagePosInds;
    std::map< int64_t, int > scanInds;
    std::map< int, std::vector<double> > scanScale, scanScale_init, scanImgNorms;
    std::map< int, std::vector< std::vector<double> > > scanLgndrs, scanLgndrs_init;
    std::map< int, std::vector< std::vector<double> > > scanAzmAvg, scanAzmAvg_init, scanAzmPCorr;
    std::map< int, std::map<std::string, double> > labTimeParams;
    std::map< int, std::pair<int, int> > labTimeMap;
    std::map< int, std::map< int, referenceStruct > > scanReferences, scanReferences_init;
    std::map< int, std::vector<double> > azmIndReference;


    std::map< int32_t, std::vector<double> > diffPs, legCoeff_map, pairCorr_map;
    std::map< int32_t, std::vector< std::vector<double> > > avgImgs_map;
    std::map< int32_t, double > counts;
    std::map< int32_t, int > Nimages;
    std::map< int32_t, std::vector<int32_t> > allPos;


    std::vector<double> smoothImgNorm, smoothImgNormSTD;

    void compareSimulations(std::vector<std::string> radicals);
    void addEntry(int scan, int64_t stagePos, int timeStamp,
                  std::vector<double>* azmAvg, 
                  std::vector<double>* legCoeffs,
                  double imgNorm);
    void addLabTimeParameter(int timeStamp,
                  std::string name, double value);
    void addReference(int scan, int64_t stagePos, int timeStamp,
                  std::vector<double>* azmAvg,
                  std::vector<double>* legCoeffs,
                  double imgNorm);
    void saveInitialData();
    void reloadInitialData();
   
    void removeLowPolynomials();
    void removeLabTimeOutliers();
    void removeOutliers();
    void removeImageOutliers();
    void removeImgNormOutliers();
    void stdParamCut(std::string paramName, double cut);
    void basicLessThanCut(std::string paramName, double cut);
    void basicGreaterThanCut(std::string paramName, double cut);

    void scaleByFit();
    void getRunMeanSTDSEM();
    void bootstrapSEM();
    void testSEMbootstrap();
    void getImageMeanSTD();
    void getImgNormMeanSTD();
    void mergeScans(bool refOnly=false, bool tdOnly=false);
    void subtractT0();
    void normalizeScansResults(); 
    void sMsNormalize();
    void smearTimeGaussian();
    void gaussianFilterQ();
    void smearTimeFFT();
    void makePairCorrs();


    /*
    int FFTsize = NradBins*2 + NautCpadding + 1;
    //int FFTsize = (NdiffInds + NautCpadding)*2 + 1;
    int FTsize = (int)(FFTsize/2) + 1;
    int indHR =  holeRat*NradBins;
    int outSize = FTsize*rMaxRat;
    */

};

#endif
