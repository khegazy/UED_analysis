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
    
    int curScan;
    int NlegBins;
    int NrefBins;
    std::string runInd;
    std::string curDate;
    std::string curRun;
    std::string simFileNameSuffix;

    // Simulation results
    std::vector<double> atmLegDiff, atmAzmDiff;
    std::vector<double> molLegDiff, molAzmDiff;

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
    };


    void getNrefBins();
    void initializeVariables();
    void calculateSMS();


  public :

    // Constructor
    mergeClass(std::string runName);
    ~mergeClass();

    std::string runName;
    double* timeDelays;

    // Output files and location
    string fileName = "mergedScans.txt";

    std::vector<double> sMsLegNorm;
    std::vector<double> sMsAzmNorm;
    std::vector<double> azmReference;
    std::vector<double> runLegRefMeans, runAzmRefMeans;
    std::vector<double> runLegRefSTD, runAzmRefSTD;
    std::vector< std::vector<double> > azimuthalAvg, azimuthalsMs;
    std::vector< std::vector<double> > legReference;
    std::vector< std::vector< std::vector<double> > > legendres, legendresMs;
    std::vector< std::vector< std::vector<double> > > smearedImg;
    std::map< int64_t, int > stagePosInds;
    std::map< int, std::vector<double> > scanCounts;
    std::map< int, std::vector< std::vector<double> > > scanLgndrs;
    std::map< int, std::vector< std::vector<double> > > scanAzmAvg;
    std::map< int, LTparamStruct > labTimeParams;
    std::map< int, std::map< int, referenceStruct > > scanReferences;


    std::map< int32_t, std::vector<double> > diffPs, legCoeff_map, pairCorr_map;
    std::map< int32_t, std::vector< std::vector<double> > > avgImgs_map;
    std::map< int32_t, double > counts;
    std::map< int32_t, int > Nimages;
    std::map< int32_t, std::vector<int32_t> > allPos;

    std::vector< std::vector<double> > runLegMeans, runAzmMeans;
    std::vector< std::vector<double> > runLegSTD, runAzmSTD;

    std::vector<double> smoothImgNorm, smoothImgNormSTD;

    void compareSimulations(std::vector<std::string> radicals);
    void addEntry(int scan, int64_t stagePos,
                  std::vector<double>* azmAvg, 
                  std::vector<double>* legCoeffs,
                  double imgNorm);
    void addLabTimeParameter(int timeStamp,
                  int scan, int64_t stagePos,
                  double imgNorm);
    void addReference(int scan, int64_t stagePos,
                  std::vector<double>* azmAvg,
                  std::vector<double>* legCoeffs,
                  double imgNorm);
   
    void removeLowPolynomials();
    void removeLabTimeOutliers();
    void removeOutliers();
    void getMeanSTD();
    void mergeScans(bool refOnly=false, bool tdOnly=false);
    void subtractT0();
    void normalize();
    void smearTimeGaussian();
    void smearTimeFFT();


    /*
    int FFTsize = NradBins*2 + NautCpadding + 1;
    //int FFTsize = (NdiffInds + NautCpadding)*2 + 1;
    int FTsize = (int)(FFTsize/2) + 1;
    int indHR =  holeRat*NradBins;
    int outSize = FTsize*rMaxRat;
    */

};

#endif
