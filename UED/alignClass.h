#ifndef ALIGNCLASS_H
#define ALIGNCLASS_H


#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <ctime>

//Home Grown
#include "/reg/neh/home/khegazy/baseScripts/tools.h"
#include "/reg/neh/home/khegazy/baseScripts/imageProcessing.h"
#include "/reg/neh/home/khegazy/baseScripts/plotClass.h"
#include "/reg/neh/home/khegazy/baseScripts/saveClass.h"
#include "/reg/neh/home/khegazy/baseScripts/saving.h"
#include "./parameters.h"


class alignClass : public parameterClass {
  public :

    // Constructor
    alignClass(std::string runName);
    ~alignClass();

    double* timeDelays;

    // Simulation results
    std::vector<double> atmLegDiff, atmAzmDiff;
    std::vector<double> molLegDiff, molAzmDiff;
    // Output files and location
    string fileName = "alignment.txt";
    string outputDir = "output/data/";

    std::vector<double> sMsLegNorm;
    std::vector<double> sMsAzmNorm;

    std::vector< std::vector<double> > azimuthalAvg;
    std::vector< std::vector< std::vector<double> > > legendres;
    std::vector< std::vector< std::vector<double> > > smearedImg;
    std::map< int32_t, int > stagePosInds;
    std::map< int, std::vector<double> > scanCounts;
    std::map< int, std::vector< std::vector<double> > > scanLgndrs;
    std::map< int, std::vector< std::vector<double> > > scanAzmAvg;


    std::map< int32_t, std::vector<double> > diffPs, legCoeff_map, pairCorr_map;
    std::map< int32_t, std::vector< std::vector<double> > > avgImgs_map;
    std::map< int32_t, double > counts;
    std::map< int32_t, int > Nimages;
    std::map< int32_t, std::vector<int32_t> > allPos;

    std::vector< std::vector<double> > runMeans;
    std::vector< std::vector<double> > runStdev;

    void compareSimulations(std::vector<std::string> radicals);
    void addEntry(int scan, int stagePos,
                  std::vector<double>* azmAvg, 
                  std::vector<double>* legCoeffs,
                  double imgNorm);


    void removeOutliers();
    void mergeScans();
    void subtractT0andNormalize();
    void smearTime();


    /*
    int FFTsize = NradBins*2 + NautCpadding + 1;
    //int FFTsize = (NdiffInds + NautCpadding)*2 + 1;
    int FTsize = (int)(FFTsize/2) + 1;
    int indHR =  holeRat*NradBins;
    int outSize = FTsize*rMaxRat;
    */

  private :

    PLOTclass* plt;
    
    int curScan;
    int legSize;
    int NrefBins;
    std::string runInd;
    std::string curDate;
    std::string curRun;
    std::string simFileNameSuffix;

    void getNrefBins();
    void initializeVariables();
    void calculateSMS();
};

#endif
