#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <ctime>

//From ROOT
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>

//From OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/objdetect/objdetect.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/contrib/contrib.hpp>

//Home Grown
#include "/reg/neh/home/khegazy/baseScripts/tools.h"
#include "/reg/neh/home/khegazy/baseScripts/imageProcessing.h"
#include "/reg/neh/home/khegazy/baseScripts/plotClass.h"
#include "/reg/neh/home/khegazy/baseScripts/saveClass.h"
#include "../simulation/diffractionPattern/simulations.h"
#include "parameters.h"
#include "alignClass.h"


using namespace std;

//////////////////////////////////////
//   Declaration of useful classes  //
//////////////////////////////////////

  PLOTclass plt;


//////////////////////////////////////////
//   Declaration of analysis variables  //
//////////////////////////////////////////

   TChain         *fChain;   //!pointer to the analyzed TTree or TChain 
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   string          *run;
   Int_t           scan;
   Int_t           imgNum;
   Int_t           imgIsRef;
   Int_t           timeStamp;
   Int_t           stagePos;
   Float_t         t0StagePos;
   Float_t         t0Time;
   Float_t         pressure;
   Float_t         pressureDer;
   Int_t           centerC;
   Int_t           centerR;
   Float_t         imgNorm;
   Float_t         readoutNoise;
   vector<vector<double> > *imgOrig;
   vector<vector<double> > *imgSubBkg;
   vector<double>  *legCoeffs;
   vector<double>  *azmAvg;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_scan;   //!
   TBranch        *b_imgNum;   //!
   TBranch        *b_imgIsRef;   //!
   TBranch        *b_timeStamp;   //!
   TBranch        *b_stagePos;   //!
   TBranch        *b_t0StagePos;   //!
   TBranch        *b_t0Time;   //!
   TBranch        *b_pressure;   //!
   TBranch        *b_pressureDer;   //!
   TBranch        *b_centerC;   //!
   TBranch        *b_centerR;   //!
   TBranch        *b_imgNorm;   //!
   TBranch        *b_readoutNoise;   //!
   TBranch        *b_imgOrig;   //!
   TBranch        *b_imgSubBkg;   //!
   TBranch        *b_legCoeffs;   //!
   TBranch        *b_azmAvg;   //!



//////////////////////////////////////
//  analysisClass for loading data  // 
//////////////////////////////////////

class analysisClass {

   public:
        analysisClass(string fileList);
        analysisClass(string fileList, string treeName);
        void initialize(string fileList, string treeName);
        ~analysisClass();
        uint64_t setupEnvironment();
        int loadEvent(uint64_t entry);

	clock_t start, stop;
};      


analysisClass::analysisClass(string fileList) {

  string treeName = "physics";
  initialize(fileList, treeName);
}


analysisClass::analysisClass(string fileList, string treeName) {

  initialize(fileList, treeName);
}
  

void analysisClass::initialize(string fileList, string treeName) {

  start=clock();

  ifstream files;
  files.open(fileList.c_str());
  if (!files.is_open()) {
    cerr<<"ERROR: Cannot open file list "<<fileList<<"!!!"<<endl;
    exit(0);
  } 
  
  ///// Initialize fChain (Chain of trees)
  fChain = new TChain(treeName.c_str());
  
  ///// Add trees from files in fileList to fChain
  string line;
  TFile* ftest;
  cout<<"!!!!!  Adding trees from files to fChain  !!!!!"<<endl;
  while (getline(files, line)) {
    cout<<"Adding tree '"<<treeName<<"' from file '"<<line<<"'"<<endl;
    ftest = TFile::Open(line.c_str());
    if (!ftest->Get(treeName.c_str())) cerr<<"WARNING: Cannot find tree "<<treeName<<" in file "<<line<<" !!!"<<endl;
    fChain->Add(line.c_str()); 
  }
  cout<<endl<<endl;

  ///// Set variable addresses to values

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("scan", &scan, &b_scan);
   fChain->SetBranchAddress("imgNum", &imgNum, &b_imgNum);
   fChain->SetBranchAddress("imgIsRef", &imgIsRef, &b_imgIsRef);
   fChain->SetBranchAddress("timeStamp", &timeStamp, &b_timeStamp);
   fChain->SetBranchAddress("stagePos", &stagePos, &b_stagePos);
   fChain->SetBranchAddress("t0StagePos", &t0StagePos, &b_t0StagePos);
   fChain->SetBranchAddress("t0Time", &t0Time, &b_t0Time);
   fChain->SetBranchAddress("pressure", &pressure, &b_pressure);
   fChain->SetBranchAddress("pressureDer", &pressureDer, &b_pressureDer);
   fChain->SetBranchAddress("centerC", &centerC, &b_centerC);
   fChain->SetBranchAddress("centerR", &centerR, &b_centerR);
   fChain->SetBranchAddress("imgNorm", &imgNorm, &b_imgNorm);
   fChain->SetBranchAddress("readoutNoise", &readoutNoise, &b_readoutNoise);
   fChain->SetBranchAddress("imgOrig", &imgOrig, &b_imgOrig);
   fChain->SetBranchAddress("imgSubBkg", &imgSubBkg, &b_imgSubBkg);
   fChain->SetBranchAddress("legCoeffs", &legCoeffs, &b_legCoeffs);
   fChain->SetBranchAddress("azmAvg", &azmAvg, &b_azmAvg);

}


analysisClass::~analysisClass() {

  delete fChain;
  stop = clock();
  cout<<endl<<"Total running time: "<<double(stop-start)/CLOCKS_PER_SEC<<endl<<endl;
}


uint64_t analysisClass::setupEnvironment() {

  /////////////////////////////////////////////////////////////////////
  //  Here goes anything needed to be setup before running the code  //
  /////////////////////////////////////////////////////////////////////

  loadEvent(0); 

  return fChain->GetEntries();
}


int analysisClass::loadEvent(uint64_t entry) {

  return fChain->GetEntry(entry);
}

  
