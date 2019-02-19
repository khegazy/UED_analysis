#include "preProcessing.h"
#include "/reg/neh/home/khegazy/baseTools/tools/parameters.h"

using namespace std;



class centerfnctr {

    public:
        int fxnType;
        int minRadBin;
        int centShellW;
        imgProc::radProcTool* radProc;
        std::vector< std::vector<double> >* img;

        //centerfnctr();

        double operator() (std::vector<double> vect) {
          if (fxnType == 0)  return imgProc::centerSymXsqr(img, vect[0], vect[1], 8, 1, minRadBin, centShellW, true);
          if (fxnType == 1)  return imgProc::centerOddLeg(img, vect[0], vect[1], 1, minRadBin, centShellW, true);
          if (fxnType == 2)  {
            if (!radProc) {
              std::cerr << "ERROR: Must specifiy radProc in centerfnctr to use this option!!!\n";
              exit(0);
            }
            return radProc->radialSliceVar(img, vect[0], vect[1], minRadBin, centShellW);
          }
          else {
            cerr << "ERROR: Did not select a correct center finding alogrithm, now exiting!!!\n\n";
            exit(0);
          }
          //N2 return imgProc::centerOddLeg(img, vect[0], vect[1], 1, 70, 40, true, NANVAL);
        }
};
         

int main(int argc, char* argv[]) {


  // Make runLists
  bool doRunLists = false;

  bool localTesting = false;


  if (argc != 2) {
    cerr << "ERROR: Must run program a list of files" 
      << " ./preProcessing.exe runList.txt!!!" << endl;
    exit(0);
  }

  std::string runListName(argv[1]);
  bool doBackground = false;


  std::string runName;
  if (doRunLists) {
    runName = "doRunLists";
  }
  else {
    auto initPos = runListName.find("runList_") + 8;
    auto iPos = runListName.find("Run", initPos);
    if (iPos != string::npos) {
      iPos += 4;
    }
    else if (runListName.find("PowerScan") != string::npos) {
      iPos = runListName.find("PowerScan") + 10;
    }
    else if (runListName.find("Background") != string::npos) {
      std::cout << "INFO: Will do background!!!\n";
      iPos = runListName.find("Background") + 11;
      doBackground = true;
    }
    else {
      cerr << "ERROR: Cannot find runName!!!\n";
      exit(0);
    }
    runName = runListName.substr(iPos, 
        runListName.find("_Scan-") - iPos);
  }

  parameterClass params(runName);
  PLOTclass plt;



  const int Nlegendres = 1;
  const int NradLegBins = 50;

  // Indices
  //const int imgSize = 935;
  const int imgSize = 895;

  if (Nlegendres != params.Nlegendres) {
    cerr << "ERROR: parameter Nlegendres does not match with parameter class!!!\n";
    exit(0);
  }
  if (NradLegBins != params.NradLegBins) {
    cerr << "ERROR: parameter NradLegBins does not match with parameter class!!!\n";
    exit(0);
  }
  if (imgSize != params.imgSize) {
    cerr << "ERROR: parameter imgSize does not match with parameter class!!!\n";
    exit(0);
  }


   // Make sure the image has odd number of bins
  if (!(imgSize%2) || !(imgSize%2)) {
    cerr << "ERROR: imgSize and imgSize must be an odd number!!!" << endl;
    exit(0);
  }

  imgProc::radProcTool radProc(params.shellWidth, params.NradAzmBins);
  
  uint ifl;
  std::vector<PLOToptions> pltOpts(2);
  std::vector<string> pltVals(2);


  std::vector<double> legCoeffs(params.NradLegBins);
  std::vector<double> rawAzmAvg(params.NradAzmBins);
  std::vector<double> azmAvg(params.NradAzmBins);
  std::vector<double> azmCounts(params.NradAzmBins);

  std::vector< std::vector<double> > oddImgImgn;
  std::vector< std::vector<double> > symImg(imgSize);
  std::vector< std::vector<double> > stdRatioImg(imgSize);
  std::vector< std::vector<double> > outlierImage;
  std::vector< std::vector<double> > outlierBkg;
  for (uint ir=0; ir<imgSize; ir++) { 
    symImg[ir].resize(imgSize, 0);
    stdRatioImg[ir].resize(imgSize, 0);
  }

  fftw_complex* fftIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*imgSize*imgSize);
  fftw_complex* fftOut = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*imgSize*imgSize);
  fftw_plan fftFref = fftw_plan_dft_2d(imgSize, imgSize, fftIn, fftOut, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan fftBref = fftw_plan_dft_2d(imgSize, imgSize, fftIn, fftOut, FFTW_BACKWARD, FFTW_MEASURE);


  std::vector< std::vector<double> > radBins;

  centerfnctr centfnctr;
  std::vector<double> center(2);
  float centerC, centerR;
  int Ncents;

  double count;

  int imgIsRef;

  std::map< std::string, float > pvVals, pvValsDer;
  std::map< std::string, long int > pvStartTimes;
  std::map< std::string, std::vector<float> > pvAll, pvAllDer;
  float throttle;

  int timeStamp;

  int radInd;
  int imgNum;
  int curScan;
  std::string fileName;
  std::string line;
  size_t ipos;
  std::string date, scan, curRun, rFileName;
  int32_t stagePos;
  float t0SP, t0Time;
  std::vector<imgProc::imgInfoStruct> imgINFO;
  float readoutNoise;

  int imgNormBinMin = (int)(params.imgNormRadMin*params.NradAzmBins);
  int imgNormBinMax = (int)(params.imgNormRadMax*params.NradAzmBins);

  TFile* file=NULL;
  TTree* tree=NULL;

  // Reference Images
  int Nref = 0;
  std::vector<double> refAzmAvg(params.NradAzmBins);

  ///// Filter Parameters  /////
  int filtFFToutSize = (int)(params.NradAzmBins/2 + 1);
  double* qSpace = (double*) fftw_malloc(sizeof(double)*(params.NradAzmBins));
  fftw_complex* rSpace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*filtFFToutSize);
  fftw_plan filtFFTf = fftw_plan_dft_r2c_1d(params.NradAzmBins, qSpace, rSpace, FFTW_MEASURE);
  fftw_plan filtFFTb = fftw_plan_dft_c2r_1d(params.NradAzmBins, rSpace, qSpace, FFTW_MEASURE);

  std::string filterName = 
      "/reg/neh/home/khegazy/analysis/filters/" + params.filterType
      + "Filter_Bins-" + to_string(filtFFToutSize) 
      //+ "_WnLow-" + to_string(params.WnLow) 
      + "_WnHigh-"+ to_string(params.WnHigh) + ".dat";
  if (!tools::fileExists(filterName)) {
    cout << "\n\nINFO: Making new filter\n";
    system(("python /reg/neh/home/khegazy/analysis/filters/makeFilters.py --Nbins "
          + to_string(filtFFToutSize)
          + " --Ftype " + params.filterType 
          + " --Order " + to_string(params.order)
          //+ " --WnLow " + to_string(params.WnLow)
          + " --WnHigh " + to_string(params.WnHigh)).c_str());
  }

  int padRange = 0;
  int curPadRange = 0;
  std::vector<double> bandPassFilter(filtFFToutSize);
  save::importDat<double>(bandPassFilter, filterName);

  // Calculate sMs normalization
  std::vector<double> atmDiff(params.NradAzmBins);
  std::vector<double> sMsNorm(params.NradAzmBins);
  if (!doRunLists) {
    save::importDat<double>(atmDiff, params.simOutputDir 
      + "/" + params.molName 
      + "_atmDiffractionPatternLineOut_Qmax-" + to_string(params.maxQazm)
      + "_Ieb-" + to_string(params.Iebeam)
      + "_scrnD-" + to_string(params.screenDist)
      + "_elE-" + to_string(params.elEnergy)
      + "_Bins[" + to_string(params.NradAzmBins) + "].dat");

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      sMsNorm[iq] = (params.maxQazm*(iq + 0.5)/params.NradAzmBins)/atmDiff[iq];
    }
  }



  ////////////////////////////////
  ////  Retrieving File Info  ////
  ////////////////////////////////

  ppFunct::getScanRunInfo(imgINFO, runListName, params.verbose);
  curRun  = imgINFO[0].run;
  curScan = imgINFO[0].scan;


  //////////////////////////////////////////////////////
  /////  Creating smaller runList files if needed  /////
  //////////////////////////////////////////////////////

  if (doRunLists) {
    ppFunct::makeRunLists(imgINFO, runName, params.preProcOutputDir);
  }


  /////////////////////////////////////
  /////  Setting image variables  /////
  /////////////////////////////////////

  string imgAddr = imgINFO[0].path + imgINFO[0].fileName;
  cv::Mat imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
  int Nrows = imgMat.rows;
  int Ncols = imgMat.cols;
  std::vector< std::vector<double> > imgOrig(Nrows);
  std::vector< std::vector<double> > imgSubBkg(Nrows);
  std::vector< std::vector<double> > imgTemp(Nrows);
  std::vector< std::vector<double> > imgLaserBkg(Nrows);
  std::vector< std::vector<double> > imgBkg(Nrows);
  std::vector< std::vector<double> > imgCent;
  float imgNorm;
  for (int ir=0; ir<Nrows; ir++) {
    imgOrig[ir].resize(Ncols, 0);
    imgSubBkg[ir].resize(Ncols, 0);
    imgTemp[ir].resize(Ncols, 0);
    imgLaserBkg[ir].resize(Ncols, 0);
    imgBkg[ir].resize(1024, 0);
  }

  //////////////////////////////////////////////
  /////  Create and save background image  /////
  //////////////////////////////////////////////

  if (doBackground) {
    if (params.verbose) std::cout << "\n\nINFO: Making background!!!\n";

    std::vector< std::vector<double> > imgVec, bkgCount;
    for (ifl=0; ifl<imgINFO.size(); ifl++) {
      ///  Get image  ///
      imgAddr   = imgINFO[ifl].path + imgINFO[ifl].fileName;
      imgMat    = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
      Nrows = imgMat.rows;
      Ncols = imgMat.cols;

      if (ifl == 0) {
        bkgCount.resize(Nrows);
        for (uint ir=0; ir<Nrows; ir++) {
          bkgCount[ir].resize(Ncols, 0);
        }
      }

      cv::Mat imgSmooth = imgMat;
      /*
      cv::Mat imgSmooth(Nrows, Ncols, imgMat.type());
      cv::medianBlur(imgMat, imgSmooth, 5);
      cv::medianBlur(imgSmooth, imgMat, 5);
      cv::medianBlur(imgMat, imgSmooth, 5);
      cv::medianBlur(imgSmooth, imgMat, 5);
      cv::GaussianBlur(imgMat, imgSmooth, cvSize(13,13), 2, 2);
      */

      if (params.imgMatType.compare("uint16") == 0) {
        imgVec = imgProc::getImgVector<uint16_t>(imgSmooth, 
                    Nrows, Nrows/2, Ncols/2, true); 
      }  
      else if (params.imgMatType.compare("uint32") == 0) {
        imgVec = imgProc::getImgVector<uint32_t>(imgSmooth, 
                    Nrows, Nrows/2, Ncols/2, true);
      }  
      else {
        cerr << "ERROR: Do not recognize imgMatType = " 
            + params.imgMatType + "!!!\n";
        exit(0);
      }

      ///  Remove Xray hits  ///
      imgProc::removeXrayHits(
          &imgVec, 
          params.XrayHighCut, params.XrayLowCut, 
          params.XraySTDcut, params.XrayWindow);

      ///  Remove Spots  ///
      int rsCount = 0;
      double rsMean = 0;
      double rsSTD  = 0;
      for (int ir=0; ir<Nrows; ir++) {
        for (int ic=0; ic<Ncols; ic++) {
          if (ir >= 510 && ir < 640
              && ic >= 450 && ic < 600) {
            continue;
          }
          if (imgVec[ir][ic] == NANVAL) {
            continue;
          }

          rsMean += imgVec[ir][ic];
          rsCount++;
        }
      }
      rsMean /= rsCount;

      for (int ir=0; ir<Nrows; ir++) {
        for (int ic=0; ic<Ncols; ic++) {
          if (ir >= 510 && ir < 640
              && ic >= 450 && ic < 600) {
            continue;
          }
          if (imgVec[ir][ic] == NANVAL) {
            continue;
          }

          rsSTD += std::pow(rsMean - imgVec[ir][ic], 2);
        }
      }
      rsSTD = std::sqrt(rsSTD/rsCount);


      for (int ir=0; ir<Nrows; ir++) {
        for (int ic=0; ic<Ncols; ic++) {
          if (ir >= 510 && ir < 640
              && ic >= 450 && ic < 600) {
            continue;
          }
          if (imgVec[ir][ic] == NANVAL) {
            continue;
          }

          if (fabs(imgVec[ir][ic] - rsMean) > params.bkgSTDcut*rsSTD) {
            imgVec[ir][ic] = NANVAL;
          }
        }
      }

      ///  Add image  ///
      if (params.verbose) std::cout << "\tAdding images\n";
      for (int ir=0; ir<Nrows; ir++) {
        for (int ic=0; ic<Ncols; ic++) {
          if (imgVec[ir][ic] != NANVAL) {
            //if (imgINFO[ifl].stagePos == 1542222 && ir>=200 && ir<210 && ic==200) {
            //  cout<<imgVec[ir][ic]<<endl;
            //}
            imgBkg[ir][ic] += imgVec[ir][ic];
            bkgCount[ir][ic] += 1;
          }
        }
      }

      for (int ir=0; ir<Nrows; ir++) {
        for (int ic=0; ic<Ncols; ic++) {
          if (imgVec[ir][ic] == NANVAL) {
            imgVec[ir][ic] = -1;
          }
        }
      }
      save::saveDat<double>(imgVec, 
          "results/testBKG_" 
          + to_string(imgINFO[ifl].stagePos)
          + "[1024,1024].dat");
    }

    ///  Scaling vector  ///
    for (int ir=0; ir<Nrows; ir++) {
      for (int ic=0; ic<Ncols; ic++) {
        if (bkgCount[ir][ic] != 0) {
          imgBkg[ir][ic] /= bkgCount[ir][ic];
        }
        else {
          imgBkg[ir][ic] = NANVAL;
        }
      }
    }

    ///  Fill in pixels with NANVALS by averaging nearest neighbors  ///
    for (int ir=0; ir<Nrows; ir++) {
      for (int ic=0; ic<Ncols; ic++) {
        if (imgBkg[ir][ic] == NANVAL) {
          std::vector<double> collection;
          int ind = 1;
          while (collection.size() < 15) {
            if (ic - ind >= 0) {
              for (int irr=-1*ind; irr<=ind; irr++) {
                if ((ir + irr >= 0) && (ir + irr < (int)imgBkg.size())) {
                  if (imgBkg[ir+irr][ic-ind] != NANVAL) {
                    collection.push_back(imgBkg[ir+irr][ic-ind]);
                  }
                }
              }
            }
            if (ic + ind < (int)imgBkg.size()) {
              for (int irr=-1*ind; irr<=ind; irr++) {
                if ((ir + irr >= 0) && (ir + irr < (int)imgBkg.size())) {
                  if (imgBkg[ir+irr][ic+ind] != NANVAL) {
                    collection.push_back(imgBkg[ir+irr][ic+ind]);
                  }
                }
              }
            }
            if (ir - ind >= 0) {
              for (int icc=-1*ind+1; icc<=ind-1; icc++) {
                if ((ic + icc >= 0) && (ic + icc < (int)imgBkg.size())) {
                  if (imgBkg[ir-ind][ic+icc] != NANVAL) {
                    collection.push_back(imgBkg[ir-ind][ic+icc]);
                  }
                }
              }
            }
            if (ir + ind < (int)imgBkg.size()) {
              for (int icc=-1*ind+1; icc<=ind-1; icc++) {
                if ((ic + icc >= 0) && (ic + icc < (int)imgBkg.size())) {
                  if (imgBkg[ir+ind][ic+icc] != NANVAL) {
                    collection.push_back(imgBkg[ir+ind][ic+icc]);
                  }
                }
              }
            }
            ind++;
          }

          double mean = std::accumulate(collection.begin(), collection.end(), 0);
          mean /= collection.size();
          imgBkg[ir][ic] = mean;
        }
      }
    }


    save::saveDat<double>(imgBkg, "results/testBKGfinal[1024,1024].dat");

    ///  Save image  ///
    if (params.verbose) std::cout << "\tSaving image\n";
    save::saveDat<double>(imgBkg, params.backgroundFolder 
        + "/backgroundImg-" + imgINFO[0].run + ".dat");

    if (params.pltVerbose) {
      delete plt.printRC(imgBkg, "plots/background-" + imgINFO[0].run);
    }

    exit(0);
  }


  //////////////////////////////////
  /////  Get background image  /////
  //////////////////////////////////

  //cerr << "WARNING!!!!! NOT IMPORTING BKG"<<endl;
  if (params.backgroundImage.compare("NULL") != 0) {
    save::importDat<double>(imgBkg, params.backgroundFolder 
        + "/" + params.backgroundImage);
    if (params.pltVerbose) {
      delete plt.printRC(imgBkg, "plots/importBkg-" + imgINFO[0].run);
    }
  }


  ///////////////////////
  /////  Import PV  /////
  ///////////////////////

  if (params.getPVs) {
    double weight;
    std::vector<double> pvOrig; 
    long int pvStartTime = 0;
    for (auto const & pv : params.pvMap) {
      std::string pvFileName = pv.second;
      std::string pvName = pv.first;

      ipos = pvFileName.find("_");
      pvStartTime = 0;
      for (int i=1; i<stoi(pvFileName.substr(ipos+1, 2));i++) {
        pvStartTime += ppFunct::monthLengths[i]*24*3600;
      }
      ipos = pvFileName.find("_", ipos+1);
      pvStartTime += stoi(pvFileName.substr(ipos+1, 2))*24*3600;
      ipos = pvFileName.find("_", ipos+1);
      ipos = pvFileName.find("_", ipos+1);
      pvStartTime += stoi(pvFileName.substr(ipos+1, 2))*3600;
      ipos = pvFileName.find("_", ipos+1);
      pvStartTime += stoi(pvFileName.substr(ipos+1, 2))*60;
      ipos = pvFileName.find("-", ipos+1);
      pvStartTime += stoi(pvFileName.substr(ipos+1, 2));

      pvStartTimes[pvName] = pvStartTime;

      ipos = pvFileName.find("-", ipos+1);
      int pvSize = stoi(pvFileName.substr(ipos+1, pvFileName.length()-5-ipos));

      pvOrig.resize(pvSize, 0); 
      save::importDat(pvOrig, params.pvFolder+pvFileName);

      pvAll[pvName].resize(pvSize, 0);
      pvAllDer[pvName].resize(pvSize, 0);
      float norm;
      int NsmearSteps = 3*params.pressureSmear/params.pvSampleTimes + 1;
      for (int i=0; i<pvSize; i++) {
        pvAll[pvName][i] = pvOrig[i];
        norm = 1;
        for (int ii=1; ii<NsmearSteps; ii++) {
          weight = exp(-1*std::pow(params.pvSampleTimes*ii, 2)
                      /(2*pow(params.pressureSmear, 2)));
          if (i+ii < pvSize) {
            pvAll[pvName][i] += pvOrig[i+ii]*weight;
            norm += weight;
          }
          if (i-ii >= 0) {
            pvAll[pvName][i] += pvOrig[i-ii]*weight;
            norm += weight;
          }
        }
        pvAll[pvName][i] /= norm;
      }

      for (int i=0; i<pvSize; i++) {
        if ((i == 0) || (i == pvSize-1)) {
          pvAllDer[pvName][i] = 0;
        }
        else {
          pvAllDer[pvName][i] = ((pvAll[pvName][i+1] - pvAll[pvName][i-1])
                              /(2*params.pvSampleTimes));
        }
      }

      //plt.print1d(pvPressure, "./plots/pressureSmooth");
      //plt.print1d(pvPressureDer, "./plots/pressureDer");
      //plt.print1d(pvOrigPressure, "./plots/pressure");
      //exit(0);
    }
  }



  //////////////////////////////
  /////  Making root file  /////
  //////////////////////////////

   //rFileName = "/reg/ued/ana/scratch/nitroBenzene/rootFiles/" + dataType 
  std::string subFolder = "";
  if (runListName.find("Power", 0) != std::string::npos) {
    subFolder = "PowerScan";
  }
  /*
  else if (runListName.find("UVpump", 0) != std::string::npos) {
    subFolder = "UVpump";
  }
  else if (runListName.find("IRalign", 0) != std::string::npos) {
    subFolder = "IRalign";
  }
  else if (runListName.find("Au", 0) != std::string::npos) {
    subFolder = "Au";
  }
  else if (runListName.find("BKG", 0) != std::string::npos) {
    subFolder = "BKG";
  }
  */

  rFileName = params.preProcOutputDir + "/" + subFolder + "/"
  //rFileName =  "testing/" + subFolder + "/"
        + "Run-" + curRun + "_"
        + "Scan-" + to_string(curScan) + ".root";
  if (localTesting) {
    rFileName = "localTest.root";
  }

  if (params.verbose) std::cout << "\n\nINFO: Making file " << rFileName << endl;

  file = TFile::Open(rFileName.c_str(), "RECREATE");
  tree = new TTree("physics","physics");

  tree->Branch("run", 	        &curRun);
  tree->Branch("scan", 	        &curScan,       "scan/I");
  tree->Branch("imgNum", 	&imgNum, 	"imgNum/I");
  tree->Branch("imgIsRef", 	&imgIsRef, 	"imgIsRef/I");
  tree->Branch("timeStamp", 	&timeStamp, 	"timeStamp/I");
  tree->Branch("stagePos", 	&stagePos, 	"stagePos/I");
  tree->Branch("t0StagePos",    &t0SP,		"t0StagePos/F");
  tree->Branch("t0Time",	&t0Time,	"t0Time/F");
  tree->Branch("throttle",      &throttle,      "throttle/F");
  tree->Branch("centerC", 	&centerC, 	"centerC/I");
  tree->Branch("centerR", 	&centerR, 	"centerR/I");
  tree->Branch("imgNorm", 	&imgNorm, 	"imgNorm/F");
  tree->Branch("readoutNoise", 	&readoutNoise, 	"readoutNoise/F");
  tree->Branch("imgOrig", 	&imgOrig);
  tree->Branch("imgSubBkg",     &imgSubBkg);
  tree->Branch("legCoeffs",     &legCoeffs);
  tree->Branch("rawAzmAvg",     &rawAzmAvg);
  tree->Branch("azmAvg",    &azmAvg);

  for (auto const & pv : params.pvMap) {
    pvVals[pv.first] = 0;
    pvValsDer[pv.first] = 0;
    tree->Branch(
        pv.first.c_str(),	  
        &pvVals[pv.first],	
        (pv.first+"/F").c_str());
    tree->Branch(
        (pv.first+"Der").c_str(),  
        &pvValsDer[pv.first],	
        (pv.first+"Der/F").c_str());
  }
  if (params.verbose) cout << "INFO: Tree and file are setup!\n";


  //////////////////////////////////
  /////  Checking if bad scan  /////
  //////////////////////////////////

  /*
  if (std::find(params.badScans.begin(), params.badScans.end(), imgINFO[ifl].scan)
        != params.badScans.end()) {
    std::cout << "INFO: This is a bad scan, saving empty tree and exiting!!!\n";
    tree->Write();
    file->Close();
    exit(1);
  }
  */
  
  //////////////////////////////////
  /////  Rough center finding  /////
  //////////////////////////////////
  std::vector<int> centersCOM;
  /*
  if (params.verbose) cout << "\n\nINFO: Start COM center finding\n"; 
  std::vector<int> centersCOM = imgProc::centerSearchCOM(imgINFO,
      params.hotPixel, params.sigma,
      params.blockCentR, params.blockCentC, params.minRad, params.maxRad,
      params.meanInd, params.COMstdScale, params.verbose, NULL); //&plt);

  if (params.verbose) std::cout << "COMcenters: " 
    << centersCOM[0] << "  " << centersCOM[1] << endl;
    */


  std::map<int, vector<float> > centers;
  ifstream infile("centers.txt");
  while (infile) {
    string s;
    if (!getline( infile, s )) break;

    istringstream ss( s );

    string sss;
    if (!getline( ss, sss, ' ' )) break;
    if (stoi(sss) != imgINFO[0].scan) continue;

    string pos;
    getline( ss, pos, ' ' );
    
    while (ss) {
      if (!getline( ss, sss, ' ' )) break;

      centers[stoi(pos)].push_back(stof(sss));
    }
  }


  /////////////////////////////////
  /////  Fine center finding  /////
  /////////////////////////////////

  if (localTesting || true) {
    /*
    centersCOM[0] = 568;
    centersCOM[1] = 495; //494;
    centerR = centersCOM[0];
    centerC = centersCOM[1];
    cout<<"centers "<<centerR<<"  "<<centerC<<endl;
    */
  }
  else {
    if (params.verbose) cout << "\n\nINFO: Starting fine center Finding!\n";
    centerR = centerC = Ncents = 0;
    for (int icnt=0; icnt<std::min(params.NavgCenters, (int)imgINFO.size()); icnt++) {

      /// Filling image std::vector ///
      imgAddr = imgINFO[icnt].path + imgINFO[icnt].fileName;
      imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
      cout<<"centering image addr: "<<imgAddr<<endl;

      Nrows = imgMat.rows;
      Ncols = imgMat.cols;
      if (params.imgMatType.compare("uint16") == 0) {
        if (params.hasLaserBkg) {
          imgCent = imgProc::getImgVector<uint16_t>(imgMat, 
                      Nrows, Nrows/2, Ncols/2, true, &imgBkg,
                      params.holeRad, params.holeR, params.holeC,
                      &params.nanMap);
        }
        else {
          imgCent = imgProc::getImgVector<uint16_t>(imgMat, 
                      Nrows, Nrows/2, Ncols/2, true, &imgBkg,
                      params.holeRad, params.holeR, params.holeC);
        }
      }  
      else if (params.imgMatType.compare("uint32") == 0) {
        if (params.hasLaserBkg) {
          imgCent = imgProc::getImgVector<uint32_t>(imgMat, 
                      Nrows, Nrows/2, Ncols/2, true, &imgBkg,
                      params.holeRad, params.holeR, params.holeC,
                      &params.nanMap);
        }
        else {
          imgCent = imgProc::getImgVector<uint32_t>(imgMat, 
                      Nrows, Nrows/2, Ncols/2, true, &imgBkg,
                      params.holeRad, params.holeR, params.holeC);
        }
      }  
      else {
        cerr << "ERROR: Do not recognize imgMatType = " 
            + params.imgMatType + "!!!\n";
        exit(0);
      }

      // Remove Xray hits
      imgProc::removeXrayHits(
          &imgCent, 
          params.XrayHighCut, params.XrayLowCut, 
          params.XraySTDcut, params.XrayWindow);

      //  Remove pixel outliers
      radProc.removeOutliers(imgCent,  
          centersCOM[0], centersCOM[1], params.imgEdgeBuffer,
          500, params.shellWidth, params.Npoly,
          params.stdIncludeLeft, params.distSTDratioLeft,
          params.stdCutLeft, params.meanBinSize,
          params.stdIncludeRight, params.distSTDratioRight,
          params.stdChangeRatio, params.stdCutRight,
          imgINFO[icnt].stagePos, params.outlierMapSTDcut,
          false, params.outlierVerbose, NULL);

    
      // Remove readout noise
      if (params.verbose) std::cout << "\tSubtract readout noise.\n";
      imgProc::removeReadOutNoise(imgCent);

      if (params.pltVerbose) {
        pltOpts[0] = minimum;	pltVals[0] = "0";
        pltOpts[1] = maximum;	pltVals[1] = "200";
        plt.printRC(imgCent, 
            "plots/center_holeRemoval" + to_string(imgINFO[icnt].imgNum),
            pltOpts, pltVals);
      }
   
      center[0] = centersCOM[0];
      center[1] = centersCOM[1];
      centfnctr.radProc     = &radProc;
      centfnctr.fxnType     = params.centerFxnType;
      centfnctr.minRadBin   = params.centerMinRadBin;
      centfnctr.centShellW  = params.centerShellWidth;
      centfnctr.img         = &imgCent;
      //centfnctr.verbose     = params.verbose;

      if (params.verbose) std::cout << "\tSearching for center (PowellMin) ... ";
      tools::powellMin<centerfnctr> (centfnctr, center, 
          params.cntrScale, params.cntrMinScale, 
          params.cntrPowellTol, params.cntrFracTol1d);
      if (params.verbose) std::cout << center[0] << "  " << center[1] << std::endl;

      //if ((curDate == "20131102") && (curScan == "LongScan1")) {
      //  center[0] = 569;
      //  center[1] = 490;
      //}
      centerR += center[0];
      centerC += center[1];
      Ncents += 1;

      // Show image center
      if (params.pltCent) {
        TH2F* cImg = plt.plotRC(imgCent, "plots/centeredImage"+to_string(icnt), maximum, "1500");
        for (int ir=(center[0])-5; ir<=center[0]+5; ir++) {
          for (int ic=(center[1])-5; ic<=center[1]+5; ic++) {
            if (sqrt(pow(ir-center[0],2)+pow(ic-center[1],2)) < 5) 
              cImg->SetBinContent(ic, ir, 100000);
          }
        }
        cImg->SetMinimum(-1);
        plt.print2d(cImg, "plots/centeredImage"+to_string(icnt));
        delete cImg;
      }

    }
    centerR /= Ncents;
    centerC /= Ncents;
    if (params.verbose) cout << "INFO: Found center " << centerR 
                  << " "<< centerC << "!\n\n";
  }


  ///////////////////////////////////////
  /////  Processing images in scan  /////
  ///////////////////////////////////////
 
  for (ifl=0; ifl<imgINFO.size(); ifl++) {

    //if (imgINFO[ifl].stagePos != 1542950) continue; 
    centerR = centers[imgINFO[ifl].stagePos][1];
    centerC = centers[imgINFO[ifl].stagePos][0];
    cout<<"centers: "<<centerR<<"  "<<centerC<<endl;

    //cout<<"\n\n\n\nSTARTING"<<endl;
    /////  Check we are in the same run/scan  /////
    if ((imgINFO[ifl].run.compare(curRun) != 0) && (imgINFO[ifl].scan != curScan)) {
      std::cerr << "ERROR new image in run/scan " << imgINFO[ifl].run 
        << "/" << imgINFO[ifl].scan << " instead of " 
        << curRun << "/" << curScan << endl;
      exit(0);
    }

    if (params.verbose) cout << "\tStage position: " 
      + to_string(imgINFO[ifl].stagePos) << endl;

    /////  Check if reference image  /////
    imgIsRef = 0;
    if (params.hasRef) {
      if (imgINFO[ifl].stagePos < params.refStagePos) {
        imgIsRef = 1;
      }
    }

    /////  Image number (ordered)  /////
    imgNum = imgINFO[ifl].imgNum;

    /////  Stage position  /////
    stagePos = imgINFO[ifl].stagePos;

    /////  Image capture time  /////
    timeStamp = imgINFO[ifl].time;

    /////  Throttle  /////
    if (imgINFO[ifl].throttle == -1) {
      throttle = params.throttle;
    }
    else {
      throttle = imgINFO[ifl].throttle;
    }

    /////  Get PV values  /////
    if (params.getPVs) {
      for (auto const & pv : params.pvMap) {
        std::string pvName = pv.first;
        pvVals[pvName] = 0;
        pvValsDer[pvName] = 0;
        int pvInd = (imgINFO[ifl].time - pvStartTimes[pvName])/
                      params.pvSampleTimes - 2;
        for (int i=0; i<params.imgShutterTime/params.pvSampleTimes; i++) {
          pvVals[pvName]    += pvAll[pvName][pvInd-i];
          pvValsDer[pvName] += pvAllDer[pvName][pvInd-i];
        }
        pvVals[pvName] /= params.imgShutterTime/params.pvSampleTimes;
        pvValsDer[pvName] /= params.imgShutterTime/params.pvSampleTimes;
      }
    }
 
    ///////  Load image  ///////
    imgAddr = imgINFO[ifl].path + imgINFO[ifl].fileName;
    if (params.verbose) cout << "INFO: Trying to open " << imgAddr << "\t .....";
    imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
    if (params.verbose) cout << "\tpassed!\n\n";

    if (params.imgMatType.compare("uint16") == 0) {
      if (params.hasLaserBkg) {
        imgOrig = imgProc::getImgVector<uint16_t>(imgMat, 
            Nrows, Nrows/2, Ncols/2, true, NULL,
            params.holeRad, params.holeR, params.holeC);
        imgSubBkg = imgProc::getImgVector<uint16_t>(imgMat, 
            Nrows, Nrows/2, Ncols/2, true, &imgBkg,
            params.holeRad, params.holeR, params.holeC,
            &params.nanMap);
      }
      else {
        imgOrig = imgProc::getImgVector<uint16_t>(imgMat, 
            Nrows, Nrows/2, Ncols/2, true, NULL,
            params.holeRad, params.holeR, params.holeC);
        imgSubBkg = imgProc::getImgVector<uint16_t>(imgMat, 
            Nrows, Nrows/2, Ncols/2, true, &imgBkg,
            params.holeRad, params.holeR, params.holeC,
            &params.nanMap);
      }
    }
    else if (params.imgMatType.compare("uint32") == 0) {
      if (params.hasLaserBkg) {
        imgOrig = imgProc::getImgVector<uint32_t>(imgMat, 
            Nrows, Nrows/2, Ncols/2, true, NULL,
            params.holeRad, params.holeR, params.holeC);
        imgSubBkg = imgProc::getImgVector<uint32_t>(imgMat, 
            Nrows, Nrows/2, Ncols/2, true, &imgBkg,
            params.holeRad, params.holeR, params.holeC,
            &params.nanMap);
      }
      else {
        imgOrig = imgProc::getImgVector<uint32_t>(imgMat, 
            Nrows, Nrows/2, Ncols/2, true, NULL,
            params.holeRad, params.holeR, params.holeC);
        imgSubBkg = imgProc::getImgVector<uint32_t>(imgMat, 
            Nrows, Nrows/2, Ncols/2, true, &imgBkg,
            params.holeRad, params.holeR, params.holeC,
            &params.nanMap);
      }
    }
    else {
      cerr << "ERROR: Do not recognize imgMatType = " 
          + params.imgMatType + "!!!\n";
      exit(0);
    }
   

    // Remove Xray hits
    imgProc::removeXrayHits(
        &imgOrig, 
        params.XrayHighCut, params.XrayLowCut, 
        params.XraySTDcut, params.XrayWindow);
    imgProc::removeXrayHits(
        &imgSubBkg, 
        params.XrayHighCut, params.XrayLowCut, 
        params.XraySTDcut, params.XrayWindow);



    //cout<<"444"<<endl;

    /////  Remove pixel outliers  /////
    /*
    radProc.removeOutliers(imgOrig,  
        centerR, centerC, params.imgEdgeBuffer, 
        params.NradAzmBins, params.shellWidth, params.Npoly,
        //imgOrig.size()/2, imgOrig.size()/2, params.shellWidth,
        params.stdIncludeLeft, params.distSTDratioLeft,
        params.stdCutLeft, params.meanBinSize,
        params.stdIncludeRight, params.distSTDratioRight,
        params.stdChangeRatio, params.stdCutRight,
        imgINFO[ifl].stagePos, params.outlierMapSTDcut,
        false, (false || params.outlierVerbose), NULL);
        */

    //cout<<"555"<<endl;
    outlierImage = radProc.removeOutliersSimple(
        imgSubBkg, centerR, centerC, params.imgEdgeBuffer,
        params.NradAzmBins, params.shellWidth, params.Npoly,
        params.outlierSTDcut, imgINFO[ifl].stagePos, params.outlierMapSTDcut,
        true, (false || params.outlierVerbose), NULL);//, &plt);

    /*
    outlierImage = radProc.removeOutliers(imgSubBkg,  
        centerR, centerC, params.imgEdgeBuffer,
        params.NradAzmBins, params.shellWidth, params.Npoly,
        //imgOrig.size()/2, imgOrig.size()/2, params.shellWidth,
        params.stdIncludeLeft, params.distSTDratioLeft,
        params.stdCutLeft, params.meanBinSize,
        params.stdIncludeRight, params.distSTDratioRight,
        params.stdChangeRatio, params.stdCutRight,
        imgINFO[ifl].stagePos, params.outlierMapSTDcut,
        true, (false || params.outlierVerbose), NULL);//, &plt);
        */
    
    if (params.pltVerbose) {
      //plt.printRC(outlierImage, "outlierSTD_" + imgINFO[ifl].run + "-" + to_string(imgINFO[ifl].scan) + "-" + to_string(imgINFO[ifl].stagePos));
      save::saveDat<double>(outlierImage, 
            "./results/outlierSTD-" + runName
                + "-" + to_string(imgINFO[ifl].scan)
                + "-" + to_string(imgINFO[ifl].stagePos) + ".dat");
    }

    ///  Find large clusters corresponding to laser spots  ///
    //std::vector< std::pair<uint, uint> > removePairs;
    /*
    std::vector< std::pair<uint, uint> > removePairs = imgProc::findClusters(
            outlierImage, centerR, centerC, params.holeRad,
            params.outlierCoreValThresh, params.outlierCoreRad, 
            params.outlierMinClusterSize, params.outlierMinPixelSize, 
            params.outlierMinDensity, params.outlierShapeVarCut, 
            params.outlierShapeEdgeCut, params.outlierClusterRad,
            params.outlierBorderValThresh, params.outlierBorderDistLimit, 
            params.outlierBorderRad, params.outlierPadRad, 
            params.outlierrMinScale, params.outlierrMaxScale,
            params.outliercMinScale, params.outliercMaxScale,
            outlierBkg, NULL); //&plt);

    ///  Remove laser background from image by setting pixels = NANVAL  ///
    for (uint ip=0; ip<removePairs.size(); ip++) {
      imgSubBkg[removePairs[ip].first][removePairs[ip].second] = NANVAL;
    }

    if (params.pltVerbose) {
      //plt.printRC(outlierBkg, "outlierBKG_" + imgINFO[ifl].run + "-" + to_string(imgINFO[ifl].scan) + "-" + to_string(imgINFO[ifl].stagePos));
      save::saveDat<double>(outlierBkg, 
            "./results/outlierBackground-" + runName
                + "-" + to_string(imgINFO[ifl].scan)
                + "-" + to_string(imgINFO[ifl].stagePos) + ".dat");
      for (int ir=0; ir<(int)imgSubBkg.size(); ir ++) {
        for (int ic=0; ic<(int)imgSubBkg[ir].size(); ic ++) {
          if (sqrt(pow(ir-params.holeR,2) + pow(ic-params.holeC,2)) < params.holeRad) continue;
          if (imgSubBkg[ir][ic] == NANVAL) imgSubBkg[ir][ic] = params.hotPixel;
        }
      }

      save::saveDat<double>(imgSubBkg, 
            "./results/image-" + runName
                + "-" + to_string(imgINFO[ifl].scan)
                + "-" + to_string(imgINFO[ifl].stagePos) 
                + "[" + to_string(imgSubBkg.size())
                + "," + to_string(imgSubBkg.size()) + "].dat");
    }

            */


    ///// Remove readout noise  /////
    /*
    if (params.verbose) std::cout << "INFO: Readout noise subtraction.\n";
    if (params.hasLaserBkg) {
      readoutNoise = imgProc::removeAvgReadOutNoise(imgSubBkg, centerR, centerC, 
                        0.94*params.NradAzmBins, params.NradAzmBins, 
                        params.imgEdgeBuffer, &params.nanMap);
    }
    else {
      readoutNoise = imgProc::removeAvgReadOutNoise(imgSubBkg, centerR, centerC, 
                        0.94*params.NradAzmBins, params.NradAzmBins,
                        params.imgEdgeBuffer);

      // Only use this when interested in reference images
      //    background fit fails when subtracting readout
      //    noise on the phosphor screen. Comment previous line.
      //readoutNoise = imgProc::removeReadOutNoise(imgSubBkg);
    }
    */


    //////////////////////////////////////////////////////
    /////  Finding and subtracting laser background  /////
    //////////////////////////////////////////////////////
    // We cluster background spots based on the assymetry of 
    //    the image. 

    /*
    if (params.hasLaserBkg && params.laserClusterRemoval) {
      if (params.verbose) std::cout << "INFO: Laser background removal.\n";
      // Fill in hole based on symmetry
      std::vector< std::vector<double> > imgLaser = imgSubBkg;
      for (int ir=centerR-params.holeRad-3; 
          ir<=centerR+params.holeRad+3; ir++) {
        for (int ic=centerC-params.holeRad-3; 
            ic<=centerC+params.holeRad+3; ic++) {
          if (std::pow(ir-centerR,2) + std::pow(ic-centerC,2) 
              < std::pow(params.holeRad, 2)) {
            imgLaser[ir][ic] = imgLaser[(imgSize-1)-ir][(imgSize-1)-ic];
          }
        }
      }


      ///  Finding asymmetric parts of the image  ///
      std::vector< std::vector<double> > oddImgReal = imgProc::asymmetrize(imgLaser, 
          centerR, centerC, imgSize, imgSize, 
          oddImgImgn, fftFref, fftIn, fftBref, fftOut);

      // Building map of ratio of noise/"signal" = asymmetric/symmetric
      for (int ir=0; ir<imgSize; ir++) { 
        for (int ic=0; ic<imgSize; ic++) {
          symImg[ir][ic] = imgLaser[ir][ic] - oddImgReal[ir][ic];
          stdRatioImg[ir][ic] = oddImgReal[ir][ic]/std::max(std::abs(symImg[ir][ic]),0.01)
                                  /sqrt(pow(ir - centerR, 2) + pow(ic - centerC, 2));
        }
      }

      ///  Find large clusters corresponding to laser spots  ///
      std::vector< std::pair<uint, uint> > removePairs = imgProc::findClusters(
              stdRatioImg, centerR, centerC, params.holeRad*1.2,
              params.coreValThresh, params.coreRad, 
              params.minClusterSize, params.minPixelSize, 
              params.minDensity, params.clusterRad,
              params.borderValThresh, 1, params.borderRad, params.padRad, imgLaserBkg);


      ///  Remove laser background from image by setting pixels = NANVAL  ///
      for (uint ip=0; ip<removePairs.size(); ip++) {
        imgSubBkg[removePairs[ip].first][removePairs[ip].second] = NANVAL;
      }
      ///  Remove laser background from image by setting pixels = NANVAL  ///
      //for (ir=0; ir<imgSize; ir++) {
      //  for (ic=0; ic<imgSize; ic++) {
      //    imgSubBkg[ir][ic] = imgOrig[ir][ic];
      //  }
      //}

      //for (ir=imgSize/2-40; ir<imgSize/2+40; ir++) {
      //  for (ic=imgSize/2-120; ic<imgSize/2; ic++) {
      //    imgSubBkg[ir][ic] = NANVAL;
      //  }
      //}
      //for (ir=imgSize/2-200; ir<imgSize/2+50; ir++) {
      //  for (ic=imgSize-225; ic<imgSize-115; ic++) {
      //    imgSubBkg[ir][ic] = NANVAL;
      //  }
      //}
      //for (ir=imgSize/2-60; ir<imgSize/2+35; ir++) {
      //  for (ic=imgSize-75; ic<imgSize; ic++) {
      //    imgSubBkg[ir][ic] = NANVAL;
      //  }
      //}
      
      if (params.pltVerbose) {
        
        std::vector< std::vector<double> > pltStdRat(imgSize);
        for (int ir=0; ir<imgSize; ir++) {
          pltStdRat[ir].resize(imgSize);
          for (int ic=0; ic<imgSize; ic++) {
            if (std::pow(ir-centerR,2) + std::pow(ic-centerC,2) 
                < std::pow(params.holeRad, 2)) {
              pltStdRat[ir][ic] = 0;
            }
            else {
              pltStdRat[ir][ic] = stdRatioImg[ir][ic];
            }
          }
        }
 
        delete plt.printRC(oddImgReal, 
            "plots/oddImgReal_"+curRun+"_"+to_string(curScan)+"_"+to_string(stagePos));//, 
            //pltOpts, pltVals);

        pltOpts[0] = minimum;	pltVals[0] = "0";
        pltOpts[1] = maximum;	pltVals[1] = "5e6";
        delete plt.printRC(imgLaserBkg, 
            "plots/imgLaserBkg_"+curRun+"_"+to_string(curScan)+"_"+to_string(stagePos));

        pltOpts[0] = minimum;	pltVals[0] = "0"; //"1e4";
        pltOpts[1] = maximum;	pltVals[1] = to_string(params.coreValThresh); //"1e6";
        //pltOpts.push_back(logz);  pltVals.push_back("");
        delete plt.printRC(pltStdRat, 
            "plots/imgStdRat_"+curRun+"_"+to_string(curScan)+"_"+to_string(stagePos),
            pltOpts, pltVals);

        for (int ir=imgSize/2-40; ir<imgSize/2+40; ir++) {
          for (int ic=imgSize/2-120; ic<imgSize/2; ic++) {
            pltStdRat[ir][ic] = 0;
          }
        }
        for (int ir=imgSize/2-200; ir<imgSize/2+50; ir++) {
          for (int ic=imgSize-225; ic<imgSize-115; ic++) {
            pltStdRat[ir][ic] = 0;
          }
        }
        for (int ir=imgSize/2-60; ir<imgSize/2+35; ir++) {
          for (int ic=imgSize-75; ic<imgSize; ic++) {
            pltStdRat[ir][ic] = 0;
          }
        }
        delete plt.printRC(pltStdRat, 
            "plots/imgStdRatRemove_"+curRun+"_"+to_string(curScan)+"_"+to_string(stagePos),
            pltOpts, pltVals);
      }
    }


    ////////////////////////////////////////
    /////  Plot results for debugging  /////
    ////////////////////////////////////////

    if (params.pltVerbose) {
      std::vector< std::vector<double> > pltOrig(Nrows);
      std::vector< std::vector<double> > pltSubBkg(Nrows);
      for (int ir=0; ir<Nrows; ir++) {
        pltOrig[ir].resize(Ncols);
        pltSubBkg[ir].resize(Ncols);
        for (int ic=0; ic<Ncols; ic++) {
          pltOrig[ir][ic] = imgOrig[ir][ic];
          pltSubBkg[ir][ic] = imgSubBkg[ir][ic];
        }
      }
      delete plt.printRC(pltOrig, 
          "plots/imgOrig_"+curRun+"_"+to_string(curScan)+"_"+to_string(stagePos));
      delete plt.printRC(pltSubBkg, 
          "plots/imgSub_"+curRun+"_"+to_string(curScan)+"_"+to_string(stagePos));
    }
    */


    //////////////////////////////////////////////////////////
    /////  Filling image variables and filling the tree  /////
    //////////////////////////////////////////////////////////


    //// TODO!!!!!!!!!!!!! MUST CENTER IMAGE BEFORE FITTING
    //////  Legendre Fit  //////
    assert(imgSize%5 == 0);
    assert(imgSize%5 == 0);
    const int lgFit_Rows = imgSize/5;
    const int lgFit_Cols = imgSize/5;
    // Check if g matrix already exists, else make new one
    string matrix_folder = "/reg/neh/home/khegazy/analysis/legendreFitMatrices/";
    string matrix_fileName = "gMatrix_row-" + to_string(lgFit_Rows)
        + "_col-" + to_string(lgFit_Cols) + "_Nrad-" + to_string(NradLegBins)
        + "_Nlg-" + to_string(Nlegendres) + ".dat";
    if (params.verbose) std::cout << "INFO: Looking for " + matrix_fileName << endl;
    if (access((matrix_folder + matrix_fileName).c_str(), F_OK) == -1) {
      cout << "INFO: Making new g matrix\n";
      system(("python " + matrix_folder + "makeLgMatrix.py --NradLegBins="
            + to_string(NradLegBins) + " --Ncols=" + to_string(lgFit_Cols)
            + " --Nrows=" + to_string(lgFit_Rows)
            + " --Nlg=" + to_string(Nlegendres)).c_str());
    }

    // Import the g matrix
    const int NgMat = lgFit_Rows*lgFit_Cols*NradLegBins*Nlegendres;
    const int Npix = lgFit_Rows*lgFit_Cols;
    double* gInp = new double[NgMat];
    FILE* inpFile = fopen((matrix_folder + matrix_fileName).c_str(), "rb");
    fread(gInp, sizeof(double), NgMat, inpFile);
    Eigen::Map< Eigen::Matrix<double, Npix, 
        Nlegendres*NradLegBins, Eigen::RowMajor> > gOrig(gInp); 

    clock_t begin = clock();
    //legCoeffs = imgProc::legendreFit(imgSubBkg, 5, Nlegendres, 
    //    NradLegBins, lgFit_Rows, lgFit_Cols, NANVAL, gOrig);
    

    clock_t end = clock();
    if (params.verbose) cout<<"TIME: "<<double(end - begin) / CLOCKS_PER_SEC<<endl;


    /////  Azimuthal average  /////

    std::fill(rawAzmAvg.begin(), rawAzmAvg.end(), 0);
    std::fill(azmCounts.begin(), azmCounts.end(), 0);
    for (int ir=0; ir<(int)imgSubBkg.size(); ir++) {
      if (ir < params.imgEdgeBuffer || imgSubBkg.size() - ir < params.imgEdgeBuffer) continue;
      for (int ic=0; ic<(int)imgSubBkg[ir].size(); ic++) {
        if (ic < params.imgEdgeBuffer || imgSubBkg[ir].size() - ic < params.imgEdgeBuffer) continue;

        if (imgSubBkg[ir][ic] != NANVAL) {
          if (params.hasLaserBkg) {
            if (params.nanMap[ir][ic] == NANVAL) {
              //continue; //FIX ME comment me out after comparing with thomas
            }
          }

          radInd = std::round(
                      std::sqrt(
                        std::pow(ir-centerR,2) + std::pow(ic-centerC,2)));

          if (radInd < params.NradAzmBins) {
            rawAzmAvg[radInd] += imgSubBkg[ir][ic];
            azmCounts[radInd] += 1;
          }
        }
        else {
          //cout<<"skipping nan"<<endl;
        }
      }
    }

    for (uint ir=0; ir<rawAzmAvg.size(); ir++) {
      if (azmCounts[ir] != 0) {
        rawAzmAvg[ir] /= azmCounts[ir];
      }
      else rawAzmAvg[ir] = NANVAL;
    }

    rawAzmAvg.resize(params.NradAzmBins, NANVAL);

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      azmAvg[iq] = rawAzmAvg[iq];
      //if (iq>=100 &&iq <=110) {
      //  cout<<"azmAvg "<<iq<<" : "<<azmAvg[iq]<<endl;
      //}
    }
    //save::saveDat<double>(azmAvg, "testcompareazmavg.dat");
      //plt.print1d(azmAvg, "azimuthalAvg_" + to_string(imgINFO[ifl].imgNum), logy, "true");

    /*
    /////  Image norm  /////
    imgNorm = 0;
    count = 0;
    for (int i=imgNormBinMin; i<imgNormBinMax; i++) {
      if (rawAzmAvg[i] == NANVAL) continue;
      imgNorm += rawAzmAvg[i];
      count++;
    }
    imgNorm /= count;

    /////  Build Reference AzmAvg  /////
    if (imgIsRef) {
      for (int iq=0; iq<params.NradAzmBins; iq++) {
        refAzmAvg[iq] = (refAzmAvg[iq]*Nref + rawAzmAvg[iq]/imgNorm)/(Nref + 1);
      }
      Nref++;
    }

    /////  Filtering  /////

    std::vector<double> tstR(filtFFToutSize);
    if (!padRange) {
      while (rawAzmAvg[padRange] == NANVAL) {
        padRange++;
      }
    }
    curPadRange = padRange;
    while (fabs(1 - rawAzmAvg[curPadRange]/rawAzmAvg[curPadRange+1]) > 0.2) {
      curPadRange++;
    }

    qSpace[0] = 0;
    for (int iq=1; iq<params.NradAzmBins; iq++) {
      qSpace[iq] = (rawAzmAvg[iq]/imgNorm)/atmDiff[iq];
      if (iq < params.suppressBins) {
        if (iq < curPadRange && iq > 0) {
          qSpace[iq] = rawAzmAvg[curPadRange]*(pow(curPadRange,3)/(atmDiff[iq]*imgNorm))/pow(iq,3);
            //params.padMaxHeight*sin(PI/2*iq/padRange);
        }
        qSpace[iq] *= pow(sin((PI/2)*iq/params.suppressBins),6);
      }
    }

    if (params.pltVerbose) {
      std::vector<double> tst(params.NradAzmBins);
      for (int iq=0; iq<params.NradAzmBins; iq++) {
        tst[iq] = qSpace[iq];
      }
      pltVals[0] = "-2";
      pltVals[1] = "2";
      delete plt.print1d(tst, "./plots/filtInput_" + to_string(stagePos),pltOpts,pltVals);
    }

    fftw_execute(filtFFTf);

    for (int ir=0; ir<filtFFToutSize; ir++) {
      rSpace[ir][0] *= bandPassFilter[ir]/sqrt(params.NradAzmBins);
      rSpace[ir][1] *= bandPassFilter[ir]/sqrt(params.NradAzmBins);
    }

    fftw_execute(filtFFTb);

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      azmAvg[iq] = qSpace[iq]*atmDiff[iq]
                          /sqrt(params.NradAzmBins);
      if (iq < params.suppressBins) {
        if (iq < padRange) {
          azmAvg[iq] = NANVAL;
        }
        else {
          azmAvg[iq] /= pow(sin((PI/2)*iq/params.suppressBins),6);
        }
      }
    }

    if (true || params.pltVerbose) {
      std::vector<double> tst1(params.NradAzmBins);
      std::vector<double> tst2(params.NradAzmBins);
      vector<TH1*> plts(2);
      for (int iq=0; iq<params.NradAzmBins; iq++) {
        tst1[iq] = rawAzmAvg[iq]/(imgNorm*atmDiff[iq]);
        tst2[iq] = azmAvg[iq]/atmDiff[iq];
      }
      pltVals[0] = "-0.2";
      pltVals[1] = "0.2";
      plts[0] = plt.plot1d(tst1, "blah_" + to_string(stagePos),pltOpts, pltVals);
      plts[1] = plt.print1d(tst2, "./plots/filtOutp_" + to_string(stagePos),pltOpts, pltVals);
      plts[1]->SetLineColor(4);
      plts[0]->SetMaximum(40);
      plts[0]->SetMinimum(0);
      plt.print1d(plts, "./plots/filtCompare_" + to_string(stagePos));
      delete plts[0];
      delete plts[1];
    }
    */

    /////  Plotting results  /////
    if (true || params.pltVerbose) {
      std::vector<double> test(NradLegBins);
      //for (int i=0; i<params.Nlegendres; i++) {
      //  for (int j=0; j<NradLegBins; j++) {
      //    test[j] = legCoeffs[i*NradLegBins + j];
      //  }
      //  plt.print1d(test, "plots/testLeg" + to_string(i) + "_" + to_string(imgINFO[ifl].imgNum));
      //}

      plt.print1d(rawAzmAvg, "plots/azimuthalAvg_" + to_string(imgINFO[ifl].imgNum));
    }
    
    // Save image and info to tree
    tree->Fill();
  }


  tree->Write();
  file->Close();
  cout<<endl<<endl<<endl;

  // Release fftw memory
  fftw_destroy_plan(filtFFTf);
  fftw_destroy_plan(filtFFTb);
  fftw_free(qSpace);
  fftw_free(rSpace);
  fftw_destroy_plan(fftFref);
  fftw_destroy_plan(fftBref);
  fftw_free(fftIn);
  fftw_free(fftOut);

return 1;
}

                                        
