#include "preProcessing.h"
#include "/reg/neh/home/khegazy/baseTools/tools/parameters.h"

using namespace std;


int main(int argc, char* argv[]) {


  // Make runLists
  bool doRunLists = false;

  bool localTesting = false;


  if (argc != 2) {
    std::cerr << "ERROR: Must run program a list of files" 
      << " ./preProcessing.exe runList.txt!!!\n";
    exit(1);
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
      std::cerr << "ERROR: Cannot find runName!!!\n";
      exit(1);
    }
    runName = runListName.substr(iPos, 
        runListName.find("_Scan-") - iPos);
  }

  parameterClass params(runName);
  PLOTclass plt;
  std::string imgCentCodeDir = "/reg/neh/home/khegazy/baseTools/UEDanalysis/preProcessing/";
  std::string dirName = "/reg/ued/ana/scratch/nitroBenzene/polarLineOutTest/";
    std::vector<double> pLO;


  //plt.printRC(params.nanMap, "testingNANmap");

  const int Nlegendres = 1;
  const int NradLegBins = 50;

  // Indices
  //const int imgSize = 935;
  const int imgSize = 895;

  if (Nlegendres != params.Nlegendres) {
    std::cerr 
      << "ERROR: parameter Nlegendres does not match with parameter class!!!\n";
    exit(1);
  }
  if (NradLegBins != params.NradLegBins) {
    std::cerr 
      << "ERROR: parameter NradLegBins does not match with parameter class!!!\n";
    exit(1);
  }
  if (imgSize != params.imgSize) {
    std::cerr 
      << "ERROR: parameter imgSize does not match with parameter class!!!\n";
    exit(1);
  }


   // Make sure the image has odd number of bins
  if (!(imgSize%2) || !(imgSize%2)) {
    std::cerr 
      << "ERROR: imgSize and imgSize must be an odd number!!!\n";
    exit(1);
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

  fftw_complex* fftIn = 
      (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*imgSize*imgSize);
  fftw_complex* fftOut = 
      (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*imgSize*imgSize);
  fftw_plan fftFref = fftw_plan_dft_2d(
      imgSize, imgSize, 
      fftIn, fftOut, 
      FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan fftBref = fftw_plan_dft_2d(
      imgSize, imgSize, 
      fftIn, fftOut, 
      FFTW_BACKWARD, FFTW_MEASURE);


  std::vector< std::vector<double> > radBins;

  imgProc::centerfnctr centfnctr;
  std::vector<double> center(2);
  std::vector<long int> inpCenter(2);
  int centerC, centerR;
  float centerC_float, centerR_float;
  float centerRstdRatio, centerCstdRatio;
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
  fftw_complex* rSpace = (fftw_complex*) fftw_malloc(
      sizeof(fftw_complex)*filtFFToutSize);
  fftw_plan filtFFTf = fftw_plan_dft_r2c_1d(
      params.NradAzmBins, 
      qSpace, rSpace, 
      FFTW_MEASURE);
  fftw_plan filtFFTb = fftw_plan_dft_c2r_1d(
      params.NradAzmBins, 
      rSpace, qSpace, 
      FFTW_MEASURE);

  Eigen::Matrix<double, Eigen::Dynamic, 2> X;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Y;
  Eigen::VectorXd w, fit;


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

  /////  Testing/Debugging variables  /////
  TH1F** xRayHitHistos = NULL;
  if (params.xRayHitDist) {
    xRayHitHistos = new TH1F*[2];
    xRayHitHistos[0] = new TH1F("xRayHits_high", "xRayHits_high", 500, 0, params.XrayHighCut);
    xRayHitHistos[1] = new TH1F("xRayHits_low", "xRayHits_low", 500, 0, params.XrayHighCut);
  }

  TH1F** radPixelHistos = NULL;
  if (params.plotRadPixDist) {
    radPixelHistos = new TH1F*[params.NradAzmBins];
    for (int i=0; i<params.NradAzmBins; i++) {
      radPixelHistos[i] = new TH1F(
            ("radPixDist_rad-" + to_string(i)).c_str(),
            ("radPixDist_rad-" + to_string(i)).c_str(),
            2000, -2000, 2000);
    }
  }



  /////  Calculate sMs normalization  /////
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


  ////////////////////////////////////
  /////  Get I0 Image Addresses  /////
  ////////////////////////////////////

  std::map<int, std::string> I0fileNames;
  std::map<int, std::vector< std::vector<double> > > I0refImgs;
  std::string lineI0;
  std::string I0runListName = "./runLists/runList_Run-"
      + runName + "_I0_Scan-" 
      + to_string(imgINFO[0].scan) + ".txt"; 
  ifstream I0runList(I0runListName);
  if (!I0runList.is_open()) {
    std::cerr << "ERROR: Cannot open file " + I0runListName + "\n";
    exit(1);
  }
  while (getline(I0runList, lineI0)) {
    ipos = lineI0.find("delay");
    ipos += 10;
    int I0stagePos = stoi(lineI0.substr(ipos, 3) + lineI0.substr(ipos+4, 4));
    I0fileNames[I0stagePos] = lineI0;
  }
  I0runList.close();


  /////////////////////////////////////
  /////  Setting image variables  /////
  /////////////////////////////////////

  string imgAddr = I0fileNames[imgINFO[0].stagePos];
  cv::Mat imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
  int Nrows = imgMat.rows;
  int Ncols = imgMat.cols;
  bool normI0ref = false;
  std::vector< std::vector<double> > imgI0(Nrows);
  std::vector< std::vector<double> > imgI0ref(Nrows);
  std::vector< std::vector<double> > imgI0refCount(Nrows);
  for (int ir=0; ir<Nrows; ir++) {
    imgI0[ir].resize(Ncols, 0);
    imgI0ref[ir].resize(Ncols, 0);
    imgI0refCount[ir].resize(Ncols, 0);
  }

  imgAddr = imgINFO[0].path + imgINFO[0].fileName;
  imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
  cv::Mat imgSmooth;
  Nrows = imgMat.rows;
  Ncols = imgMat.cols;
  std::vector<double> imgRadSTD(params.meanInds.size());
  //std::vector<double> radPixMeans(params.NradAzmBins);
  std::vector< std::vector<double> > imgRadSTDs(imgINFO.size());
  std::vector< std::vector<double> > imgOrig(Nrows);
  std::vector< std::vector<double> > imgSubBkg(Nrows);
  std::vector< std::vector<double> > imgTemp(Nrows);
  std::vector< std::vector<double> > imgLaserBkg(Nrows);
  std::vector< std::vector<double> > imgBkg(Nrows);
  std::vector< std::vector<double> > imgCent;
  //std::vector< std::vector<double> > radPixDist(params.NradAzmBins);


  float imgNorm;
  for (int ir=0; ir<Nrows; ir++) {
    imgOrig[ir].resize(Ncols, 0);
    imgSubBkg[ir].resize(Ncols, 0);
    imgTemp[ir].resize(Ncols, 0);
    imgLaserBkg[ir].resize(Ncols, 0);
    imgBkg[ir].resize(1024, 0);
  }
  for (uint i=0; i<imgINFO.size(); i++) {
    imgRadSTDs[i].resize(params.meanInds.size());
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

      imgSmooth = imgMat;
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
        std::cerr << "ERROR: Do not recognize imgMatType = " 
            + params.imgMatType + "!!!\n";
        exit(1);
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

          double mean = std::accumulate(
              collection.begin(), 
              collection.end(), 0);
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
    }
  }



  /// For testing centers
  std::map<int, vector<float> > centersComp;
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

      centersComp[stoi(pos)].push_back(stof(sss));
    }
  }

    /////  Find Center  /////
  //imgProc::centerfnctr centfnctr;
  double meanR, meanC;
  double stdR, stdC;
  double centerValSTD = 2;
  std::vector<double> resultsR(params.meanInds.size());
  std::vector<double> resultsC(params.meanInds.size());
  std::vector< std::vector<int> > centers(imgINFO.size());
  std::vector< std::vector<double> > allMedVals(imgINFO.size());
  std::vector< std::vector< std::vector<int> > > allIndsR(imgINFO.size());
  std::vector< std::vector< std::vector<int> > > allIndsC(imgINFO.size());
  std::vector< std::vector< std::vector<double> > > allCentVals(imgINFO.size());
  for (ifl=0; ifl<imgINFO.size(); ifl++) {
    centers[ifl].resize(2);
    imgAddr = imgINFO[ifl].path + imgINFO[ifl].fileName;
    imgMat = cv::imread(imgAddr.c_str(), CV_LOAD_IMAGE_ANYDEPTH); 

    Nrows = imgMat.rows;
    Ncols = imgMat.cols;
    if (params.imgMatType.compare("uint16") == 0) {
      imgCent = imgProc::getImgVector<uint16_t>(imgMat, //imgSmooth, 
                  Nrows, Nrows/2, Ncols/2, true, &imgBkg,
                  params.holeRad, params.holeR, params.holeC,
                  &params.nanMap);
    }  
    else if (params.imgMatType.compare("uint32") == 0) {
      imgCent = imgProc::getImgVector<uint32_t>(imgMat, //imgSmooth, 
                  Nrows, Nrows/2, Ncols/2, true, &imgBkg,
                  params.holeRad, params.holeR, params.holeC,
                  &params.nanMap);
    }  
    else {
      std::cerr << "ERROR: Do not recognize imgMatType = " 
          + params.imgMatType + "!!!\n";
      exit(1);
    }

    /*
    pltOpts[0] = minimum; pltOpts[1] = maximum;
    pltVals[0] = "1"; pltVals[1] = "1000";
    plt.printRC(imgCent,"testBlock",pltOpts, pltVals);
    exit(0);
    */
    /*
    // Remove Xray hits
    imgProc::removeXrayHits(
        &imgCent, 
        params.XrayHighCut, params.XrayLowCut, 
        params.XraySTDcut, params.XrayWindow);

    for (uint ir=0; ir<Nrows; ir++) {
      for (uint ic=0; ic<Ncols; ic++) {
        imgMat.at<imgMat.type()>(ir, ic) = imgCent[ir][ic];
      }
    }
    */

    int GSsize = (int)(params.sigma*5);
    GSsize += 1 - GSsize%2;
    cv::Mat imgSmooth(Nrows, Ncols, imgMat.type());
    cv::medianBlur(imgMat, imgSmooth, 5);
    cv::medianBlur(imgSmooth, imgMat, 5);
    cv::GaussianBlur(
        imgMat, imgSmooth, 
        cvSize(GSsize, GSsize), 
        params.sigma, params.sigma);

    if (params.imgMatType.compare("uint16") == 0) {
      imgCent = imgProc::getImgVector<uint16_t>(imgSmooth, 
                  Nrows, Nrows/2, Ncols/2, true, NULL,
                  params.holeRad+5, params.holeR, params.holeC,
                  &params.nanMap);
    }  
    else if (params.imgMatType.compare("uint32") == 0) {
      imgCent = imgProc::getImgVector<uint32_t>(imgSmooth, 
                  Nrows, Nrows/2, Ncols/2, true, NULL,
                  params.holeRad+5, params.holeR, params.holeC,
                  &params.nanMap);
    }  
    else {
      std::cerr << "ERROR: Do not recognize imgMatType = " 
          + params.imgMatType + "!!!\n";
      exit(1);
    }

    std::vector<int> centersCOM = imgProc::centerSearchCOM(imgCent,
        params.blockCentR, params.blockCentC, params.minRad, params.maxRad,
        params.meanInd, params.COMstdScale, true, NULL); //&plt);

    cout<<"COM: "<<centersCOM[0]<<"  "<<centersCOM[1]<<endl;
    centerR_float = 0;
    centerC_float = 0;
    Ncents  = 0;
    double medVal;
    std::vector<double> Rs, centVals;
    std::vector<double> circVals(4);
    std::vector<double> orderedInds(4);
    std::vector<int> indsR, indsC, centerRemoveInds;
    allIndsR[ifl].resize(params.meanInds.size());
    allIndsC[ifl].resize(params.meanInds.size());
    allCentVals[ifl].resize(params.meanInds.size());
    allMedVals[ifl].resize(params.meanInds.size());
    for (int i=0; i<params.meanInds.size(); i++) {
      circVals[0] = imgCent[centersCOM[0]+params.meanInds[i]][centersCOM[1]];
      circVals[1] = imgCent[centersCOM[0]-params.meanInds[i]][centersCOM[1]];
      circVals[2] = imgCent[centersCOM[0]][centersCOM[1]+params.meanInds[i]];
      circVals[3] = imgCent[centersCOM[0]][centersCOM[1]-params.meanInds[i]];
      iota(orderedInds.begin(), orderedInds.end(), 0);
      sort(orderedInds.begin(), orderedInds.end(),
          [&circVals](int i1, int i2)
          {return circVals[i1] < circVals[i2];});
      medVal = (circVals[orderedInds[1]] + circVals[orderedInds[2]])/2.;

      // Find indices with values close to medVal
      indsR.clear();
      indsC.clear();
      centVals.clear();
      for (int ir=0; ir<Nrows; ir++) {
        for (int ic=0; ic<Ncols; ic++) {
          if (fabs(imgCent[ir][ic] - medVal) <= 4*centerValSTD) {
            indsR.push_back(ir);
            indsC.push_back(ic);
            centVals.push_back(imgCent[ir][ic]);
          }
        }
      }

      // Remove indices far from cluster
      for (int j=0; j<3; j++) {
        meanR = 0;
        Rs.resize(indsR.size());
        for (uint k=0; k<indsR.size(); k++) {
          Rs[k] = std::sqrt(
                      std::pow(indsR[k] - centerR_float, 2)
                      + std::pow(indsC[k] - centerC_float, 2));
          meanR += Rs[k];
        }
        meanR /= indsR.size();

        stdR = 0;
        for (uint k=0; k<indsR.size(); k++) {
          stdR += std::pow(meanR - Rs[k], 2);
        }
        stdR = std::sqrt(stdR/indsR.size());

        centerRemoveInds.clear();
        for (uint k=0; k<indsR.size(); k++) {
          if (std::fabs(Rs[k] - meanR) > 4*stdR) {
            centerRemoveInds.push_back(k);
          }
        }

        for (int k=centerRemoveInds.size()-1; k>=0; k--) {
          indsR.erase(indsR.begin()+centerRemoveInds[k]);
          indsC.erase(indsC.begin()+centerRemoveInds[k]);
          centVals.erase(centVals.begin()+centerRemoveInds[k]);
        }
      }

      imgRadSTDs[ifl][i] = stdR;

      allIndsR[ifl][i].resize(indsR.size());
      allIndsC[ifl][i].resize(indsC.size());
      allCentVals[ifl][i].resize(indsR.size());
      for (uint k=0; k<indsR.size(); k++) {
        allIndsR[ifl][i][k] = indsR[k];
        allIndsC[ifl][i][k] = indsC[k];
        allCentVals[ifl][i][k] = centVals[k];
      }
      allMedVals[ifl][i] = medVal;

      center[0] = centersCOM[0];
      center[1] = centersCOM[1];
      centfnctr.radProc     = &radProc;
      centfnctr.fxnType     = params.centerFxnType;
      centfnctr.meanVal     = medVal;
      centfnctr.valSTD      = centerValSTD;
      centfnctr.vals        = &centVals;
      centfnctr.indsR       = &indsR;
      centfnctr.indsC       = &indsC;
      //centfnctr.verbose     = params.verbose;

      if (params.verbose) std::cout << "\tSearching for center (PowellMin) ... ";
      tools::powellMin<imgProc::centerfnctr> (centfnctr, center,
          params.cntrScale, params.cntrMinScale,
          params.cntrPowellTol, params.cntrFracTol1d);
      if (params.verbose) std::cout << center[0] << "  " << center[1] << std::endl;

      double minLoss = 1e90;
      int minR, minC;
      for (int ir=-5; ir<=5; ir++) {
        for (int ic=-5; ic<=5; ic++) {
          std::vector<double> curCent{center[0]+ir, center[1]+ic};

          if (minLoss > centfnctr(curCent)) {
            minLoss = centfnctr(curCent);
            minR = curCent[0];
            minC = curCent[1];
          }
        }
      }

      resultsR[i] = minR;
      resultsC[i] = minC;
    }

    meanR = std::accumulate(resultsR.begin(), resultsR.end(), 0);
    meanR /= (float)resultsR.size();
    meanC = std::accumulate(resultsC.begin(), resultsC.end(), 0);
    meanC /= (float)resultsC.size();

    stdR = 0;
    stdC = 0;
    for (uint k=0; k<params.meanInds.size(); k++) {
      stdR += std::pow(meanR-resultsR[k], 2);
      stdC += std::pow(meanC-resultsC[k], 2);
    }
    stdR = std::sqrt(stdR/params.meanInds.size());
    stdC = std::sqrt(stdC/params.meanInds.size());

    centerR_float = 0;
    centerC_float = 0;
    Ncents  = 0;
    for (int k=0; k<params.meanInds.size(); k++) {
      if ((std::fabs(resultsR[k]-meanR) <= 2*stdR) 
          && (std::fabs(resultsC[k]-meanC) <= 2*stdC)) {

        centerR_float += resultsR[k];
        centerC_float += resultsC[k];
        Ncents += 1;
      }
    }

    centers[ifl][0] = std::round(centerR_float/Ncents);
    centers[ifl][1] = std::round(centerC_float/Ncents);

    // Show image center
    if (params.pltCent) {
      TH2F* cImg = plt.plotRC(
          imgCent, 
          "plots/centeredImage"
            + to_string(imgINFO[ifl].stagePos), 
          maximum, "1500");
      centerR = std::round(centerR_float);
      centerC = std::round(centerC_float);
      for (int ir=(centerR)-5; ir<=centerR+5; ir++) {
        for (int ic=(centerC)-5; ic<=centerC+5; ic++) {
          if (sqrt(pow(ir-centerR,2)+pow(ic-centerC,2)) < 5) 
            cImg->SetBinContent(ic, ir, 100000);
        }
      }
      cImg->SetMinimum(-1);
      plt.print2d(cImg, 
          "plots/centeredImage"
            + to_string(imgINFO[ifl].stagePos));
      delete cImg;
    }
  }

  // Center mean
  double centerRmeanTemp = 0;
  double centerCmeanTemp = 0;
  for (uint i=0; i<centers.size(); i++) {
    centerRmeanTemp += centers[i][0];
    centerCmeanTemp += centers[i][1];
  }
  centerRmeanTemp /= centers.size();
  centerCmeanTemp /= centers.size();

  // Center std
  double centerRstd = 0;
  double centerCstd = 0;
  for (uint i=0; i<centers.size(); i++) {
    centerRstd += std::pow(centerRmeanTemp-centers[i][0], 2);
    centerCstd += std::pow(centerCmeanTemp-centers[i][1], 2);
  }
  centerRstd = std::sqrt(centerRstd/centers.size());
  centerCstd = std::sqrt(centerCstd/centers.size());

  // Center mean
  double centerRmean = 0;
  double centerCmean = 0;
  Ncents      = 0;
  for (uint i=0; i<centers.size(); i++) {
    if ((fabs(centers[i][0] - centerRmeanTemp) <= 2*centerRstd)
        && (fabs(centers[i][1] - centerCmeanTemp) <= 2*centerCstd)) {

      centerRmean += centers[i][0];
      centerCmean += centers[i][1];
      Ncents++;
    }
  }
  centerRmean /= Ncents;
  centerCmean /= Ncents;


  if (params.verbose) {
    std::cout << "INFO: Searching for global min: "
      << centerRmean << " / " << centerRstd << "    "
      << centerCmean << " / " << centerCstd << std::endl;
  }
  for (ifl=0; ifl<imgINFO.size(); ifl++) {
    for (int i=0; i<params.meanInds.size(); i++) {
      centfnctr.radProc     = &radProc;
      centfnctr.fxnType     = params.centerFxnType;
      centfnctr.meanVal     = allMedVals[ifl][i];
      centfnctr.valSTD      = centerValSTD;
      centfnctr.vals        = &allCentVals[ifl][i];
      centfnctr.indsR       = &allIndsR[ifl][i];
      centfnctr.indsC       = &allIndsC[ifl][i];
      //centfnctr.verbose     = params.verbose;


      double minLoss = 1e90;
      int minR, minC;
      for (int ir=-5; ir<=5; ir++) {
        for (int ic=-5; ic<=5; ic++) {
          std::vector<double> curCent{centerRmean+ir, centerCmean+ic};

          if (minLoss > centfnctr(curCent)) {
            minLoss = centfnctr(curCent);
            minR = curCent[0];
            minC = curCent[1];
          }
        }
      }

      resultsR[i] = minR;
      resultsC[i] = minC;
    }

    meanR = std::accumulate(resultsR.begin(), resultsR.end(), 0);
    meanR /= (float)resultsR.size();
    meanC = std::accumulate(resultsC.begin(), resultsC.end(), 0);
    meanC /= (float)resultsC.size();

    stdR = 0;
    stdC = 0;
    for (uint k=0; k<params.meanInds.size(); k++) {
      stdR += std::pow(meanR-resultsR[k], 2);
      stdC += std::pow(meanC-resultsC[k], 2);
    }
    stdR = std::sqrt(stdR/params.meanInds.size());
    stdC = std::sqrt(stdC/params.meanInds.size());

    centerR_float = 0;
    centerC_float = 0;
    Ncents  = 0;
    for (int k=0; k<params.meanInds.size(); k++) {
      if ((std::fabs(resultsR[k]-meanR) <= 2.5*stdR) 
          && (std::fabs(resultsC[k]-meanC) <= 2.5*stdC)) {

        centerR_float += resultsR[k];
        centerC_float += resultsC[k];
        Ncents += 1;
      }
    }

    centers[ifl][0] = std::round(centerR_float/Ncents);
    centers[ifl][1] = std::round(centerC_float/Ncents);

    /*
    if (fabs(centers[ifl][0] - centersComp[imgINFO[ifl].stagePos][1])
          /centersComp[imgINFO[ifl].stagePos][1] > 0.01 
        || fabs(centers[ifl][1] - centersComp[imgINFO[ifl].stagePos][0])
          /centersComp[imgINFO[ifl].stagePos][0] > 0.01) {

      cout<<"centers do not match\t"<<scan<<"  "<<imgINFO[ifl].stagePos<<"\t"
        <<centers[ifl][0]<<" / "<<centersComp[imgINFO[ifl].stagePos][1]<<" / "
        <<fabs(centers[ifl][0]-centersComp[imgINFO[ifl].stagePos][1])
          /centersComp[imgINFO[ifl].stagePos][1]<<"\t"<<centers[ifl][1]<<" / "
        <<centersComp[imgINFO[ifl].stagePos][0]<<" / "
        <<fabs(centers[ifl][1]-centersComp[imgINFO[ifl].stagePos][0])
          /centersComp[imgINFO[ifl].stagePos][0]<<endl;
    }
    */

  }

  // Center mean
  centerRmeanTemp = 0;
  centerCmeanTemp = 0;
  for (uint i=0; i<centers.size(); i++) {
    centerRmeanTemp += centers[i][0];
    centerCmeanTemp += centers[i][1];
  }
  centerRmeanTemp /= centers.size();
  centerCmeanTemp /= centers.size();

  // Center std
  centerRstd = 0;
  centerCstd = 0;
  for (uint i=0; i<centers.size(); i++) {
    centerRstd += std::pow(centerRmeanTemp-centers[i][0], 2);
    centerCstd += std::pow(centerCmeanTemp-centers[i][1], 2);
  }
  centerRstd = std::sqrt(centerRstd/centers.size());
  centerCstd = std::sqrt(centerCstd/centers.size());

  // Center mean
  centerRmean = 0;
  centerCmean = 0;
  Ncents      = 0;
  for (uint i=0; i<centers.size(); i++) {
    if ((fabs(centers[i][0] - centerRmeanTemp) <= 2.5*centerRstd)
        && (fabs(centers[i][1] - centerCmeanTemp) <= 2.5*centerCstd)) {

      centerRmean += centers[i][0];
      centerCmean += centers[i][1];
      Ncents++;
    }
  }

  centerRmean /= Ncents;
  centerCmean /= Ncents;
  cout<<"centers: "<<centerRmean<<"  "<<centerCmean<<endl;

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

  tree->Branch("run", 	          &curRun);
  tree->Branch("scan", 	          &curScan,         "scan/I");
  tree->Branch("imgNum", 	  &imgNum, 	    "imgNum/I");
  tree->Branch("imgIsRef", 	  &imgIsRef, 	    "imgIsRef/I");
  tree->Branch("timeStamp", 	  &timeStamp, 	    "timeStamp/I");
  tree->Branch("stagePos", 	  &stagePos, 	    "stagePos/I");
  tree->Branch("t0StagePos",      &t0SP,	    "t0StagePos/F");
  tree->Branch("t0Time",	  &t0Time,	    "t0Time/F");
  tree->Branch("throttle",        &throttle,        "throttle/F");
  tree->Branch("centerC", 	  &centerC, 	    "centerC/I");
  tree->Branch("centerR", 	  &centerR, 	    "centerR/I");
  tree->Branch("centerCstdRatio", &centerCstdRatio, "centerCstdRatio/F");
  tree->Branch("centerRstdRatio", &centerRstdRatio, "centerRstdRatio/F");
  tree->Branch("imgNorm", 	  &imgNorm, 	    "imgNorm/F");
  tree->Branch("readoutNoise", 	  &readoutNoise,    "readoutNoise/F");
  tree->Branch("imgOrig", 	  &imgOrig);
  tree->Branch("imgSubBkg",       &imgSubBkg);
  tree->Branch("legCoeffs",       &legCoeffs);
  tree->Branch("rawAzmAvg",       &rawAzmAvg);
  tree->Branch("azmAvg",          &azmAvg);
  //tree->Branch("radPixDist",      &radPixDist);
  //tree->Branch("radPixMeans",     &radPixMeans);
  tree->Branch("imgRadSTD",       &imgRadSTD);

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


  ///////////////////////////////////////
  /////  Processing images in scan  /////
  ///////////////////////////////////////
 
  for (ifl=0; ifl<imgINFO.size(); ifl++) {
    std::cout << "INFO: processing image: " << imgINFO[ifl].fileName << endl;

    //if (imgINFO[ifl].stagePos > 1542900 && imgINFO[ifl].stagePos < 1554000) continue; 

    /////  Check we are in the same run/scan  /////
    if ((imgINFO[ifl].run.compare(curRun) != 0) 
        && (imgINFO[ifl].scan != curScan)) {

      std::cerr << "ERROR new image in run/scan " << imgINFO[ifl].run 
        << "/" << imgINFO[ifl].scan << " instead of " 
        << curRun << "/" << curScan << endl;
      exit(1);
    }

    if (params.verbose) cout << "\tStage position: " 
      + to_string(imgINFO[ifl].stagePos) << endl;


    /////  Check if reference image  /////
    imgIsRef = 0;
    if (params.hasRef) {
      if (imgINFO[ifl].stagePos < params.refStagePosCut) {
        imgIsRef = 1;
      }
    }


    /////  Image number (ordered)  /////
    imgNum = imgINFO[ifl].imgNum;

    /////  Stage position  /////
    stagePos = imgINFO[ifl].stagePos;

    /////  Image capture time  /////
    timeStamp = imgINFO[ifl].time;

    /////  Image radial stds  /////
    for (uint i=0; i<imgRadSTD.size(); i++) {
      imgRadSTD[i] = imgRadSTDs[ifl][i];
    }

    /////  Center  /////
    
    //centerR = std::round(centersComp[stagePos][1]);
    //centerC = std::round(centersComp[stagePos][0]);
    centerR = centers[ifl][0];
    centerC = centers[ifl][1];
    centerRstdRatio = std::fabs(centerR - centerRmean)/centerRstd;
    centerCstdRatio = std::fabs(centerC - centerCmean)/centerCstd;
    centerR = std::round(centerRmean);
    centerC = std::round(centerCmean);

    /*
    if ((centerRstdRatio > params.centerSTDcut) 
        || (centerCstdRatio > params.centerSTDcut)) {

      centerR = centerRmean;
      centerC = centerCmean;
    }
    */



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
    imgAddr = I0fileNames[imgINFO[ifl].stagePos];
    if (params.verbose) cout << "INFO: Trying to open " << imgAddr << "\t .....";
    imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
    if (params.verbose) cout << "\tpassed!\n\n";

    if (params.imgMatType.compare("uint16") == 0) {
      imgI0 = imgProc::getImgVector<uint16_t>(imgMat, 
          imgI0.size(), imgI0.size()/2, imgI0.size()/2, true, NULL,
          -1, params.holeR, params.holeC);
    }
    else if (params.imgMatType.compare("uint32") == 0) {
      imgI0 = imgProc::getImgVector<uint32_t>(imgMat, 
          imgI0.size(), imgI0.size()/2, imgI0.size()/2, true, NULL,
          -1, params.holeR, params.holeC);
    }
    else {
      std::cerr << "ERROR: Do not recognize imgMatType = " 
          + params.imgMatType + "!!!\n";
      exit(1);
    }
 
    imgAddr = imgINFO[ifl].path + imgINFO[ifl].fileName;
    if (params.verbose) cout << "INFO: Trying to open " << imgAddr << "\t .....";
    imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
    if (params.verbose) cout << "\tpassed!\n\n";

    if (params.imgMatType.compare("uint16") == 0) {
      imgOrig = imgProc::getImgVector<uint16_t>(imgMat, 
          Nrows, Nrows/2, Ncols/2, true, NULL,
          params.holeRad, params.holeR, params.holeC);
      imgSubBkg = imgProc::getImgVector<uint16_t>(imgMat, 
          Nrows, Nrows/2, Ncols/2, true, &imgBkg,
          params.holeRad, params.holeR, params.holeC,
          &params.nanMap);
    }
    else if (params.imgMatType.compare("uint32") == 0) {
      imgOrig = imgProc::getImgVector<uint32_t>(imgMat, 
          Nrows, Nrows/2, Ncols/2, true, NULL,
          params.holeRad, params.holeR, params.holeC);
      imgSubBkg = imgProc::getImgVector<uint32_t>(imgMat, 
          Nrows, Nrows/2, Ncols/2, true, &imgBkg,
          params.holeRad, params.holeR, params.holeC,
          &params.nanMap);
    }
    else {
      std::cerr << "ERROR: Do not recognize imgMatType = " 
          + params.imgMatType + "!!!\n";
      exit(1);
    }
   

    // Remove Xray hits
    imgProc::removeXrayHits(
        &imgOrig, 
        params.XrayHighCut, params.XrayLowCut, 
        params.XraySTDcut, params.XrayWindow);
    imgProc::removeXrayHits(
        &imgSubBkg, 
        params.XrayHighCut, params.XrayLowCut, 
        params.XraySTDcut, params.XrayWindow,
        xRayHitHistos);

    if (params.pltVerbose) {
      plt.printRC(imgSubBkg, "./plots/imgSubBkg_original");

      std::vector< vector<double> > testplot(imgSubBkg.size());
      for (int ir=0; ir<imgSubBkg.size(); ir++) {
        testplot[ir].resize(imgSubBkg[ir].size(), 0);
        for (int ic=0; ic<imgSubBkg[ir].size(); ic++) {
          testplot[ir][ic] = imgSubBkg[ir][ic];
          if (fabs(ir-centerR) < 7 && fabs(ic-centerC) < 7) {
            testplot[ir][ic] = -1;
          }
        }
      }

      plt.printRC(testplot, "testingCenter", minimum, "0");
    }




    /////  Remove pixel outliers  /////
    outlierImage = radProc.removeOutliers(imgSubBkg,  
        centerR, centerC, params.imgEdgeBuffer,
        params.NradAzmBins, params.shellWidth, params.Npoly,
        params.stdIncludeLeft, params.distSTDratioLeft,
        params.stdCutLeft, params.meanBinSize,
        params.stdIncludeRight, params.distSTDratioRight,
        params.stdChangeRatio, params.stdCutRight,
        imgINFO[ifl].stagePos, params.outlierMapSTDcut,
        true, params.outlierVerbose, 
        NULL, radPixelHistos);

    outlierImage = radProc.removeOutliersSimple(
        imgSubBkg, centerR, centerC, params.imgEdgeBuffer,
        params.NradAzmBins, params.shellWidth, false, params.Npoly,
        params.outlierSTDcut, imgINFO[ifl].stagePos, params.outlierMapSTDcut,
        true, (false || params.outlierVerbose), NULL);//, &plt);
        
    /*
    radProc.removeOutliers(imgOrig,  
        centerR, centerC, params.imgEdgeBuffer, 
        params.NradAzmBins, params.shellWidth, params.Npoly,
        params.stdIncludeLeft, params.distSTDratioLeft,
        params.stdCutLeft, params.meanBinSize,
        params.stdIncludeRight, params.distSTDratioRight,
        params.stdChangeRatio, params.stdCutRight,
        imgINFO[ifl].stagePos, params.outlierMapSTDcut,
        false, params.outlierVerbose, 
        NULL, NULL);

 
        */

    
    if (params.pltVerbose) {
      //plt.printRC(outlierImage, "outlierSTD_" + imgINFO[ifl].run + "-" + to_string(imgINFO[ifl].scan) + "-" + to_string(imgINFO[ifl].stagePos));
      save::saveDat<double>(outlierImage, 
            "./results/outlierSTD-" + runName
                + "-" + to_string(imgINFO[ifl].scan)
                + "-" + to_string(imgINFO[ifl].stagePos) + ".dat");
    }

    ///  Removing small clusters  ///
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
    if (params.verbose) std::cout << "INFO: Readout noise subtraction.\n";
    readoutNoise = imgProc::removeAvgReadOutNoise(imgSubBkg, centerR, centerC, 
    //readoutNoise = imgProc::removeMedianReadOutNoise(imgSubBkg, centerR, centerC, 
                      params.readoutAzmBinStart, params.readoutAzmBinEnd, 
                      params.imgEdgeBuffer, &params.nanMap);

    // Only use this when interested in reference images
    //    background fit fails when subtracting readout
    //    noise on the phosphor screen. Comment previous line.
    //readoutNoise = imgProc::removeReadOutNoise(imgSubBkg);


    //////////////////////////////////////////////////////
    /////  Finding and subtracting laser background  /////
    //////////////////////////////////////////////////////
    // We cluster background spots based on the assymetry of 
    //    the image. 

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
    //for (int k=0; k<params.NradAzmBins; k++) {
    //  radPixDist[radInd].clear();
    //}
    for (int ir=0; ir<(int)imgSubBkg.size(); ir++) {
      if (ir < params.imgEdgeBuffer 
          || imgSubBkg.size() - ir < params.imgEdgeBuffer) continue;
      for (int ic=0; ic<(int)imgSubBkg[ir].size(); ic++) {
        if (ic < params.imgEdgeBuffer 
            || imgSubBkg[ir].size() - ic < params.imgEdgeBuffer) continue;

        if (imgSubBkg[ir][ic] != NANVAL) {
          if (params.nanMap[ir][ic] == NANVAL) {
            continue; 
          }

          radInd = std::round(
                      std::sqrt(
                        std::pow(ir-centerR,2) + std::pow(ic-centerC,2)));

          if (radInd < params.NradAzmBins) {
            rawAzmAvg[radInd] += imgSubBkg[ir][ic];
            azmCounts[radInd] += 1;
            //radPixDist[radInd].push_back(imgSubBkg[ir][ic]); 
          }
        }
      }
    }
    /*
    for (int k=0; k<params.NradAzmBins; k++) {
      radPixMeans[k] = std::accumulate(radPixDist[k].begin(), 
                                radPixDist[k].end(), 0)
                                /((float)radPixDist[k].size());
    }
    */
    for (uint ir=0; ir<rawAzmAvg.size(); ir++) {
      if (azmCounts[ir] != 0) {
        rawAzmAvg[ir] /= azmCounts[ir];
      }
      else rawAzmAvg[ir] = NANVAL;
    }

    /////  Image norm  /////
    imgNorm = 0;
    count = 0;
    for (int i=imgNormBinMin; i<imgNormBinMax; i++) {
      if (rawAzmAvg[i] == NANVAL) continue;
      imgNorm += rawAzmAvg[i]/atmDiff[i];
      count++;
    }
    imgNorm /= count;

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      if (rawAzmAvg[iq] != NANVAL) {
        rawAzmAvg[iq] /= imgNorm;
      }
      azmAvg[iq] = rawAzmAvg[iq];
    }

    /////  Plotting azm lines  /////
    cout<<"NORM: "<<stagePos<<"  "<<imgNorm<<endl;
    pLO = radProc.getPolarLineOut(&imgSubBkg,
    //pLO = radProc.getPolarLineOut(&imgOrig,
            centerR, centerC, 168, 45, 200);
    for (int k=0; k<pLO.size(); k++) {
      pLO[k] /= imgNorm;
    }
    save::saveDat<double>(pLO, dirName 
        + "scan-" + to_string(curScan)
        + "_stagePos-" + to_string(stagePos)
        + "_pixRange-168-213_Bins[" + to_string(pLO.size()) + "].dat"); 

    pLO = radProc.getPolarLineOut(&imgSubBkg,
    //pLO = radProc.getPolarLineOut(&imgOrig,
            centerR, centerC, 222, 44, 200);
    for (int k=0; k<pLO.size(); k++) {
      pLO[k] /= imgNorm;
    }
    save::saveDat<double>(pLO, dirName 
        + "scan-" + to_string(curScan)
        + "_stagePos-" + to_string(stagePos)
        + "_pixRange-222-266_Bins[" + to_string(pLO.size()) + "].dat"); 

    //pLO = radProc.getPolarLineOut(&imgOrig,
    pLO = radProc.getPolarLineOut(&imgSubBkg,
            centerR, centerC, 266, 45, 200);
    for (int k=0; k<pLO.size(); k++) {
      pLO[k] /= imgNorm;
    }
    save::saveDat<double>(pLO, dirName 
        + "scan-" + to_string(curScan)
        + "_stagePos-" + to_string(stagePos)
        + "_pixRange-266-311_Bins[" + to_string(pLO.size()) + "].dat"); 

    /////  Build References  /////
    if (imgIsRef) {
      for (int iq=0; iq<params.NradAzmBins; iq++) {
        refAzmAvg[iq] = (refAzmAvg[iq]*Nref + rawAzmAvg[iq])/(Nref + 1);
      }
      Nref++;
      I0refImgs[stagePos].resize(imgI0.size());
      for (uint ir=0; ir<imgI0.size(); ir++) {
        I0refImgs[stagePos][ir].resize(imgI0[ir].size(), NANVAL);
        for (uint ic=0; ic<imgI0.size(); ic++) {
          if (imgI0[ir][ic] != NANVAL) {
            imgI0ref[ir][ic] += imgI0[ir][ic];
            imgI0refCount[ir][ic]++;
          }
          I0refImgs[stagePos][ir][ic] = imgI0[ir][ic];
        }
      }
    }
    else if (!normI0ref) {
      for (uint ir=0; ir<imgI0.size(); ir++) {
        for (uint ic=0; ic<imgI0.size(); ic++) {
          if (imgI0refCount[ir][ic]) {
            imgI0ref[ir][ic] /= imgI0refCount[ir][ic];
          }
          else {
            imgI0ref[ir][ic] = NANVAL;
          }
        }
      }

      for (auto & itr : I0refImgs) {
        for (uint ir=0; ir<imgI0.size(); ir++) {
          for (uint ic=0; ic<imgI0.size(); ic++) {
            if ((itr.second[ir][ic] != NANVAL) && (imgI0ref[ir][ic] != NANVAL)) {
              itr.second[ir][ic] -= imgI0ref[ir][ic];
            }
          }
        }

        save::saveDat<double>(itr.second,
            params.preProcI0OutputDir
            + "data-" + runName
            + "_I0refSubtracted_scan-" + to_string(curScan)
            + "_stagePos-" + to_string(itr.first)
            + "_bins[" + to_string(itr.second.size())
            + "," + to_string(itr.second[0].size()) + "].dat");
      }

      normI0ref = true;
    }


    /////  I0 Results  /////
    save::saveDat<double>(imgI0,
        params.preProcI0OutputDir
        + "data-" + runName
        + "_I0_scan-" + to_string(curScan)
        + "_stagePos-" + to_string(stagePos)
        + "_bins[" + to_string(imgI0.size())
        + "," + to_string(imgI0[0].size()) + "].dat");

    if (!imgIsRef) {
      for (uint ir=0; ir<imgI0.size(); ir++) {
        for (uint ic=0; ic<imgI0.size(); ic++) {
          if ((imgI0[ir][ic] != NANVAL) && (imgI0ref[ir][ic] != NANVAL)) {
            imgI0[ir][ic] -= imgI0ref[ir][ic];
          }
        }
      }

      save::saveDat<double>(imgI0,
          params.preProcI0OutputDir
          + "data-" + runName
          + "_I0refSubtracted_scan-" + to_string(curScan)
          + "_stagePos-" + to_string(stagePos)
          + "_bins[" + to_string(imgI0.size())
          + "," + to_string(imgI0[0].size()) + "].dat");
    }

    /////  Low Pass Filtering  /////
    std::vector<double> tstR(filtFFToutSize);
    if (!padRange) {
      while (rawAzmAvg[padRange] == NANVAL) {
        padRange++;
      }
    }
    curPadRange = padRange;
    while (((rawAzmAvg[curPadRange+1] - rawAzmAvg[curPadRange])
            /rawAzmAvg[curPadRange]) > 0.2) {
    //while (fabs(1 - rawAzmAvg[curPadRange]/rawAzmAvg[curPadRange+1]) > 0.2) {
      curPadRange++;
    }


    // Change DC so pattern tends to 0
    double filtBaseline = 0;
    for (int iq=(int)(0.98*params.NradAzmBins); iq<params.NradAzmBins; iq++) {
      filtBaseline += rawAzmAvg[iq];
    }
    filtBaseline /= params.NradAzmBins - (int)(0.98*params.NradAzmBins);

    qSpace[0] = 0;
    for (int iq=1; iq<params.NradAzmBins; iq++) {
      qSpace[iq] = (rawAzmAvg[iq] - filtBaseline)/atmDiff[iq];
      if (iq < params.suppressBins) {
        if (iq < curPadRange && iq > 0) {
          qSpace[iq] = (rawAzmAvg[curPadRange] - filtBaseline)
              *(pow(curPadRange,4)/atmDiff[iq])/pow(iq,4);
            //params.padMaxHeight*sin(PI/2*iq/padRange);
        }
        qSpace[iq] *= pow(sin((PI/2)*iq/params.suppressBins), 6);
      }
    }


    if (params.pltFilterVerbose) {
      std::vector<double> tst(params.NradAzmBins);
      for (int iq=0; iq<params.NradAzmBins; iq++) {
        tst[iq] = qSpace[iq];
      }
      pltVals[0] = "-2";
      pltVals[1] = "2";
      delete plt.print1d(tst, 
          "./plots/filtInput_" + to_string(stagePos),
          pltOpts, pltVals);
    }

    fftw_execute(filtFFTf);

    for (int ir=0; ir<filtFFToutSize; ir++) {
      rSpace[ir][0] *= bandPassFilter[ir]/sqrt(params.NradAzmBins);
      rSpace[ir][1] *= bandPassFilter[ir]/sqrt(params.NradAzmBins);
    }

    fftw_execute(filtFFTb);

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      if (azmAvg[iq] != NANVAL) {
        azmAvg[iq] = qSpace[iq]*atmDiff[iq]/sqrt(params.NradAzmBins);

        if (iq < params.suppressBins) {
          azmAvg[iq] /= pow(sin((PI/2)*iq/params.suppressBins), 6);
        }
        azmAvg[iq] += filtBaseline;

        if (iq < padRange) {
          azmAvg[iq] = NANVAL;
        }
      }
    }

    ///  Send high Q to 0 for filtered results  ///
    filtBaseline = 0;
    for (int iq=params.readoutAzmBinStart; iq<params.readoutAzmBinEnd; iq++) {
      filtBaseline += azmAvg[iq];
    }
    filtBaseline /= params.readoutAzmBinEnd - params.readoutAzmBinStart;

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      if (azmAvg[iq] != NANVAL) {
        azmAvg[iq] -= filtBaseline;
      }
    }

    ///  Image norm filtered result ///
    double filtImgNorm = 0;
    count = 0;
    for (int i=imgNormBinMin; i<imgNormBinMax; i++) {
      if (azmAvg[i] == NANVAL) continue;
      filtImgNorm += azmAvg[i]/atmDiff[i];
      count++;
    }
    filtImgNorm /= count;

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      if (azmAvg[iq] != NANVAL) {
        azmAvg[iq] /= filtImgNorm;
      }
    }

    ///  Updating and applying image norm  ///
    imgNorm *= filtImgNorm;
    /*
    for (uint i=0; i<radPixMeans.size(); i++) {
      radPixMeans[i] /= imgNorm;
    }
    for (uint i=0; i<radPixDist.size(); i++) {
      for (uint j=0; j<radPixDist[i].size(); j++) {
        radPixDist[i][j] /= imgNorm;
      }
    }
    */

    ///  Plotting Filter Results  ///
    if (params.pltFilterVerbose) {
      std::vector<double> tst1(params.NradAzmBins, 0);
      std::vector<double> tst2(params.NradAzmBins, 0);
      vector<TH1*> plts(2);
      for (int iq=1; iq<params.NradAzmBins; iq++) {
        tst1[iq] = rawAzmAvg[iq]/atmDiff[iq];
        tst2[iq] = azmAvg[iq]/atmDiff[iq];
      //cout<<"tst2: "<<iq<<"  "<<tst2[iq]<<endl;
      }
      pltVals[0] = "-0.2";
      pltVals[1] = "0.2";
      plts[0] = plt.plot1d(tst1, "blah_" + to_string(stagePos),pltOpts, pltVals);
      plts[1] = plt.print1d(tst2, "./plots/filtOutp_" + to_string(stagePos),pltOpts, pltVals);
      plts[1]->SetLineColor(4);
      plts[0]->SetMaximum(40);
      plt.print1d(plts, "./plots/filtCompare_" + to_string(stagePos));
      plts[0]->SetMinimum(1);
      plt.print1d(
          plts, "./plots/filtCompareLogy_" + to_string(stagePos), 
          logy, "true");
      delete plts[0];
      delete plts[1];
    }

    /////  Plotting results  /////
    if (params.pltVerbose) {
      //plt.print1d(rawAzmAvg, "plots/azimuthalAvg_" + to_string(imgINFO[ifl].imgNum));
      plt.print1d(rawAzmAvg, "tazimuthalAvg_" + to_string(stagePos), logy, "true");
    }
 
    if (imgIsRef || stagePos == 1565700) {
      std::vector< std::vector<double> > temp(imgSubBkg.size());
      for (uint kr=0; kr<imgSubBkg.size(); kr++) {
        temp[kr].resize(imgSubBkg[kr].size());
        for (uint kc=0; kc<imgSubBkg[kr].size(); kc++) {
          if (imgSubBkg[kr][kc] != NANVAL) {
            temp[kr][kc] = imgSubBkg[kr][kc]/(imgNorm*filtImgNorm);
          }
          else {
            temp[kr][kc] = NANVAL;
          }
        }
      }
      save::saveDat<double>(temp,
          "/reg/ued/ana/scratch/nitroBenzene/refTesting/run-"
          + runName + "_scan-" + to_string(curScan)
          + "_stagePos-" + to_string(stagePos) + "_bins[1024,1024].dat");
    }


    // Save image and info to tree
    tree->Fill();
  }


  /////  Cleaning Up  /////
  tree->Write();
  file->Close();

  if (params.xRayHitDist) {
    std::vector<PLOToptions> xOpts(2);
    std::vector<std::string> xVals(2);
    xOpts[0] = logy;    xVals[0] = "true";
    xOpts[1] = xLabel;  xVals[1] = "Pixel Value";
    delete plt.print1d(xRayHitHistos[0], 
        "./plots/xRayHitDist_highCut_scan-"+to_string(curScan), 
        xOpts, xVals);
    delete plt.print1d(xRayHitHistos[1], 
        "./plots/xRayHitDist_lowCut_scan-"+to_string(curScan), 
        xOpts, xVals);
    delete[] xRayHitHistos;
  }

  if (params.plotRadPixDist) {
    std::vector<PLOToptions> xOpts(3);
    std::vector<std::string> xVals(3);
    xOpts[0] = logy;    xVals[0] = "true";
    xOpts[1] = xLabel;  xVals[1] = "Pixel Value";
    xOpts[2] = xBinSpan;   
    for (int i=0; i<params.NradAzmBins; i++) {
      xVals[2] = to_string(radPixelHistos[i]->FindFirstBinAbove(0.5) - 3)
        + "," + to_string(radPixelHistos[i]->FindLastBinAbove(0.5) + 3);
      radPixelHistos[i]->SetMinimum(0.5);
      delete plt.print1d(radPixelHistos[i], 
        "./plots/radPixDist_rad-" + to_string(i) 
          + "_scan-"+to_string(curScan), 
        xOpts, xVals);
    }
    delete[] radPixelHistos;
  }


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

                                        
