#include "preProcessing.h"
#include "/reg/neh/home/khegazy/baseTools/tools/parameters.h"

using namespace std;


int main(int argc, char* argv[]) {

  // Check for arguments
  if (argc < 2) {
    std::cerr << "ERROR: Must run program a list of files" 
      << " ./preProcessing.exe runList.txt!!!\n";
    exit(1);
  }



  // Select how to run code
  bool localTesting   = false;
  bool doRunLists     = false;
  bool doBackground   = false;
  bool computeCenters = false;


  // Get command line inputs
  std::string runListName(argv[1]);
  if (argc > 2) {
    for (int iarg=2; iarg<argc; iarg+=2) {
      if (strcmp(argv[iarg], "-computeCenters") == 0)
        {computeCenters = atoi(argv[iarg+1]);}
      else {
        std::cerr << "ERROR!!! Option " << argv[iarg] << " does not exist!\n";
        exit(0);
      }
    }
  }


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
  std::vector<double> pLO;

  //plt.printRC(params.nanMap, "testingNANmap");

  const int Nlegendres = 1;
  const int NradLegBins = 50;

  // Indices
  //const int imgCut = 935;
  const int imgCut = 895;

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
  if (imgCut < params.imgSize) {
    std::cerr 
      << "ERROR: parameter imgSize must be >= imgCut!!!\n";
    exit(1);
  }


   // Make sure the image has odd number of bins
  if (!(imgCut%2) || !(params.imgSize%2)) {
    std::cerr 
      << "ERROR: imgCut and params.imgSize must be an odd number!!!\n";
    exit(1);
  }

  imgProc::radProcTool radProc(params.shellWidth, params.NradAzmBins);
  
  uint ifl;
  std::vector<PLOToptions> pltOpts(2);
  std::vector<string> pltVals(2);


  std::vector<double> legCoeffs(params.NradLegBins);
  std::vector<double> rawAzmAvg(params.NradAzmBins);
  std::vector<int> rawAzmAvg_nanMap(params.NradAzmBins);
  std::vector<double> azmAvg(params.NradAzmBins);
  std::vector<int> azmAvg_nanMap(params.NradAzmBins);
  std::vector<double> azmCounts(params.NradAzmBins);

  std::vector< std::vector<double> > oddImgImgn;
  std::vector< std::vector<double> > symImg(imgCut);
  std::vector< std::vector<double> > stdRatioImg(imgCut);
  std::vector< std::vector<double> > outlierImage;
  std::vector< std::vector<double> > outlierBkg;
  for (uint ir=0; ir<imgCut; ir++) { 
    symImg[ir].resize(imgCut, 0);
    stdRatioImg[ir].resize(imgCut, 0);
  }

  fftw_complex* fftIn = 
      (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*imgCut*imgCut);
  fftw_complex* fftOut = 
      (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*imgCut*imgCut);
  fftw_plan fftFref = fftw_plan_dft_2d(
      imgCut, imgCut, 
      fftIn, fftOut, 
      FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan fftBref = fftw_plan_dft_2d(
      imgCut, imgCut, 
      fftIn, fftOut, 
      FFTW_BACKWARD, FFTW_MEASURE);


  std::vector< std::vector<double> > radBins;

  imgProc::centerfnctr centfnctr;
  std::vector<double> center(2);
  std::vector<long int> inpCenter(2);
  int centerC, centerR;
  int I0centerC, I0centerR;
  float centerC_float, centerR_float;
  float centerRstdRatio, centerCstdRatio;
  int Ncents;

  double count;

  int imgIsRef;

  std::map< std::string, std::vector<double> > pvVals, pvValsDer;
  std::map< std::string, long int > pvStartTimes;
  std::map< std::string, std::vector<float> > pvAll, pvAllDer;
  float throttle;

  int timeStamp;

  int radInd;
  int curScan;
  std::string fileName;
  std::string line;
  size_t ipos, fpos;
  std::string date, scan, curRun;
  std::string results_folder, results_file_prefix, radial_base_folder; 
  float t0SP, t0Time;
  std::vector<imgProc::imgInfoStruct> imgINFO;
  std::map<int, std::string> I0fileNames;
  float readoutNoise;

  int imgNormBinMin = (int)(params.imgNormRadMin*params.NradAzmBins);
  int imgNormBinMax = (int)(params.imgNormRadMax*params.NradAzmBins);

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
      + "_order-" + to_string(params.order)
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
    std::string num = to_string(params.maxQazm);
    num = num.substr(0, 5);
    save::importDat<double>(atmDiff, params.simOutputDir 
      + "/" + params.molName 
      + "_sim_atmDiffraction-lineout_align-random_Qmax-" + num
      + "_Bins[" + to_string(params.NradAzmBins) + "].dat");

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      sMsNorm[iq] = (params.maxQazm*(iq + 0.5)/params.NradAzmBins)/atmDiff[iq];
    }
  }



  ////////////////////////////////
  ////  Retrieving File Info  ////
  ////////////////////////////////

  if (params.verbose) cout<<"Retrieving parameters"<<endl;
  ppFunct::getScanRunInfo(imgINFO, runListName, params.verbose);
  
  if (params.hasI0) {
    ppFunct::getI0RunInfo(I0fileNames, runListName, params.verbose);
  }

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

  if (params.verbose) cout<<"Setting image variables"<<endl;
  std::string imgAddr;
  cv::Mat imgMat;
  int Nrows, Ncols;
  bool normI0ref = false;
  int I0refCount;
  std::vector< std::vector<double> > imgI0;
  std::vector< std::vector<double> > imgI0ref;
  //std::vector< std::vector<double> > imgI0refCount;
  std::map<int, std::vector< std::vector<double> > > I0refImgs;

  if (params.hasI0) {
    imgAddr = I0fileNames[imgINFO[0].imgNum];
    imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
    Nrows = imgMat.rows;
    Ncols = imgMat.cols;
    imgI0.resize(Nrows);
    imgI0ref.resize(Nrows);
    //imgI0refCount.resize(Nrows);
    for (int ir=0; ir<Nrows; ir++) {
      imgI0[ir].resize(Ncols, 0);
      imgI0ref[ir].resize(Ncols, 0);
      //imgI0refCount[ir].resize(Ncols, 0);
    }
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
    std::vector< std::vector<int> > bkg_nanMap;
    for (ifl=0; ifl<imgINFO.size(); ifl++) {
      ///  Get image  ///
      imgAddr   = imgINFO[ifl].path + imgINFO[ifl].fileName;
      imgMat    = cv::imread(imgAddr.c_str(), CV_LOAD_IMAGE_ANYDEPTH); 
      Nrows = imgMat.rows;
      Ncols = imgMat.cols;
      std::vector< std::vector<int> > img_nanMap(Nrows);
      for (uint ir=0; ir<Nrows; ir++) {
        img_nanMap[ir].resize(Ncols, 0);
      }

      if (ifl == 0) {
        bkgCount.resize(Nrows);
        bkg_nanMap.resize(Nrows);;
        for (uint ir=0; ir<Nrows; ir++) {
          bkgCount[ir].resize(Ncols, 0);
          bkg_nanMap[ir].resize(Ncols, 0);
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
      img_nanMap = imgProc::removeXrayHits(
          &imgVec, img_nanMap,
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
          if (img_nanMap[ir][ic]) {
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
          if (img_nanMap[ir][ic]) {
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
          if (img_nanMap[ir][ic]) {
            continue;
          }

          if (fabs(imgVec[ir][ic] - rsMean) > params.bkgSTDcut*rsSTD) {
            img_nanMap[ir][ic] = 1;
          }
        }
      }

      ///  Add image  ///
      if (params.verbose) std::cout << "\tAdding images\n";
      for (int ir=0; ir<Nrows; ir++) {
        for (int ic=0; ic<Ncols; ic++) {
          if (img_nanMap[ir][ic] == 0) {
            imgBkg[ir][ic] += imgVec[ir][ic];
            bkgCount[ir][ic] += 1;
          }
        }
      }

      for (int ir=0; ir<Nrows; ir++) {
        for (int ic=0; ic<Ncols; ic++) {
          if (img_nanMap[ir][ic]) {
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
          bkg_nanMap[ir][ic] = 1;
        }
      }
    }

    ///  Fill in pixels with NANVALS by averaging nearest neighbors  ///
    for (int ir=0; ir<Nrows; ir++) {
      for (int ic=0; ic<Ncols; ic++) {
        if (bkg_nanMap[ir][ic]) {
          std::vector<double> collection;
          int ind = 1;
          while (collection.size() < 15) {
            if (ic - ind >= 0) {
              for (int irr=-1*ind; irr<=ind; irr++) {
                if ((ir + irr >= 0) && (ir + irr < (int)imgBkg.size())) {
                  if (bkg_nanMap[ir+irr][ic-ind] == 0) {
                    collection.push_back(imgBkg[ir+irr][ic-ind]);
                  }
                }
              }
            }
            if (ic + ind < (int)imgBkg.size()) {
              for (int irr=-1*ind; irr<=ind; irr++) {
                if ((ir + irr >= 0) && (ir + irr < (int)imgBkg.size())) {
                  if (bkg_nanMap[ir+irr][ic+ind] == 0) {
                    collection.push_back(imgBkg[ir+irr][ic+ind]);
                  }
                }
              }
            }
            if (ir - ind >= 0) {
              for (int icc=-1*ind+1; icc<=ind-1; icc++) {
                if ((ic + icc >= 0) && (ic + icc < (int)imgBkg.size())) {
                  if (bkg_nanMap[ir-ind][ic+icc] == 0) {
                    collection.push_back(imgBkg[ir-ind][ic+icc]);
                  }
                }
              }
            }
            if (ir + ind < (int)imgBkg.size()) {
              for (int icc=-1*ind+1; icc<=ind-1; icc++) {
                if ((ic + icc >= 0) && (ic + icc < (int)imgBkg.size())) {
                  if (bkg_nanMap[ir+ind][ic+icc] == 0) {
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

    delete plt.printRC(imgBkg, "plots/background-" + imgINFO[0].run);

    exit(0);
  }


  //////////////////////////////////
  /////  Get background image  /////
  //////////////////////////////////

  //cerr << "WARNING!!!!! NOT IMPORTING BKG"<<endl;
  if (params.backgroundImage.compare("NULL") != 0) {
    if (params.verbose) cout << "INFO: Importing Background\n";
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


  /////  I0 image norm variables  /////
  float I0norm;
  std::vector<double> I0norms(imgINFO.size());


  /////////////////////////////
  /////  Centering Image  /////
  /////////////////////////////

  std::vector< std::vector<int> > centers(imgINFO.size());
  std::vector< std::vector<int> > I0centers(imgINFO.size());
  for (ifl=0; ifl<imgINFO.size(); ifl++) {
    centers[ifl].resize(2);
    I0centers[ifl].resize(2);
  }

  double centerRmeanTemp = 0;
  double centerCmeanTemp = 0;
  double centerRmean = 0;
  double centerCmean = 0;
  double centerRstd = 0;
  double centerCstd = 0;



  if (computeCenters || params.computeCenters || params.scanAvgCenter) {

    /////  Find Center  /////
    //imgProc::centerfnctr centfnctr;
    double meanR, meanC;
    double stdR, stdC;
    double centerValSTD = 2;
    std::vector<double> resultsR(params.meanInds.size());
    std::vector<double> resultsC(params.meanInds.size());
    std::vector< std::vector<double> > allMedVals(imgINFO.size());
    std::vector< std::vector< std::vector<int> > > allIndsR(imgINFO.size());
    std::vector< std::vector< std::vector<int> > > allIndsC(imgINFO.size());
    std::vector< std::vector< std::vector<double> > > allCentVals(imgINFO.size());
    int rad;
    std::vector< std::vector<int> > centImg_nanMap(Nrows);
    for (int ir=0; ir<Nrows; ir++) {
      centImg_nanMap[ir].resize(Ncols, 0);
    }
    for (ifl=0; ifl<imgINFO.size(); ifl++) {
      imgAddr = imgINFO[ifl].path + imgINFO[ifl].fileName;
      imgMat = cv::imread(imgAddr.c_str(), CV_LOAD_IMAGE_ANYDEPTH); 

      Nrows = imgMat.rows;
      Ncols = imgMat.cols;
      if (params.imgMatType.compare("uint16") == 0) {
        imgCent = imgProc::getImgVector<uint16_t>(imgMat, //imgSmooth, 
                    Nrows, Nrows/2, Ncols/2, true, &imgBkg);
                    //params.holeRad, params.holeR, params.holeC,
                    //&params.nanMap);
      }  
      else if (params.imgMatType.compare("uint32") == 0) {
        imgCent = imgProc::getImgVector<uint32_t>(imgMat, //imgSmooth, 
                    Nrows, Nrows/2, Ncols/2, true, &imgBkg);
                    //params.holeRad, params.holeR, params.holeC,
                    //&params.nanMap);
      }  
      else {
        std::cerr << "ERROR: Do not recognize imgMatType = " 
            + params.imgMatType + "!!!\n";
        exit(1);
      }

      // Create initial nanMap for centering images
      for (int ir=0; ir<Nrows; ir++) {
        for (int ic=0; ic<Ncols; ic++) {
          rad = (int)(sqrt(std::pow(ir-params.holeR, 2)
              + std::pow(ic-params.holeC, 2)));

          if (params.nanMap[ir][ic] || (rad < params.holeRad)) {
            centImg_nanMap[ir][ic] = 1;
          }
          else {
            centImg_nanMap[ir][ic] = 1;
          }
        }
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
                    Nrows, Nrows/2, Ncols/2, true, NULL);
                    //params.holeRad+5, params.holeR, params.holeC,
                    //&params.nanMap);
      }  
      else if (params.imgMatType.compare("uint32") == 0) {
        imgCent = imgProc::getImgVector<uint32_t>(imgSmooth, 
                    Nrows, Nrows/2, Ncols/2, true, NULL);
                    //params.holeRad+5, params.holeR, params.holeC,
                    //&params.nanMap);
      }  
      else {
        std::cerr << "ERROR: Do not recognize imgMatType = " 
            + params.imgMatType + "!!!\n";
        exit(1);
      }

      // Create initial nanMap for centering images
      for (int ir=0; ir<Nrows; ir++) {
        for (int ic=0; ic<Ncols; ic++) {
          rad = (int)(sqrt(std::pow(ir-params.holeR, 2)
              + std::pow(ic-params.holeC, 2)));

          if (params.nanMap[ir][ic] || (rad < params.holeRad+5)) {
            centImg_nanMap[ir][ic] = 1;
          }
          else {
            centImg_nanMap[ir][ic] = 0;
          }
        }
      }

      if (params.pltVerbose) {
        plt.printRC(imgCent, "testingCenter");
      }
      std::vector<int> centersCOM(2);
      if (params.doCOM_center) {
        centersCOM = imgProc::centerSearchCOM(imgCent,
            params.centR_estimate, params.centC_estimate,
            params.minRad, params.maxRad,
            params.meanInd, params.COMstdScale, true, &plt);
      }
      else {
        centersCOM[0] = params.centR_estimate;
        centersCOM[1] = params.centC_estimate;
      }


      cout<<"COM: "<<centersCOM[0]<<"  "<<centersCOM[1]<<endl;
      centerR_float = 0;
      centerC_float = 0;
      Ncents  = 0;
      double meanCentVal;
      int meanCentCount;
      std::vector<double> Rs, centVals;
      std::vector<double> circVals(4);
      std::vector<double> orderedInds(4);
      std::vector<int> indsR, indsC, centerRemoveInds;
      allIndsR[ifl].resize(params.meanInds.size());
      allIndsC[ifl].resize(params.meanInds.size());
      allCentVals[ifl].resize(params.meanInds.size());
      allMedVals[ifl].resize(params.meanInds.size());
      for (int i=0; i<params.meanInds.size(); i++) {
        meanCentVal = 0;
        meanCentCount = 0;
        if (centImg_nanMap[centersCOM[0]+params.meanInds[i]][centersCOM[1]] == 0) {
          meanCentVal += imgCent[centersCOM[0]+params.meanInds[i]][centersCOM[1]];
          meanCentCount += 1;
        }
        if (centImg_nanMap[centersCOM[0]-params.meanInds[i]][centersCOM[1]] == 0) {
         meanCentVal += imgCent[centersCOM[0]-params.meanInds[i]][centersCOM[1]];
         meanCentCount += 1;
        }
        if (centImg_nanMap[centersCOM[0]][centersCOM[1]+params.meanInds[i]] == 0) {
          meanCentVal += imgCent[centersCOM[0]][centersCOM[1]+params.meanInds[i]];
          meanCentCount += 1;
        }
        if (centImg_nanMap[centersCOM[0]][centersCOM[1]-params.meanInds[i]] == 0) {
          meanCentVal += imgCent[centersCOM[0]][centersCOM[1]-params.meanInds[i]];
          meanCentCount += 1;
        }
        //iota(orderedInds.begin(), orderedInds.end(), 0);
        //sort(orderedInds.begin(), orderedInds.end(),
        //    [&circVals](int i1, int i2)
        //    {return circVals[i1] < circVals[i2];});
        //meanCentVal = (circVals[orderedInds[1]] + circVals[orderedInds[2]])/2.;
        meanCentVal /= meanCentCount;

        // Find indices with values close to meanCentVal
        indsR.clear();
        indsC.clear();
        centVals.clear();
        for (int ir=0; ir<Nrows; ir++) {
          for (int ic=0; ic<Ncols; ic++) {
            if (centImg_nanMap[ir][ic] == 0) {
              if (fabs(imgCent[ir][ic] - meanCentVal) <= 4*centerValSTD) {
                indsR.push_back(ir);
                indsC.push_back(ic);
                centVals.push_back(imgCent[ir][ic]);
              }
            }
          }
        }

        // Remove indices far from cluster
        center[0] = centersCOM[0];
        center[1] = centersCOM[1];
        for (int j=0; j<3; j++) {
          
          // Find best center for given inds
          centfnctr.radProc     = &radProc;
          centfnctr.fxnType     = params.centerFxnType;
          centfnctr.meanVal     = meanCentVal;
          centfnctr.valSTD      = centerValSTD;
          centfnctr.vals        = &centVals;
          centfnctr.indsR       = &indsR;
          centfnctr.indsC       = &indsC;
          //centfnctr.verbose     = params.verbose;

          if (params.verbose) std::cout << "\t\t Searching for center (PowellMin) ... ";
          tools::powellMin<imgProc::centerfnctr> (centfnctr, center,
              params.cntrScale, params.cntrMinScale,
              params.cntrPowellTol, params.cntrFracTol1d);
          if (params.verbose) std::cout << center[0] << "  " << center[1] << std::endl;


          // Prune inds that are far from the center
          meanR = 0;
          Rs.resize(indsR.size());
          for (uint k=0; k<indsR.size(); k++) {
            Rs[k] = std::sqrt(
                        std::pow(indsR[k] - center[0], 2)
                        + std::pow(indsC[k] - center[1], 2));
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
        allMedVals[ifl][i] = meanCentVal;

        centfnctr.radProc     = &radProc;
        centfnctr.fxnType     = params.centerFxnType;
        centfnctr.meanVal     = meanCentVal;
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

    /////  Calculate the mean and std of the centers for global min  /////
    
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


    }

    std::string I0centers_filename = params.centerDir
        + "I0-centers_run-" + curRun + "_scan-" + to_string(curScan)
        + "_shape[" + to_string(I0centers.size())
        + "," + to_string(I0centers[0].size()) + "].dat";

    if (params.hasI0) {
      if (params.I0centers || !tools::fileExists(I0centers_filename)) {

        // Save centers
        save::saveDat<int>(centers, params.centerDir 
            + "rawImg-centers_run-" + curRun + "_scan-" + to_string(curScan)
            + "_shape[" + to_string(centers.size()) 
            + "," + to_string(centers[0].size()) + "].dat");
       
        /////   I0 Centers   /////
        std::vector<int> I0stagePos;
        for (ifl=0; ifl<imgINFO.size(); ifl++) {
          imgAddr = I0fileNames[imgINFO[ifl].imgNum];
          system(("python I0centriod.py --fileName=" + imgAddr
              + " --minPixVal=" + to_string(params.I0minPixVal)
              + " --approxR=" + to_string(params.I0approxR)
              + " --approxC=" + to_string(params.I0approxC)).c_str());
          I0stagePos.push_back(imgINFO[ifl].stagePos);
        }

        // Combine centers
        std::vector< std::vector<int> > I0centers;
        std::vector<int> singleI0Cent(2);
        for (ifl=0; ifl<imgINFO.size(); ifl++) {
          imgAddr = I0fileNames[imgINFO[ifl].imgNum];
          auto indI = imgAddr.rfind("/") + 1;
          std::string centFileName = params.centerDir + "temp/centers_"
              + imgAddr.substr(indI, imgAddr.size() - (indI + 4)) + "[2].dat";
          save::importDat<int>(singleI0Cent, centFileName);
          I0centers.push_back(singleI0Cent);
        }
        save::saveDat<int>(I0centers, I0centers_filename);

        // Combine results
        std::vector< std::vector< std::vector<double> > > I0results;
        std::vector< std::vector<double> > singleI0results(4);
        for (uint i=0; i<singleI0results.size(); i++) {
          singleI0results[i].resize(params.I0ellRats.size());
        }
        for (ifl=0; ifl<imgINFO.size(); ifl++) {
          imgAddr = I0fileNames[imgINFO[ifl].imgNum];
          auto indI = imgAddr.rfind("/") + 1;
          std::string resultsFileName = params.centerDir + "temp/results_"
              + imgAddr.substr(indI, imgAddr.size() - (indI + 4)) + "[4,4].dat";
          save::importDat<double>(singleI0results, resultsFileName);
          I0results.push_back(singleI0results);
        }
        save::saveDat<double>(I0results, params.centerDir 
            + "I0-results_run-" + curRun + "_scan-" + to_string(curScan)
            + "_shape[" + to_string(I0results.size()) 
            + "," + to_string(I0results[0].size())
            + "," + to_string(I0results[0][0].size()) + "].dat");

        // Combine I0 norms
        std::vector<double> singleI0norm(1);
        for (ifl=0; ifl<imgINFO.size(); ifl++) {
          imgAddr = I0fileNames[imgINFO[ifl].imgNum];
          auto indI = imgAddr.rfind("/") + 1;
          std::string resultsFileName = params.centerDir + "temp/norm_"
              + imgAddr.substr(indI, imgAddr.size() - (indI + 4)) + "[1].dat";
          save::importDat<double>(singleI0norm, resultsFileName);
          I0norms[ifl] = singleI0norm[0];
        }
        save::saveDat<double>(I0norms, params.centerDir 
            + "I0-norms_run-" + curRun + "_scan-" + to_string(curScan)
            + "_shape[" + to_string(I0norms.size()) + "].dat");
        cout<<"NORM "<<I0norms[0]<<endl;


        // Save stage positions
        std::string fName = params.centerDir
            + "I0-stagePos_run-" + curRun
            + "_shape[" + to_string(I0stagePos.size()) + "].dat";
        if (!tools::fileExists(fName)) {
          save::saveDat<int>(I0stagePos, fName);
        }
      }
      else {
        save::importDat<int>(I0centers, I0centers_filename);
        I0norms.resize(imgINFO.size());
        save::importDat<double>(I0norms, params.centerDir
            + "I0-norms_run-" + curRun + "_scan-" + to_string(curScan)
            + "_shape[" + to_string(I0norms.size()) + "].dat");
        /*  
        std::string fName = params.centerDir
            + "I0-stagePos_run-" + curRun
            + "_shape[" + to_string(I0stagePos.size()) + "].dat";
        save::importDat<int>(I0stagePos, fName);
        save::importDat<double>(I0results, params.centerDir 
            + "I0-results_run-" + curRun + "_scan-" + to_string(curScan)
            + "_shape[" + to_string(I0results.size()) 
            + "," + to_string(I0results[0].size())
            + "," + to_string(I0results[0][0].size()) + "].dat");
        */
      }
    }

  }
  else {
    std::string I0fileName = params.centerDir
        + "avgI0Center_run-" + curRun + "_shape[2].dat";
    std::string rawFileName = params.centerDir
        + "avgRawCenter_run-" + curRun + "_shape[2].dat";
    if (!tools::fileExists(I0fileName) || !tools::fileExists(rawFileName)) {
      int NrawCents = 0;
      float raw_centerR = 0;
      float raw_centerC = 0;
      std::string fileName;
      std::string prefRI("rawImg-centers_run-" + curRun);
      DIR* dir = opendir(params.centerDir.c_str());
      struct dirent* ent;
      while ((ent = readdir(dir)) != NULL) {
        string curFileName(ent->d_name);
        if (prefRI.compare(curFileName.substr(0, prefRI.length())) == 0) {
          fileName = curFileName;
          auto sInd = fileName.find("scan");
          auto fInd = fileName.find("_", sInd);
          int tmpScan = stoi(fileName.substr(sInd+5, fInd - (sInd+5)));
          // Check if scan is bad
          if (std::find(params.badScans.begin(), params.badScans.end(), tmpScan)
              == params.badScans.end()) {
            
            auto badImgItr = params.badImages.find(tmpScan);
            if (badImgItr != params.badImages.end()) {
              if (std::find(badImgItr->second.begin(), badImgItr->second.end(),
                    imgINFO[ifl].stagePos) != badImgItr->second.end()) {
                continue;
              }
            }
 
            save::importDat<int>(centers, params.centerDir 
                + "rawImg-centers_run-" + curRun + "_scan-" + to_string(tmpScan)
                + "_shape[" + to_string(centers.size()) 
                + "," + to_string(centers[0].size()) + "].dat");
            for (ifl=0; ifl<imgINFO.size(); ifl++) {
              raw_centerR += centers[ifl][0];
              raw_centerC += centers[ifl][1];
              NrawCents += 1;
            }
          }
        }
      }
      std::vector<int> avgRawCent(2);
      avgRawCent[0] = tools::round(raw_centerR/NrawCents);
      avgRawCent[1] = tools::round(raw_centerC/NrawCents);
      save::saveDat<int>(avgRawCent, rawFileName);

      int NI0Cents = 0;
      float I0_centerR = 0;
      float I0_centerC = 0;
      prefRI = "I0-results_run-" + curRun;
      std::vector< std::vector< std::vector<double> > > I0results(imgINFO.size());
      for (uint i=0; i<imgINFO.size(); i++) {
        I0results[i].resize(params.I0ellRats.size());
        for (uint ii=0; ii<params.I0ellRats.size(); ii++) {
          I0results[i][ii].resize(4);
        }
      }
      DIR* dirI0 = opendir(params.centerDir.c_str());
      while ((ent = readdir(dirI0)) != NULL) {
        string curFileName(ent->d_name);
        if (prefRI.compare(curFileName.substr(0, prefRI.length())) == 0) {
          fileName = curFileName;
          auto sInd = fileName.find("scan");
          auto fInd = fileName.find("_", sInd);
          int tmpScan = stoi(fileName.substr(sInd+5, fInd - (sInd+5)));
          // Check if scan is bad
          if (std::find(params.badScans.begin(), params.badScans.end(), tmpScan)
              == params.badScans.end()) {
            
            auto badImgItr = params.badImages.find(tmpScan);
            if (badImgItr != params.badImages.end()) {
              if (std::find(badImgItr->second.begin(), badImgItr->second.end(),
                  imgINFO[ifl].stagePos) != badImgItr->second.end()) {
                continue;
              }
            }
            
            save::importDat<double>(I0results, params.centerDir 
                + "I0-results_run-" + curRun + "_scan-" + to_string(tmpScan)
                + "_shape[" + to_string(I0results.size()) 
                + "," + to_string(I0results[0].size())
                + "," + to_string(I0results[0][0].size()) + "].dat");
            for (ifl=0; ifl<imgINFO.size(); ifl++) {
              for (uint i=0; i<params.I0ellRats.size(); i++) {
                if (std::isnan(I0results[ifl][i][0])) {
                  continue;
                }

                I0_centerR += I0results[ifl][i][1];
                I0_centerC += I0results[ifl][i][3];
                NI0Cents += 1;
              }
            }
          }
        }
      }
      std::vector<int> avgI0Cent(2);
      avgI0Cent[0] = tools::round(I0_centerR/NI0Cents);
      avgI0Cent[1] = tools::round(I0_centerC/NI0Cents);
      save::saveDat<int>(avgI0Cent, I0fileName);
    }

    /////  Import average centers and I0 centers to make true center  /////

    std::vector<int> I0avgCent(2);
    std::vector<int> rawAvgCent(2);
    save::importDat<int>(I0avgCent, I0fileName);
    save::importDat<int>(rawAvgCent, rawFileName);
    save::importDat<int>(I0centers, params.centerDir
        + "I0-centers_run-" + curRun + "_scan-" + to_string(curScan)
        + "_shape[" + to_string(centers.size())
        + "," + to_string(centers[0].size()) + "].dat");
    
    for (uint ifl=0; ifl<imgINFO.size(); ifl++) {
      centers[ifl][0] = rawAvgCent[0] + (I0centers[ifl][0] - I0avgCent[0]);
      centers[ifl][1] = rawAvgCent[1] + (I0centers[ifl][1] - I0avgCent[1]);
    }

    // Import I0 norms
    save::importDat<double>(I0norms, params.centerDir
        + "I0-norms_run-" + curRun + "_scan-" + to_string(curScan)
        + "_shape[" + to_string(I0norms.size()) + "].dat");

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

  //////////////////////////////////////
  /////  Making Variables to Save  /////
  //////////////////////////////////////

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

  results_folder = params.preProcOutputDir + "/" + subFolder
        + "/Run-" + curRun + "/";

  if (params.verbose) std::cout << "\n\nINFO: Saving to "
      + results_file_prefix << endl;
  if (!tools::fileExists(results_folder)) {
    mkdir(results_folder.c_str(), 0777);
    mkdir((results_folder + "/radialGrouping").c_str(), 0777);
  }

  results_file_prefix = results_folder + "/Scan-" + to_string(curScan);

  std::vector<int>    res_imgNum;
  std::vector<int>    res_imgIsRef;
  std::vector<int>    res_timeStamp;
  std::vector<int>    res_stagePos;
  std::vector<float>  res_throttle;
  std::vector<int>    res_centerC;
  std::vector<int>    res_centerR;
  std::vector<float>  res_centerCstdRatio;
  std::vector<float>  res_centerRstdRatio;
  std::vector<int>    res_I0centerC;
  std::vector<int>    res_I0centerR;
  std::vector<float>  res_I0norm;
  std::vector<float>  res_imgNorm;
  std::vector<float>  res_readoutNoise;
  std::vector<float>  res_zeroHighQ;

  std::vector<std::vector<std::vector<double> > > res_imgOrig(imgINFO.size());
  std::vector<std::vector<std::vector<int> > > res_imgOrig_nanMap(imgINFO.size());
  std::vector<std::vector<std::vector<double> > > res_imgSubBkg(imgINFO.size());
  std::vector<std::vector<std::vector<int> > > res_imgSubBkg_nanMap(imgINFO.size());
  //std::vector<std::vector<float> > res_legCoeffs;
  std::vector<std::vector<double> > res_rawAzmAvg(imgINFO.size());
  std::vector<std::vector<int> > res_rawAzmAvg_nanMap(imgINFO.size());
  std::vector<std::vector<double> > res_azmAvg(imgINFO.size());
  std::vector<std::vector<int> > res_azmAvg_nanMap(imgINFO.size());
  std::vector<std::vector<double> > res_imgRadSTD(imgINFO.size());


  for (auto const & pv : params.pvMap) {
    pvVals[pv.first].resize(imgINFO.size(), 0);
    pvValsDer[pv.first].resize(imgINFO.size(), 0);
  }



  ///////////////////////////////////////
  /////  Processing images in scan  /////
  ///////////////////////////////////////
 
  for (ifl=0; ifl<imgINFO.size(); ifl++) {
    std::cout << "INFO: processing image: " << imgINFO[ifl].fileName << endl;


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
    res_imgIsRef.push_back(imgIsRef);


    /////  Image number (ordered)  /////
    res_imgNum.push_back(imgINFO[ifl].imgNum);

    /////  Stage position  /////
    res_stagePos.push_back(imgINFO[ifl].stagePos);

    /////  Image capture time  /////
    res_timeStamp.push_back(imgINFO[ifl].time);

    /////  Image radial stds  /////
    res_imgRadSTD.push_back(imgRadSTDs[ifl]);

    /////  I0 image norm  /////
    res_I0norm.push_back(I0norms[ifl]);

    /////  Center  /////
    
    centerR = centers[ifl][0];
    centerC = centers[ifl][1];
    centerRstdRatio = std::fabs(centerR - centerRmean)/centerRstd;
    centerCstdRatio = std::fabs(centerC - centerCmean)/centerCstd;
    res_centerR.push_back(std::round(centers[ifl][0]));
    res_centerC.push_back(std::round(centers[ifl][1]));
    res_centerRstdRatio.push_back(centerRstdRatio);
    res_centerCstdRatio.push_back(centerCstdRatio);

    centerR = std::round(centers[ifl][0]);
    centerC = std::round(centers[ifl][1]);
    if (params.hasI0) {
      I0centerR = I0centers[ifl][0];
      I0centerC = I0centers[ifl][1];
    }
    else if (params.scanAvgCenter) {
      I0centerR = -1;
      I0centerC = -1;
      centerR = std::round(centerRmean);
      centerC = std::round(centerCmean);
    }
    else if (params.scanAvgCenter && params.I0centers) {
      std::cerr << "Must specify only one center finding option: scanAvgCenter or I0centers./n";
      exit(0);
    }
    else {
      std::cerr << "Must specify only one center finding option: scanAvgCenter or I0centers./n";
      exit(0);
    }
    // TODO: write an averaging option


    /*
    if ((centerRstdRatio > params.centerSTDcut) 
        || (centerCstdRatio > params.centerSTDcut)) {

      centerR = centerRmean;
      centerC = centerCmean;
    }
    */
    res_I0centerR.push_back(std::round(I0centers[ifl][0]));
    res_I0centerC.push_back(std::round(I0centers[ifl][1]));




    /////  Throttle  /////
    if (imgINFO[ifl].throttle == -1) {
      res_throttle.push_back(params.throttle);
    }
    else {
      res_throttle.push_back(imgINFO[ifl].throttle);
    }


    /////  Get PV values  /////
    if (params.getPVs) {
      for (auto const & pv : params.pvMap) {
        std::string pvName = pv.first;
        int pvInd = (imgINFO[ifl].time - pvStartTimes[pvName])/
                      params.pvSampleTimes - 2;
        for (int i=0; i<params.imgShutterTime/params.pvSampleTimes; i++) {
          pvVals[pvName][ifl]    += pvAll[pvName][pvInd-i];
          pvValsDer[pvName][ifl] += pvAllDer[pvName][pvInd-i];
        }
        pvVals[pvName][ifl] /= params.imgShutterTime/params.pvSampleTimes;
        pvValsDer[pvName][ifl] /= params.imgShutterTime/params.pvSampleTimes;
      }
    }
 


    ///////  Load image  ///////
    if (params.hasI0) {
      imgAddr = I0fileNames[imgINFO[ifl].imgNum];
      if (params.verbose) cout << "INFO: Trying to open " << imgAddr << "\t .....";
      imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
      if (params.verbose) cout << "\tpassed!\n\n";

      if (params.imgMatType.compare("uint16") == 0) {
        imgI0 = imgProc::getImgVector<uint16_t>(imgMat, 
            imgI0.size(), imgI0.size()/2, imgI0.size()/2, true, NULL);
            //-1, params.holeR, params.holeC);
      }
      else if (params.imgMatType.compare("uint32") == 0) {
        imgI0 = imgProc::getImgVector<uint32_t>(imgMat, 
            imgI0.size(), imgI0.size()/2, imgI0.size()/2, true, NULL);
            //-1, params.holeR, params.holeC);
      }
      else {
        std::cerr << "ERROR: Do not recognize imgMatType = " 
            + params.imgMatType + "!!!\n";
        exit(1);
      }
    }
 
    imgAddr = imgINFO[ifl].path + imgINFO[ifl].fileName;
    if (params.verbose) cout << "INFO: Trying to open " << imgAddr << "\t .....";
    imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
    if (params.verbose) cout << "\tpassed!\n\n";

    if (params.imgMatType.compare("uint16") == 0) {
      imgOrig = imgProc::getImgVector<uint16_t>(imgMat, 
          Nrows, Nrows/2, Ncols/2, true, NULL);
          //params.holeRad, params.holeR, params.holeC);
      imgSubBkg = imgProc::getImgVector<uint16_t>(imgMat, 
          Nrows, Nrows/2, Ncols/2, true, &imgBkg);
          //params.holeRad, params.holeR, params.holeC,
          //&params.nanMap);
    }
    else if (params.imgMatType.compare("uint32") == 0) {
      imgOrig = imgProc::getImgVector<uint32_t>(imgMat, 
          Nrows, Nrows/2, Ncols/2, true, NULL);
          //params.holeRad, params.holeR, params.holeC);
      imgSubBkg = imgProc::getImgVector<uint32_t>(imgMat, 
          Nrows, Nrows/2, Ncols/2, true, &imgBkg);
          //params.holeRad, params.holeR, params.holeC,
          //&params.nanMap);
    }
    else {
      std::cerr << "ERROR: Do not recognize imgMatType = " 
          + params.imgMatType + "!!!\n";
      exit(1);
    }

    // Making new nanMaps //
    cout<<"starting to make the map"<<endl;
    std::vector< std::vector<int> > imgOrig_nanMap(imgOrig.size());
    std::vector< std::vector<int> > imgSubBkg_nanMap(imgOrig.size());
    for (uint ir=0; ir<imgOrig.size(); ir++) {
      imgOrig_nanMap[ir].resize(imgOrig[ir].size());
      imgSubBkg_nanMap[ir].resize(imgOrig[ir].size());
      for (uint ic=0; ic<imgOrig[ir].size(); ic++) {
        imgOrig_nanMap[ir][ic] = params.nanMap[ir][ic];
        imgSubBkg_nanMap[ir][ic] = params.nanMap[ir][ic];
      }
    }
   

    // Remove Xray hits
    cout<<"start removing xray"<<endl;
    imgOrig_nanMap = imgProc::removeXrayHits(
          &imgOrig, imgOrig_nanMap,
          params.XrayHighCut, params.XrayLowCut, 
          params.XraySTDcut, params.XrayWindow);
    cout<<"start removing xray"<<endl;
    imgSubBkg_nanMap = imgProc::removeXrayHits(
          &imgSubBkg, imgSubBkg_nanMap,
          params.XrayHighCut, params.XrayLowCut, 
          params.XraySTDcut, params.XrayWindow,
          xRayHitHistos);
    cout<<"end removing xray"<<endl;

    if (params.pltVerbose) {
      plt.printRC(imgSubBkg, "./plots/imgSubBkg_original");
      std::vector< std::vector<double> > tst;
      for (int ir=0; ir<imgSubBkg_nanMap.size(); ir++) {
        tst[ir].resize(imgSubBkg_nanMap[ir].size());
        for (int ic=0; ic<imgSubBkg_nanMap[ir].size(); ic++) {
          tst[ir][ic] = imgSubBkg_nanMap[ir][ic];
        }
      }
      plt.printRC(tst, "./plots/imgSubBkg_nanMap_original");

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
    }




    /////  Remove pixel outliers  /////
    cout<<"start removing outliers"<<endl;
    /*
    outlierImage = radProc.removeOutliers(
        imgSubBkg, imgSubBkg_nanMap, 
        centerR, centerC, params.imgEdgeBuffer,
        params.NradAzmBins, params.NshellOutlierLoops,
        params.shellWidth, params.Npoly,
        params.stdChangeRatioLeft, params.stdChangeRatioRight,
        params.stdAccRatioLeft, params.stdAccRatioRight,
        params.stdCutLeft, params.stdCutRight, 
        params.fracShellSTDcutLeft, params.fracShellSTDcutRight,
        imgINFO[ifl].stagePos, params.outlierMapSTDcut, 
        true, params.outlierVerbose, NULL, radPixelHistos);
    */

    /*
    outlierImage = radProc.removeOutliers_stdInclude(
        imgSubBkg, imgSubBkg_nanMap, 
        centerR, centerC, params.imgEdgeBuffer,
        params.NradAzmBins, params.NshellOutlierLoops,
        params.shellWidth, params.Npoly, 
        params.stdIncludeLeft, params.distSTDratioLeft,
        params.stdCutLeft, params.fracShellSTDcutLeft,
        params.stdIncludeRight, params.distSTDratioRight,
        params.stdCutRight, params.fracShellSTDcutRight,
        params.stdChangeRatio, imgINFO[ifl].stagePos, 
        params.outlierMapSTDcut, true, params.outlierVerbose, 
        &plt, radPixelHistos);
    */


    cout<<"start removing outliers"<<endl;
    if (params.outlierSimpleSTDcut > 0) {
      outlierImage = radProc.removeOutliersSimple(
          imgSubBkg, imgSubBkg_nanMap,
          centerR, centerC, params.imgEdgeBuffer,
          params.NradAzmBins, params.shellWidth,
          false, params.Npoly, params.outlierSimpleSTDcut, 
          imgINFO[ifl].stagePos, params.outlierMapSTDcut,
          true, (false || params.outlierVerbose), NULL);//, &plt);
    }
    cout<<"end removing outliers"<<endl;
        
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
    readoutNoise = imgProc::removeAvgReadOutNoise(
        imgSubBkg, imgSubBkg_nanMap, centerR, centerC, 
    //readoutNoise = imgProc::removeMedianReadOutNoise(imgSubBkg, centerR, centerC, 
        params.readoutAzmBinStart, params.readoutAzmBinEnd, 
        params.imgEdgeBuffer);

    // Only use this when interested in reference images
    //    background fit fails when subtracting readout
    //    noise on the phosphor screen. Comment previous line.
    //readoutNoise = imgProc::removeReadOutNoise(imgSubBkg, imgSubBkg_nanMap);

    res_readoutNoise.push_back(readoutNoise);

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
            imgLaser[ir][ic] = imgLaser[(imgCut-1)-ir][(imgCut-1)-ic];
          }
        }
      }


      ///  Finding asymmetric parts of the image  ///
      std::vector< std::vector<double> > oddImgReal = imgProc::asymmetrize(imgLaser, 
          centerR, centerC, imgCut, imgCut, 
          oddImgImgn, fftFref, fftIn, fftBref, fftOut);

      // Building map of ratio of noise/"signal" = asymmetric/symmetric
      for (int ir=0; ir<imgCut; ir++) { 
        for (int ic=0; ic<imgCut; ic++) {
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
        imgSubBkg_nanMap[removePairs[ip].first][removePairs[ip].second] = 1;
      }
      ///  Remove laser background from image by setting pixels = NANVAL  ///
      //for (ir=0; ir<imgCut; ir++) {
      //  for (ic=0; ic<imgCut; ic++) {
      //    imgSubBkg[ir][ic] = imgOrig[ir][ic];
      //  }
      //}

      //for (ir=imgCut/2-40; ir<imgCut/2+40; ir++) {
      //  for (ic=imgCut/2-120; ic<imgCut/2; ic++) {
      //    imgSubBkg[ir][ic] = NANVAL;
      //  }
      //}
      //for (ir=imgCut/2-200; ir<imgCut/2+50; ir++) {
      //  for (ic=imgCut-225; ic<imgCut-115; ic++) {
      //    imgSubBkg[ir][ic] = NANVAL;
      //  }
      //}
      //for (ir=imgCut/2-60; ir<imgCut/2+35; ir++) {
      //  for (ic=imgCut-75; ic<imgCut; ic++) {
      //    imgSubBkg[ir][ic] = NANVAL;
      //  }
      //}
      
      if (params.pltVerbose) {
        
        std::vector< std::vector<double> > pltStdRat(imgCut);
        for (int ir=0; ir<imgCut; ir++) {
          pltStdRat[ir].resize(imgCut);
          for (int ic=0; ic<imgCut; ic++) {
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
            "plots/oddImgReal_"+curRun+"_"+to_string(curScan)
            +"_"+to_string(imgINFO[ifl].stagePos));//, pltOpts, pltVals);

        pltOpts[0] = minimum;	pltVals[0] = "0";
        pltOpts[1] = maximum;	pltVals[1] = "5e6";
        delete plt.printRC(imgLaserBkg, 
            "plots/imgLaserBkg_"+curRun+"_"+to_string(curScan)
            +"_"+to_string(imgINFO[ifl].stagePos));

        pltOpts[0] = minimum;	pltVals[0] = "0"; //"1e4";
        pltOpts[1] = maximum;	pltVals[1] = to_string(params.coreValThresh); //"1e6";
        //pltOpts.push_back(logz);  pltVals.push_back("");
        delete plt.printRC(pltStdRat, 
            "plots/imgStdRat_"+curRun+"_"+to_string(curScan)
            +"_"+to_string(imgINFO[ifl].stagePos), pltOpts, pltVals);

        for (int ir=imgCut/2-40; ir<imgCut/2+40; ir++) {
          for (int ic=imgCut/2-120; ic<imgCut/2; ic++) {
            pltStdRat[ir][ic] = 0;
          }
        }
        for (int ir=imgCut/2-200; ir<imgCut/2+50; ir++) {
          for (int ic=imgCut-225; ic<imgCut-115; ic++) {
            pltStdRat[ir][ic] = 0;
          }
        }
        for (int ir=imgCut/2-60; ir<imgCut/2+35; ir++) {
          for (int ic=imgCut-75; ic<imgCut; ic++) {
            pltStdRat[ir][ic] = 0;
          }
        }
        delete plt.printRC(pltStdRat, 
            "plots/imgStdRatRemove_"+curRun+"_"+to_string(curScan)
            +"_"+to_string(imgINFO[ifl].stagePos), pltOpts, pltVals);
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
          "plots/imgOrig_"+curRun+"_"+to_string(curScan)
          +"_"+to_string(imgINFO[ifl].stagePos));
      delete plt.printRC(pltSubBkg, 
          "plots/imgSub_"+curRun+"_"+to_string(curScan)
          +"_"+to_string(imgINFO[ifl].stagePos));
    }


    //// TODO!!!!!!!!!!!!! MUST CENTER IMAGE BEFORE FITTING
    //////  Legendre Fit  //////
    assert(imgCut%5 == 0);
    assert(imgCut%5 == 0);
    const int lgFit_Rows = imgCut/5;
    const int lgFit_Cols = imgCut/5;
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
    std::fill(rawAzmAvg_nanMap.begin(), rawAzmAvg_nanMap.end(), 0);
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

        if (imgSubBkg_nanMap[ir][ic] == 0) {
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
      else rawAzmAvg_nanMap[ir] = 1;
    }

    /////  Image norm  /////
    imgNorm = 0;
    count = 0;
    for (int i=imgNormBinMin; i<imgNormBinMax; i++) {
      if (rawAzmAvg_nanMap[i]) continue;
      imgNorm += rawAzmAvg[i]/atmDiff[i];
      count++;
    }
    imgNorm /= count;

    for (uint ir=0; ir<imgSubBkg.size(); ir++) {
      for (uint ic=0; ic<imgSubBkg[ir].size(); ic++) {
        imgSubBkg[ir][ic] /= imgNorm;
      }
    }

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      rawAzmAvg[iq] /= imgNorm;
      azmAvg[iq] = rawAzmAvg[iq];
      azmAvg_nanMap[iq] = rawAzmAvg_nanMap[iq];
    }

    /////  Plotting azm lines  /////
    cout<<"NORM: "<<imgINFO[ifl].stagePos<<"  "<<imgNorm<<endl;
    std::string dirName = "/reg/ued/ana/scratch/"
        + params.molName + "/" + params.experiment + "/polarLineOutTest/";
    pLO = radProc.getPolarLineOut(&imgSubBkg,
    //pLO = radProc.getPolarLineOut(&imgOrig,
            centerR, centerC, 168, 45, 200);
    cout<<"made it"<<endl;
    for (int k=0; k<pLO.size(); k++) {
      pLO[k] /= imgNorm;
    }
    save::saveDat<double>(pLO, dirName 
        + "scan-" + to_string(curScan)
        + "_stagePos-" + to_string(imgINFO[ifl].stagePos)
        + "_pixRange-168-213_Bins[" + to_string(pLO.size()) + "].dat"); 

    pLO = radProc.getPolarLineOut(&imgSubBkg,
    //pLO = radProc.getPolarLineOut(&imgOrig,
            centerR, centerC, 222, 44, 200);
    for (int k=0; k<pLO.size(); k++) {
      pLO[k] /= imgNorm;
    }
    save::saveDat<double>(pLO, dirName 
        + "scan-" + to_string(curScan)
        + "_stagePos-" + to_string(imgINFO[ifl].stagePos)
        + "_pixRange-222-266_Bins[" + to_string(pLO.size()) + "].dat"); 

    //pLO = radProc.getPolarLineOut(&imgOrig,
    pLO = radProc.getPolarLineOut(&imgSubBkg,
            centerR, centerC, 266, 45, 200);
    for (int k=0; k<pLO.size(); k++) {
      pLO[k] /= imgNorm;
    }
    save::saveDat<double>(pLO, dirName 
        + "scan-" + to_string(curScan)
        + "_stagePos-" + to_string(imgINFO[ifl].stagePos)
        + "_pixRange-266-311_Bins[" + to_string(pLO.size()) + "].dat"); 

    /////  Build References  /////
    if (imgIsRef && params.hasRef) {
      // I0 does not have any nan values in current configuration
      for (int iq=0; iq<params.NradAzmBins; iq++) {
        refAzmAvg[iq] = (refAzmAvg[iq]*Nref + rawAzmAvg[iq])/(Nref + 1);
      }
      Nref++;
      I0refImgs[imgINFO[ifl].stagePos].resize(imgI0.size());
      I0refCount += 1;
      for (uint ir=0; ir<imgI0.size(); ir++) {
        I0refImgs[imgINFO[ifl].stagePos][ir].resize(imgI0[ir].size(), 0);
        for (uint ic=0; ic<imgI0.size(); ic++) {
          imgI0ref[ir][ic] += imgI0[ir][ic];
          I0refImgs[imgINFO[ifl].stagePos][ir][ic] = imgI0[ir][ic];
        }
      }
    }
    else if (!normI0ref) {
      for (uint ir=0; ir<imgI0.size(); ir++) {
        for (uint ic=0; ic<imgI0.size(); ic++) {
          imgI0ref[ir][ic] /= I0refCount;//imgI0refCount[ir][ic];
        }
      }

      for (auto & itr : I0refImgs) {
        for (uint ir=0; ir<imgI0.size(); ir++) {
          for (uint ic=0; ic<imgI0.size(); ic++) {
            itr.second[ir][ic] -= imgI0ref[ir][ic];
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
    if (params.hasI0) {
      save::saveDat<double>(imgI0,
          params.preProcI0OutputDir
          + "data-" + runName
          + "_I0_scan-" + to_string(curScan)
          + "_stagePos-" + to_string(imgINFO[ifl].stagePos)
          + "_bins[" + to_string(imgI0.size())
          + "," + to_string(imgI0[0].size()) + "].dat");

      if (!imgIsRef) {
        for (uint ir=0; ir<imgI0.size(); ir++) {
          for (uint ic=0; ic<imgI0.size(); ic++) {
            imgI0[ir][ic] -= imgI0ref[ir][ic];
          }
        }

        save::saveDat<double>(imgI0,
            params.preProcI0OutputDir
            + "data-" + runName
            + "_I0refSubtracted_scan-" + to_string(curScan)
            + "_stagePos-" + to_string(imgINFO[ifl].stagePos)
            + "_bins[" + to_string(imgI0.size())
            + "," + to_string(imgI0[0].size()) + "].dat");
      }
    }

    /////  Low Pass Filtering  /////
    if (params.verbose) cout << "INFO: Starting Low Pass Filter\n";

    std::vector<double> tstR(filtFFToutSize);
    if (!padRange) {
      while (rawAzmAvg_nanMap[padRange]) {
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


    if (params.verbose) cout << "\tForward FFT\n";
    if (params.pltFilterVerbose) {
      std::vector<double> tst(params.NradAzmBins);
      for (int iq=0; iq<params.NradAzmBins; iq++) {
        tst[iq] = qSpace[iq];
      }
      pltVals[0] = "-2";
      pltVals[1] = "2";
      delete plt.print1d(tst, 
          "./plots/filtInput_" + to_string(imgINFO[ifl].stagePos),
          pltOpts, pltVals);
    }

    fftw_execute(filtFFTf);

    if (params.verbose) cout << "\tBackwards FFT\n";
    for (int ir=0; ir<filtFFToutSize; ir++) {
      rSpace[ir][0] *= bandPassFilter[ir]/sqrt(params.NradAzmBins);
      rSpace[ir][1] *= bandPassFilter[ir]/sqrt(params.NradAzmBins);
    }

    fftw_execute(filtFFTb);

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      if (azmAvg_nanMap[iq] == 0) {
        azmAvg[iq] = qSpace[iq]*atmDiff[iq]/sqrt(params.NradAzmBins);

        if (iq < params.suppressBins) {
          azmAvg[iq] /= pow(sin((PI/2)*iq/params.suppressBins), 6);
        }
        azmAvg[iq] += filtBaseline;

        if (iq < padRange) {
          azmAvg_nanMap[iq] = 1;
        }
      }
    }

    ///  Send high Q to 0 for filtered results  ///
    if (params.verbose) cout << "INFO: Send high Q to 0\n";
    filtBaseline = 0;
    for (int iq=params.readoutAzmBinStart; iq<params.readoutAzmBinEnd; iq++) {
      filtBaseline += azmAvg[iq];
    }
    filtBaseline /= params.readoutAzmBinEnd - params.readoutAzmBinStart;

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      azmAvg[iq] -= filtBaseline;
    }
    for (uint ir=0; ir<imgSubBkg.size(); ir++) {
      for (uint ic=0; ic<imgSubBkg[ir].size(); ic++) {
        imgSubBkg[ir][ic] -= filtBaseline;
      }
    }
    res_zeroHighQ.push_back(filtBaseline);

    ///  Image norm filtered result ///
    double filtImgNorm = 0;
    count = 0;
    for (int i=imgNormBinMin; i<imgNormBinMax; i++) {
      if (azmAvg_nanMap[i]) continue;
      filtImgNorm += azmAvg[i]/atmDiff[i];
      count++;
    }
    filtImgNorm /= count;

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      if (azmAvg_nanMap[iq] == 0) {
        azmAvg[iq] /= filtImgNorm;
      }
    }

    ///  Updating and applying image norm  ///
    imgNorm *= filtImgNorm;
    for (uint ir=0; ir<imgSubBkg.size(); ir++) {
      for (uint ic=0; ic<imgSubBkg[ir].size(); ic++) {
        imgSubBkg[ir][ic] /= filtImgNorm;
      }
    }


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
      plts[0] = plt.plot1d(tst1, "blah_" + to_string(imgINFO[ifl].stagePos),
          pltOpts, pltVals);
      plts[1] = plt.print1d(tst2, "./plots/filtOutp_" + to_string(imgINFO[ifl].stagePos),
          pltOpts, pltVals);
      plts[1]->SetLineColor(4);
      plts[0]->SetMaximum(40);
      plt.print1d(plts, "./plots/filtCompare_" + to_string(imgINFO[ifl].stagePos));
      plts[0]->SetMinimum(1);
      plt.print1d(
          plts, "./plots/filtCompareLogy_" + to_string(imgINFO[ifl].stagePos), 
          logy, "true");
      delete plts[0];
      delete plts[1];
    }

    /////  Plotting results  /////
    if (params.pltVerbose) {
      //plt.print1d(rawAzmAvg, "plots/azimuthalAvg_" + to_string(imgINFO[ifl].imgNum));
      plt.print1d(rawAzmAvg, "tazimuthalAvg_" + to_string(imgINFO[ifl].stagePos),
          logy, "true");
    }
 

    //////////////////////////////////////////////////////////
    /////  Filling image variables and filling the tree  /////
    //////////////////////////////////////////////////////////

    if (params.verbose) cout << "INFO: Filling result vectors\n";

    res_imgNorm.push_back(imgNorm);

    res_imgOrig[ifl].resize(imgOrig.size());
    res_imgOrig_nanMap[ifl].resize(imgOrig.size());
    res_imgSubBkg[ifl].resize(imgOrig.size());
    res_imgSubBkg_nanMap[ifl].resize(imgOrig.size());
    for (int ir=0; ir<(int)imgOrig.size(); ir++) {
      res_imgOrig[ifl][ir].resize(imgOrig[ir].size(), 0);
      res_imgOrig_nanMap[ifl][ir].resize(imgOrig[ir].size(), 0);
      res_imgSubBkg[ifl][ir].resize(imgSubBkg[ir].size(), 0);
      res_imgSubBkg_nanMap[ifl][ir].resize(imgSubBkg[ir].size(), 0);
      for (int ic=0; ic<1024; ic++) {
        res_imgOrig[ifl][ir][ic] = imgOrig[ir][ic];
        res_imgOrig_nanMap[ifl][ir][ic] = imgOrig_nanMap[ir][ic];
        res_imgSubBkg[ifl][ir][ic] = imgSubBkg[ir][ic];
        res_imgSubBkg_nanMap[ifl][ir][ic] = imgSubBkg_nanMap[ir][ic];
      }
    }
    /*
    res_imgSubBkg[ifl].resize(params.imgSize);
    res_imgSubBkg_nanMap[ifl].resize(params.imgSize);
    int offset = params.imgSize/2;
    for (int ir=-1*offset; ir<=offset; ir++) {
      res_imgSubBkg[ifl][ir+offset].resize(imgSubBkg[ir+centerR].size(), 0);
      res_imgSubBkg_nanMap[ifl][ir+offset].resize(imgSubBkg[ir+centerR].size(), 0);
      for (int ic=-1*offset; ic<=offset; ic++) {
        res_imgSubBkg[ifl][ir+offset][ic+offset] = 
            imgSubBkg[ir+centerR][ic+centerC];
        res_imgSubBkg_nanMap[ifl][ir+offset][ic+offset] = 
            imgSubBkg_nanMap[ir+centerR][ic+centerC];

      }
    }
    */

    res_azmAvg[ifl].resize(azmAvg.size());
    res_azmAvg_nanMap[ifl].resize(azmAvg.size());
    res_rawAzmAvg[ifl].resize(rawAzmAvg.size());
    res_rawAzmAvg_nanMap[ifl].resize(rawAzmAvg.size());
    for (uint ir=0; ir<rawAzmAvg.size(); ir++) {
      res_azmAvg[ifl][ir] = azmAvg[ir];
      res_azmAvg_nanMap[ifl][ir] = azmAvg_nanMap[ir];
      res_rawAzmAvg[ifl][ir] = rawAzmAvg[ir];
      res_rawAzmAvg_nanMap[ifl][ir] = rawAzmAvg_nanMap[ir];
    }

    
   
    /*
    /////  Saving Indices/Values at Specific Radii  /////
    std::map<int, std::vector< std::vector<int> > > img_inds;
    std::map<int, std::vector< std::vector<int> > > rel_inds;
    std::map<int, std::vector<double> > img_vals;
    std::map<int, std::vector<int> >    img_nans;
    for (int ir=0; ir<(int)imgSubBkg.size(); ir++) {
      if (ir < params.imgEdgeBuffer 
          || imgSubBkg.size() - ir < params.imgEdgeBuffer) continue;
      for (int ic=0; ic<(int)imgSubBkg[ir].size(); ic++) {
        if (ic < params.imgEdgeBuffer 
            || imgSubBkg[ir].size() - ic < params.imgEdgeBuffer) continue;

        radInd = std::round(std::sqrt(
            std::pow(ir-centerR,2) + std::pow(ic-centerC,2)));

        std::vector<int> iinds(2);
        iinds[0] = ir; iinds[1] = ic;
        img_inds[radInd].push_back(iinds);
        
        std::vector<int> rinds(2);
        rinds[0] = ir - centerR; rinds[1] = ic - centerC;
        rel_inds[radInd].push_back(rinds);

        img_vals[radInd].push_back(imgSubBkg[ir][ic]);
        img_nans[radInd].push_back(imgSubBkg_nanMap[ir][ic]);
      }
    }

    radial_base_folder = results_folder + "/radialGrouping/Scan-"
      + to_string(curScan);
  
    if (!tools::fileExists(radial_base_folder)) {
      mkdir(radial_base_folder.c_str(), 0777);
      mkdir((radial_base_folder+"/indices").c_str(), 0777);
      mkdir((radial_base_folder+"/values").c_str(), 0777);
    }
      + "imgNum-" + to_string(imgINFO[ifl].imgNum) + "_";
    for (auto& itr : img_vals) {
      save::saveDat<int>(img_inds[itr.first],
          radial_base_folder + "/indices/" 
          + "imgNum-" + to_string(imgINFO[ifl].imgNum) + "_"
          + "rad-" + to_string(itr.first) + "_imgInds_Shape["
          + to_string(img_inds[itr.first].size()) + ",2].dat");
      save::saveDat<int>(rel_inds[itr.first],
          radial_base_folder + "/indices/" 
          + "imgNum-" + to_string(imgINFO[ifl].imgNum) + "_"
          + "rad-" + to_string(itr.first) + "_relInds_Shape["
          + to_string(rel_inds[itr.first].size()) + ",2].dat");
      save::saveDat<double>(img_vals[itr.first],
          radial_base_folder + "/values/"
          + "imgNum-" + to_string(imgINFO[ifl].imgNum) + "_"
          + "rad-" + to_string(itr.first) + "_imgVals_Shape["
          + to_string(img_nans[itr.first].size()) + "].dat");
      save::saveDat<int>(img_nans[itr.first],
          radial_base_folder + "/values/"
          + "imgNum-" + to_string(imgINFO[ifl].imgNum) + "_"
          + "rad-" + to_string(itr.first) + "_imgNans_Shape["
          + to_string(img_nans[itr.first].size()) + "].dat");
    }
    */
  }


  ////////////////////////////
  /////  Saving Results  /////
  ////////////////////////////


  results_file_prefix = results_folder + "/Scan-" + to_string(curScan) + "_";

  save::saveDat<int>(res_imgNum,
      results_file_prefix +"imgNum_Shape["
      + to_string(res_imgNum.size()) + "].dat");
  save::saveDat<int>(res_imgIsRef,
      results_file_prefix +"imgIsRef_Shape["
      + to_string(res_imgIsRef.size()) + "].dat");
  save::saveDat<int>(res_timeStamp,
      results_file_prefix +"timeStamp_Shape["
      + to_string(res_timeStamp.size()) + "].dat");
  save::saveDat<int>(res_stagePos,
      results_file_prefix +"stagePos_Shape["
      + to_string(res_stagePos.size()) + "].dat");
  save::saveDat<float>(res_throttle,
      results_file_prefix +"throttle_Shape["
      + to_string(res_throttle.size()) + "].dat");
  save::saveDat<int>(res_centerC,
      results_file_prefix +"centerC_Shape["
      + to_string(res_centerC.size()) + "].dat");
  save::saveDat<int>(res_centerR,
      results_file_prefix +"centerR_Shape["
      + to_string(res_centerR.size()) + "].dat");
  save::saveDat<float>(res_centerCstdRatio,
      results_file_prefix +"centerCstdRatio_Shape["
      + to_string(res_centerCstdRatio.size()) + "].dat");
  save::saveDat<float>(res_centerRstdRatio,
      results_file_prefix +"centerRstdRatio_Shape["
      + to_string(res_centerRstdRatio.size()) + "].dat");
  save::saveDat<int>(res_I0centerC,
      results_file_prefix +"I0centerC_Shape["
      + to_string(res_I0centerC.size()) + "].dat");
  save::saveDat<int>(res_I0centerR,
      results_file_prefix +"I0centerR_Shape["
      + to_string(res_I0centerR.size()) + "].dat");
  save::saveDat<float>(res_I0norm,
      results_file_prefix +"I0norm_Shape["
      + to_string(res_I0norm.size()) + "].dat");
  save::saveDat<float>(res_imgNorm,
      results_file_prefix +"imgNorm_Shape["
      + to_string(res_imgNorm.size()) + "].dat");
  save::saveDat<float>(res_readoutNoise,
      results_file_prefix +"readoutNoise_Shape["
      + to_string(res_readoutNoise.size()) + "].dat");
  save::saveDat<float>(res_zeroHighQ,
      results_file_prefix +"zeroHighQ_Shape["
      + to_string(res_zeroHighQ.size()) + "].dat");

  save::saveDat<double>(res_imgOrig,
      results_file_prefix +"imgOrig_Shape["
      + to_string(res_imgOrig.size()) + ","
      + to_string(res_imgOrig[0].size()) + ","
      + to_string(res_imgOrig[0][0].size()) + "].dat");
  save::saveDat<int>(res_imgOrig_nanMap,
      results_file_prefix +"imgOrig_nanMap_Shape["
      + to_string(res_imgOrig.size()) + ","
      + to_string(res_imgOrig[0].size()) + ","
      + to_string(res_imgOrig[0][0].size()) + "].dat");
  save::saveDat<double>(res_imgSubBkg,
      results_file_prefix +"procImg_Shape["
      + to_string(res_imgSubBkg.size()) + ","
      + to_string(res_imgSubBkg[0].size()) + ","
      + to_string(res_imgSubBkg[0][0].size()) + "].dat");
  save::saveDat<int>(res_imgSubBkg_nanMap,
      results_file_prefix +"procImg_nanMap_Shape["
      + to_string(res_imgSubBkg.size()) + ","
      + to_string(res_imgSubBkg[0].size()) + ","
      + to_string(res_imgSubBkg[0][0].size()) + "].dat");

  //std::vector<std::vector<float> > res_legCoeffs;
  save::saveDat<double>(res_rawAzmAvg,
      results_file_prefix +"rawAzmAvg_Shape["
      + to_string(res_rawAzmAvg.size()) + ","
      + to_string(res_rawAzmAvg[0].size()) + "].dat");
  save::saveDat<int>(res_rawAzmAvg_nanMap,
      results_file_prefix +"rawAzmAvg_nanMap_Shape["
      + to_string(res_rawAzmAvg.size()) + ","
      + to_string(res_rawAzmAvg[0].size()) + "].dat");
  save::saveDat<double>(res_azmAvg,
      results_file_prefix +"azmAvg_Shape["
      + to_string(res_rawAzmAvg.size()) + ","
      + to_string(res_rawAzmAvg[0].size()) + "].dat");
  save::saveDat<int>(res_azmAvg_nanMap,
      results_file_prefix +"azmAvg_nanMap_Shape["
      + to_string(res_rawAzmAvg.size()) + ","
      + to_string(res_rawAzmAvg[0].size()) + "].dat");

  save::saveDat<double>(res_imgRadSTD,
      results_file_prefix +"imgRadSTD_Shape["
      + to_string(res_imgRadSTD.size()) + ","
      + to_string(res_imgRadSTD[0].size()) + "].dat");




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

                                        
