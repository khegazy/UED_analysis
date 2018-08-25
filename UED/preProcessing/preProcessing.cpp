#include "preProcessing.h"
#include "../../parameters.h"

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
    if (iPos == string::npos) {
      std::cout << "INFO: Will do background!!!\n";
      iPos = runListName.find("Background") + 11;
      doBackground = true;
    }
    else {
      iPos += 4;
    }
    runName = runListName.substr(iPos, 
        runListName.find("_Scan-") - iPos);
  }

  cout<<"RunName: "<<runName<<endl;
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
  cout<<imgSize<<"  "<<params.imgSize<<endl;
  if (imgSize != params.imgSize) {
    cerr << "ERROR: parameter imgSize does not match with parameter class!!!\n";
    exit(0);
  }


   // Make sure the image has odd number of bins
  if (!(imgSize%2) || !(imgSize%2)) {
    cerr << "ERROR: imgSize and imgSize must be an odd number!!!" << endl;
    exit(0);
  }

  imgProc::radProcTool radProc(params.indicesPath);
  
  uint ifl, iffl, runScanStartItr;
  std::vector<PLOToptions> pltOpts(2);
  std::vector<string> pltVals(2);


  std::vector<double> legCoeffs(params.NradLegBins);
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
  int centerC, centerR, Ncents;

  double count;

  int imgIsRef;
  int Nrefs = 2;

  float pressure;
  float pressureDer;

  int timeStamp;

  int radInd;
  int imgNum;
  int curScan;
  std::string fileName;
  std::string line;
  size_t ipos, fpos;
  std::string date, scan, curRun, rFileName;
  int32_t stagePos;
  float t0SP, t0Time;
  std::vector<imgProc::imgInfoStruct> imgINFO;
  float readoutNoise;

  TFile* file=NULL;
  TTree* tree=NULL;

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
  for (uint ir=0; ir<Nrows; ir++) {
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
    if (params.verbose) std::cout << "INFO: Making background!!!\n";

    double count = 0;
    std::vector< std::vector<double> > imgVec;
    for (ifl=0; ifl<imgINFO.size(); ifl++) {
      ///  Get image  ///
      imgAddr = imgINFO[ifl].path + imgINFO[ifl].fileName;
      imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
      cout<<"image addr: "<<imgAddr<<endl;
      imgProc::threshold(imgMat, params.hotPixel);

      Nrows = imgMat.rows;
      Ncols = imgMat.cols;
      if (params.imgMatType.compare("uint16") == 0) {
        imgVec = imgProc::getImgVector<uint16_t>(imgMat, 
                    Nrows, Nrows/2, Ncols/2); 
      }  
      else if (params.imgMatType.compare("uint32") == 0) {
        imgVec = imgProc::getImgVector<uint32_t>(imgMat, 
                    Nrows, Nrows/2, Ncols/2);
      }  
      else {
        cerr << "ERROR: Do not recognize imgMatType = " 
            + params.imgMatType + "!!!\n";
        exit(0);
      }

      ///  Add image  ///
      if (params.verbose) std::cout << "\tAdding images\n";
      for (uint ir=0; ir<imgVec.size(); ir++) {
        std::transform(imgVec[ir].begin(), imgVec[ir].end(), 
            imgBkg[ir].begin(), imgBkg[ir].begin(),
            std::plus<double>());
      }
      count++;
    }

    ///  Scaling vector  ///
    for (uint ir=0; ir<imgVec.size(); ir++) {
      tools::vScale<double>(imgBkg[ir], 1.0/count);
    }

    ///  Save image  ///
    if (params.verbose) std::cout << "\tSaving image\n";
    save::saveDat<double>(imgBkg, params.backgroundFolder 
        + "/backgroundImg-" + imgINFO[0].run + ".dat");

    if (params.pltVerbose) {
      delete plt.printRC(imgBkg, "plots/background-" + imgINFO[0].run);
    }

    exit(1);
  }


  //////////////////////////////////
  /////  Get background image  /////
  //////////////////////////////////

  cout<<"getting background"<<endl;
  if (params.backgroundImage.compare("NULL") != 0) {
    cout<<"fileName: "<<params.backgroundFolder+ "/" + params.backgroundImage<<endl;
    save::importDat<double>(imgBkg, params.backgroundFolder 
        + "/" + params.backgroundImage);
    if (params.pltVerbose) {
      delete plt.printRC(imgBkg, "plots/importBkg-" + imgINFO[0].run);
    }
  }
  cout<<"done"<<endl;


  ///////////////////////
  /////  Import PV  /////
  ///////////////////////

  long int pvStartTime = 0;
  ipos = params.pressureFileName.find("_");
  for (int i=1; i<stoi(params.pressureFileName.substr(ipos+1, 2));i++) {
    pvStartTime += ppFunct::monthLengths[i]*24*3600;
  }
  ipos = params.pressureFileName.find("_", ipos+1);
  pvStartTime += stoi(params.pressureFileName.substr(ipos+1, 2))*24*3600;
  ipos = params.pressureFileName.find("_", ipos+1);
  ipos = params.pressureFileName.find("_", ipos+1);
  pvStartTime += stoi(params.pressureFileName.substr(ipos+1, 2))*3600;
  ipos = params.pressureFileName.find("_", ipos+1);
  pvStartTime += stoi(params.pressureFileName.substr(ipos+1, 2))*60;
  ipos = params.pressureFileName.find("-", ipos+1);
  pvStartTime += stoi(params.pressureFileName.substr(ipos+1, 2));

  ipos = params.pressureFileName.find("-", ipos+1);
  int pvSize = stoi(params.pressureFileName.substr(ipos+1, params.pressureFileName.length()-5-ipos));

  std::vector<double> pvOrigPressure(pvSize); 
  save::importDat(pvOrigPressure, params.pressureFileName);

  std::vector<double> pvPressure(pvSize, 0);
  std::vector<double> pvPressureDer(pvSize, 0);
  float norm;
  int NsmearSteps = 3*params.pressureSmear/params.pvSampleTimes + 1;
  for (int i=0; i<pvSize; i++) {
    pvPressure[i] = pvOrigPressure[i];
    norm = 1;
    for (int ii=1; ii<NsmearSteps; ii++) {
      if (i+ii < pvSize) {
        pvPressure[i] += pvOrigPressure[i+ii]
                         *exp(-1*std::pow(params.pvSampleTimes*ii, 2)
                             /(2*pow(params.pressureSmear, 2)));
        norm += exp(-1*std::pow(params.pvSampleTimes*ii, 2)
                    /(2*std::pow(params.pressureSmear, 2)));
      }
      if (i-ii >= 0) {
        pvPressure[i] += pvOrigPressure[i-ii]
                         *exp(-1*std::pow(params.pvSampleTimes*ii, 2)
                             /(2*std::pow(params.pressureSmear, 2)));
        norm += exp(-1*std::pow(params.pvSampleTimes*ii, 2)
                    /(2*std::pow(params.pressureSmear, 2)));
      }
    }
    pvPressure[i] /= norm;
  }

  for (int i=0; i<pvSize; i++) {
    if ((i == 0) || (i == pvSize-1)) {
      pvPressureDer[i] = 0;
    }
    else {
      pvPressureDer[i] = ((pvPressure[i+1] - pvPressure[i-1])
                          /(2*params.pvSampleTimes));
    }
  }

  //plt.print1d(pvPressure, "./plots/pressureSmooth");
  //plt.print1d(pvPressureDer, "./plots/pressureDer");
  //plt.print1d(pvOrigPressure, "./plots/pressure");
  //exit(0);



  //////////////////////////////
  /////  Making root file  /////
  //////////////////////////////

   //rFileName = "/reg/ued/ana/scratch/nitroBenzene/rootFiles/" + dataType 
  std::string subFolder;
  if (runListName.find("UVpump", 0) != std::string::npos) {
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

  rFileName =  params.preProcOutputDir + "/" + subFolder + "/"
  //rFileName =  "testing/" + subFolder + "/"
        + "Run-" + curRun + "_"
        + "Scan-" + to_string(curScan) + ".root";

  if (params.verbose) std::cout << "INFO: Making file " << rFileName << endl;

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
  tree->Branch("pressure",	&pressure,	"pressure/F");
  tree->Branch("pressureDer",	&pressureDer,	"pressureDer/F");
  tree->Branch("centerC", 	&centerC, 	"centerC/I");
  tree->Branch("centerR", 	&centerR, 	"centerR/I");
  tree->Branch("imgNorm", 	&imgNorm, 	"imgNorm/F");
  tree->Branch("readoutNoise", 	&readoutNoise, 	"readoutNoise/F");
  tree->Branch("imgOrig", 	&imgOrig);
  tree->Branch("imgSubBkg",     &imgSubBkg);
  tree->Branch("legCoeffs",     &legCoeffs);
  tree->Branch("azmAvg",        &azmAvg);
  if (params.verbose) cout << "INFO: Tree and file are setup!\n\n";


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

  if (params.verbose) cout << "INFO: Start COM center finding\n"; 
  std::vector<int> centersCOM = imgProc::centerSearchCOM(imgINFO,
      params.hotPixel, params.sigma,
      params.blockCentR, params.blockCentC, params.minRad, params.maxRad,
      params.meanInd, params.COMstdScale, params.verbose, NULL); //&plt);

  if (params.verbose) std::cout << "COMcenters: " 
    << centersCOM[0] << "  " << centersCOM[1] << endl;


  /////////////////////////////////
  /////  Fine center finding  /////
  /////////////////////////////////

  //centerR = centersCOM[0];
  //centerC = centersCOM[1];
  ///// Finding image center /////
  if (params.verbose) cout << "INFO: Starting fine center Finding!\n\n";
  centerR = centerC = Ncents = 0;
  for (int icnt=0; icnt<std::min(params.NavgCenters, (int)imgINFO.size()); icnt++) {

    /// Filling image std::vector ///
    imgAddr = imgINFO[icnt].path + imgINFO[icnt].fileName;
    imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
    cout<<"centering image addr: "<<imgAddr<<endl;
    imgProc::threshold(imgMat, params.hotPixel);

    Nrows = imgMat.rows;
    Ncols = imgMat.cols;
    if (params.imgMatType.compare("uint16") == 0) {
      imgCent = imgProc::getImgVector<uint16_t>(imgMat, 
                  Nrows, Nrows/2, Ncols/2, &imgBkg,
                  params.holeRad, params.holeR, params.holeC);
    }  
    else if (params.imgMatType.compare("uint32") == 0) {
      imgCent = imgProc::getImgVector<uint32_t>(imgMat, 
                  Nrows, Nrows/2, Ncols/2, &imgBkg,
                  params.holeRad, params.holeR, params.holeC);
    }  
    else {
      cerr << "ERROR: Do not recognize imgMatType = " 
          + params.imgMatType + "!!!\n";
      exit(0);
    }

    //  Remove pixel outliers
    radProc.removeOutliers(imgCent,  
        centersCOM[0], centersCOM[1], params.buffer,
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


  ///////////////////////////////////////
  /////  Processing images in scan  /////
  ///////////////////////////////////////
 
  for (ifl=0; ifl<imgINFO.size(); ifl++) {

    //if (imgINFO[ifl].stagePos != 1543400) continue; 

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

    /////  Get chamber pressure  /////
    pressure = 0;
    pressureDer = 0;
    int pvInd = (imgINFO[ifl].time - pvStartTime)/params.pvSampleTimes - 2;
    for (int i=0; i<params.imgShutterTime/params.pvSampleTimes; i++) {
      pressure += pvPressure[pvInd-i];
      pressureDer += pvPressureDer[pvInd-i];
    }
    pressure /= params.imgShutterTime/params.pvSampleTimes;
    pressureDer /= params.imgShutterTime/params.pvSampleTimes;
 
    ///////  Load image  ///////
    imgAddr = imgINFO[ifl].path + imgINFO[ifl].fileName;
    if (params.verbose) cout << "INFO: Trying to open " << imgAddr << "\t .....";
    imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
    imgProc::threshold(imgMat, params.hotPixel);
    if (params.verbose) cout << "\tpassed!\n\n";

    if (params.imgMatType.compare("uint16") == 0) {
      imgOrig = imgProc::getImgVector<uint16_t>(imgMat, 
          Nrows, Nrows/2, Ncols/2, NULL,
          params.holeRad, params.holeR, params.holeC);
      imgSubBkg = imgProc::getImgVector<uint16_t>(imgMat, 
          Nrows, Nrows/2, Ncols/2, &imgBkg,
          params.holeRad, params.holeR, params.holeC);
    }
    else if (params.imgMatType.compare("uint32") == 0) {
      imgOrig = imgProc::getImgVector<uint32_t>(imgMat, 
          Nrows, Nrows/2, Ncols/2, NULL,
          params.holeRad, params.holeR, params.holeC);
      imgSubBkg = imgProc::getImgVector<uint32_t>(imgMat, 
          Nrows, Nrows/2, Ncols/2, &imgBkg,
          params.holeRad, params.holeR, params.holeC);
    }
    else {
      cerr << "ERROR: Do not recognize imgMatType = " 
          + params.imgMatType + "!!!\n";
      exit(0);
    }
    //cout<<"444"<<endl;

    /////  Remove pixel outliers  /////
    cout<<"IMGORIG  "<<stagePos<<endl;
    radProc.removeOutliers(imgOrig,  
        centerR, centerC, params.buffer, 
        params.NradAzmBins, params.shellWidth, params.Npoly,
        //imgOrig.size()/2, imgOrig.size()/2, params.shellWidth,
        params.stdIncludeLeft, params.distSTDratioLeft,
        params.stdCutLeft, params.meanBinSize,
        params.stdIncludeRight, params.distSTDratioRight,
        params.stdChangeRatio, params.stdCutRight,
        imgINFO[ifl].stagePos, params.outlierMapSTDcut,
        false, (false || params.outlierVerbose), NULL);

    //cout<<"555"<<endl;
    cout<<"IMGSUB "<<stagePos<<endl;
    outlierImage = radProc.removeOutliers(imgSubBkg,  
        centerR, centerC, params.buffer,
        params.NradAzmBins, params.shellWidth, params.Npoly,
        //imgOrig.size()/2, imgOrig.size()/2, params.shellWidth,
        params.stdIncludeLeft, params.distSTDratioLeft,
        params.stdCutLeft, params.meanBinSize,
        params.stdIncludeRight, params.distSTDratioRight,
        params.stdChangeRatio, params.stdCutRight,
        imgINFO[ifl].stagePos, params.outlierMapSTDcut,
        true, (false || params.outlierVerbose), NULL);//, &plt);
    
    //cout<<"666"<<endl;
    if (params.pltVerbose) {
      //plt.printRC(outlierImage, "outlierSTD_" + imgINFO[ifl].run + "-" + to_string(imgINFO[ifl].scan) + "-" + to_string(imgINFO[ifl].stagePos));
      save::saveDat<double>(outlierImage, 
            "plots/data/outlierSTD-" + runName
                + "-" + to_string(imgINFO[ifl].scan)
                + "-" + to_string(imgINFO[ifl].stagePos) + ".dat");
    }

    ///  Find large clusters corresponding to laser spots  ///
    //std::vector< std::pair<uint, uint> > removePairs;
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

    if (params.pltVerbose) {
      //plt.printRC(outlierBkg, "outlierBKG_" + imgINFO[ifl].run + "-" + to_string(imgINFO[ifl].scan) + "-" + to_string(imgINFO[ifl].stagePos));
      save::saveDat<double>(outlierBkg, 
            "plots/data/outlierBackground-" + runName
                + "-" + to_string(imgINFO[ifl].scan)
                + "-" + to_string(imgINFO[ifl].stagePos) + ".dat");
      for (int ir=0; ir<imgSubBkg.size(); ir ++) {
        for (int ic=0; ic<imgSubBkg[ir].size(); ic ++) {
          if (sqrt(pow(ir-params.holeR,2) + pow(ic-params.holeC,2)) < params.holeRad) continue;
          if (imgSubBkg[ir][ic] == NANVAL) imgSubBkg[ir][ic] = params.hotPixel;
        }
      }

      save::saveDat<double>(imgSubBkg, 
            "plots/data/image-" + runName
                + "-" + to_string(imgINFO[ifl].scan)
                + "-" + to_string(imgINFO[ifl].stagePos) 
                + "[" + to_string(imgSubBkg.size())
                + "," + to_string(imgSubBkg.size()) + "].dat");
    }


    ///  Remove laser background from image by setting pixels = NANVAL  ///
    for (uint ip=0; ip<removePairs.size(); ip++) {
      imgSubBkg[removePairs[ip].first][removePairs[ip].second] = NANVAL;
    }


    ///// Remove readout noise  /////
    if (params.verbose) std::cout << "INFO: Readout noise subtraction.\n";
    /*
    imgProc::removeReadOutNoise(imgOrig);
    readoutNoise = imgProc::removeReadOutNoise(imgSubBkg);
    */
    if (params.hasLaserBkg) {
      readoutNoise = imgProc::removeAvgReadOutNoise(imgSubBkg, centerR, centerC, 
                        0.94*params.NradAzmBins, params.NradAzmBins, 
                        params.buffer, &params.nanMap);
    }
    else {
      readoutNoise = imgProc::removeAvgReadOutNoise(imgSubBkg, centerR, centerC, 
                        0.94*params.NradAzmBins, params.NradAzmBins,
                        params.buffer);
    }


    //////////////////////////////////////////////////////
    /////  Finding and subtracting laser background  /////
    //////////////////////////////////////////////////////
    // We cluster background spots based on the assymetry of 
    //    the image. 

    if (params.hasLaserBkg) {
      if (params.verbose) std::cout << "INFO: Laser background removal.\n";
      // Fill in hole based on symmetry
      std::vector< std::vector<double> > imgLaser = imgSubBkg;
      /*
      for (int ir=centerR-params.holeRad-3; 
          ir<=centerR+params.holeRad+3; ir++) {
        for (int ic=centerC-params.holeRad-3; 
            ic<=centerC+params.holeRad+3; ic++) {
          if (std::pow(ir-centerR,2) + std::pow(ic-centerC,2) 
              < std::pow(params.holeRad, 2)) {
            imgOrig[ir][ic] = imgOrig[(imgSize-1)-ir][(imgSize-1)-ic];
          }
        }
      }
      */

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
      /*
      ///  Remove laser background from image by setting pixels = NANVAL  ///
      for (ir=0; ir<imgSize; ir++) {
        for (ic=0; ic<imgSize; ic++) {
          imgSubBkg[ir][ic] = imgOrig[ir][ic];
        }
      }

      for (ir=imgSize/2-40; ir<imgSize/2+40; ir++) {
        for (ic=imgSize/2-120; ic<imgSize/2; ic++) {
          imgSubBkg[ir][ic] = NANVAL;
        }
      }
      for (ir=imgSize/2-200; ir<imgSize/2+50; ir++) {
        for (ic=imgSize-225; ic<imgSize-115; ic++) {
          imgSubBkg[ir][ic] = NANVAL;
        }
      }
      for (ir=imgSize/2-60; ir<imgSize/2+35; ir++) {
        for (ic=imgSize-75; ic<imgSize; ic++) {
          imgSubBkg[ir][ic] = NANVAL;
        }
      }
      */
      
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

    std::fill(azmAvg.begin(), azmAvg.end(), 0);
    std::fill(azmCounts.begin(), azmCounts.end(), 0);
    for (int ir=0; ir<imgSubBkg.size(); ir++) {
      if (ir < params.buffer || imgSubBkg.size() - ir < params.buffer) continue;
      for (int ic=0; ic<imgSubBkg[ir].size(); ic++) {
        if (ic < params.buffer || imgSubBkg[ir].size() - ic < params.buffer) continue;

        if (imgSubBkg[ir][ic] != NANVAL) {
          radInd = (int)sqrt(pow(ir-centerR,2) + pow(ic-centerC,2));
          if (radInd < params.NradAzmBins) {
            azmAvg[radInd] += imgSubBkg[ir][ic];
            azmCounts[radInd] += 1;
          }
        }
        else {
          //cout<<"skipping nan"<<endl;
        }
      }
    }

    for (uint ir=0; ir<azmAvg.size(); ir++) {
      if (azmCounts[ir] != 0) {
        azmAvg[ir] /= azmCounts[ir];
      }
      else azmAvg[ir] = NANVAL;
    }

    azmAvg.resize(params.NradAzmBins, NANVAL);


    /////  Image norm  /////
    imgNorm = 0;
    count = 0;
    for (int i=0; i<params.NradAzmBins; i++) {
      if (azmAvg[i] == NANVAL) continue;
      imgNorm += azmAvg[i];
      count++;
    }
    imgNorm /= count;


    /////  Plotting results  /////
    if (params.pltVerbose) {
      std::vector<double> test(NradLegBins);
      for (int i=0; i<params.Nlegendres; i++) {
        for (int j=0; j<NradLegBins; j++) {
          test[j] = legCoeffs[i*NradLegBins + j];
        }
        plt.print1d(test, "plots/testLeg" + to_string(i) + "_" + to_string(imgINFO[ifl].imgNum));
      }

      plt.print1d(azmAvg, "plots/azimuthalAvg_" + to_string(imgINFO[ifl].imgNum));
    }
    
    // Save image and info to tree
    tree->Fill();
  }


  tree->Write();
  file->Close();
  cout<<endl<<endl<<endl;

  // Release fftw memory
  fftw_destroy_plan(fftFref);
  fftw_destroy_plan(fftBref);
  fftw_free(fftIn);
  fftw_free(fftOut);

return 1;
}

                                        
