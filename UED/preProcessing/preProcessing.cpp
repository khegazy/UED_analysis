#include "preProcessing.h"
#include "../parameters.h"
#include "config.h"

using namespace std;

class centerfnctr {

    public:
        int symVleg;
        std::vector< std::vector<double> >* img;
        int count;
        double nanVal;
        int minRadBin;
        int centShellW;
        double operator() (std::vector<double> vect) {
          count++;
          if (symVleg == 0)  return imgProc::centerSymXsqr(img, vect[0], vect[1], 8, 1, minRadBin, centShellW, true, nanVal);
          if (symVleg == 1)  return imgProc::centerOddLeg(img, vect[0], vect[1], 1, minRadBin, centShellW, true, nanVal);
          else {
            cerr << "ERROR: Did not select a correct center finding alogrithm, now exiting!!!\n\n";
            exit(0);
          }
          //N2 return imgProc::centerOddLeg(img, vect[0], vect[1], 1, 70, 40, true, params.nanVal);
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


  experimentParameters  expP = experimentParameters();
  std::string runName;
  if (doRunLists) {
    runName = "doRunLists";
  }
  else {
    auto iPos = runListName.find("Run") + 4;
    runName = runListName.substr(iPos, 
        runListName.find("_Scan-") - iPos);
  }

  cout<<"RunName: "<<runName<<endl;
  parameterClass params(runName);
  PLOTclass plt;



  const int Nlegendres = 1;
  const int NradBins = 50;

  // Indices
  //const int imgSize = 935;
  const int imgSize = 895;
  const int Npix = std::pow(imgSize, 2);

  if (Nlegendres != params.Nlegendres) {
    cerr << "ERROR: parameter Nlegendres does not match with parameter class!!!\n";
    exit(0);
  }
  if (NradBins != params.NradBins) {
    cerr << "ERROR: parameter NradBinss does not match with parameter class!!!\n";
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
  uint ifl, iffl, runScanStartItr;
  std::vector<PLOToptions> pltOpts(2);
  std::vector<string> pltVals(2);


  std::vector<double> legCoeffs(params.NradBins);
  std::vector<double> legCounts(params.NradBins);

  std::vector< std::vector<double> > imgOrig(imgSize);
  std::vector< std::vector<double> > imgSubBkg(imgSize);
  std::vector< std::vector<double> > imgTemp(imgSize);
  std::vector< std::vector<double> > imgBkg(imgSize);
  std::vector< std::vector<double> > imgCent;
  float imgNorm;
  for (uint ir=0; ir<imgSize; ir++) {
    imgOrig[ir].resize(imgSize, 0);
    imgSubBkg[ir].resize(imgSize, 0);
    imgTemp[ir].resize(imgSize, 0);
    imgBkg[ir].resize(imgSize, 0);
  }

  std::vector< std::vector<double> > oddImgImgn;
  std::vector< std::vector<double> > symImg(imgSize);
  std::vector< std::vector<double> > stdRatioImg(imgSize);
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

  int Nrows, Ncols;
  double count;

  int imgIsRef;
  int Nrefs = 2;

  int imgNum;
  int curScan;
  std::string fileName;
  std::string line;
  size_t ipos, fpos;
  std::string date, scan, curRun, rFileName;
  int32_t stagePos;
  float t0SP, t0Time;
  std::vector<imgProc::imgInfoStruct> imgINFO;

  TFile* file=NULL;
  TTree* tree=NULL;

  ////////////////////////////////
  ////  Retrieving File Info  ////
  ////////////////////////////////

  ppFunct::getScanRunInfo(imgINFO, runListName, params.verbose);


  //////////////////////////////////////////////////////
  /////  Creating smaller runList files if needed  /////
  //////////////////////////////////////////////////////

  if (doRunLists) {
    ppFunct::makeRunLists(imgINFO);
  }

 
  //////////////////////////////////
  /////  Rough center finding  /////
  //////////////////////////////////

  std::vector<int> centersCOM = imgProc::centerSearchCOM(imgINFO,
      params.hotPixel, params.sigma,
      params.blockCentR, params.blockCentC, params.minRad, params.maxRad,
      params.meanInd, params.stdScale, params.verbose);

  if (params.verbose) std::cout << "COMcenters: " 
    << centersCOM[0] << "  " << centersCOM[1] << endl;

  ////////////////////////////////////////////////////////////
  /////  Create root file for each run-scan combination  /////
  ////////////////////////////////////////////////////////////

  curRun=""; curScan=-1;
  file=NULL; tree=NULL;
  for (ifl=0; ifl<imgINFO.size(); ifl++) {

    
    ////////////////////////////////////////
    // For a new run/scan do setup        //
    //    -Run/scan numbers               //
    //    -Number of images in scan       //
    //    -Close previous root file       //
    //    -Create new root file           //
    //    -Find image center              //
    ////////////////////////////////////////
    
    if ((curScan != imgINFO[ifl].scan) || (curRun != imgINFO[ifl].run)) {

      ///////  Get information about new run-scan combination  ///////
      runScanStartItr = ifl;
      curScan = imgINFO[ifl].scan;
      curRun = imgINFO[ifl].run;

      ///////  Create new output ROOT file  ///////
      if (file && tree) {
     	tree->Write();
  	file->Close();
      }

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

      cout << "  Making file " << rFileName << endl;

      file = TFile::Open(rFileName.c_str(), "RECREATE");
      tree = new TTree("physics","physics");

      tree->Branch("run", 	&curRun);
      tree->Branch("scan", 	&curScan,       "scan/I");
      tree->Branch("imgNum", 	&imgNum, 	"imgNum/I");
      tree->Branch("imgIsRef", 	&imgIsRef, 	"imgIsRef/I");
      tree->Branch("stagePos", 	&stagePos, 	"stagePos/I");
      tree->Branch("t0StagePos",&t0SP,		"t0StagePos/F");
      tree->Branch("t0Time",	&t0Time,	"t0Time/F");
      tree->Branch("centerC", 	&centerC, 	"centerC/I");
      tree->Branch("centerR", 	&centerR, 	"centerR/I");
      tree->Branch("imgNorm", 	&imgNorm, 	"imgNorm/F");
      tree->Branch("imgOrig", 	&imgOrig);
      tree->Branch("imgBkg", 	&imgBkg);
      tree->Branch("imgSubBkg", &imgSubBkg);
      tree->Branch("legCoeffs", &legCoeffs);
      if (params.verbose) cout << "INFO: Tree and file are setup!\n\n";


      ///// Finding image center /////
      if (params.verbose) cout << "INFO: Starting center Finding!\n\n";
      centerR = centerC = Ncents = 0;
      for (int icnt=ifl; icnt<std::min(ifl+params.NavgCenters, (uint)imgINFO.size()); icnt++) {

        /// Filling image std::vector ///
        string imgAddr = imgINFO[icnt].path + imgINFO[icnt].fileName;
        cv::Mat imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
        cout<<"image addr: "<<imgAddr<<endl;
        imgProc::threshold(imgMat, params.hotPixel);

        Nrows = imgMat.rows;
        Ncols = imgMat.cols;
        imgCent.resize(Nrows);
        int index;
        if (params.imgMatType.compare("uint16") == 0) {
          imgCent = imgProc::getImgVector<uint16_t>(imgMat, 
                      Nrows, Nrows/2, Ncols/2, 
                      params.holeRad, params.holeR, 
                      params.holeC, params.nanVal, false);
        }  
        else if (params.imgMatType.compare("uint32") == 0) {
          imgCent = imgProc::getImgVector<uint16_t>(imgMat, 
                      Nrows, Nrows, Ncols, 
                      params.holeRad, params.holeR, 
                      params.holeC, params.nanVal, false);
        }  
        else {
          cerr << "ERROR: Do not recognize imgMatType = " 
              + params.imgMatType + "!!!\n";
          exit(0);
        }
          


        // Remove readout noise
        if (params.verbose) std::cout << "\tSubtract readout noise.\n";
        imgProc::removeReadOutNoise(imgCent);

        // Center finding: fine tuned
        if (params.pltVerbose) {
          plt.printRC(imgCent, "center_holeRemoval" + to_string(imgINFO[icnt].imgNum), maximum, "200");
        }
        center[0] = centersCOM[0];
        center[1] = centersCOM[1];
        centfnctr.symVleg = params.symVLeg;
        centfnctr.nanVal = params.nanVal;
        centfnctr.minRadBin = params.minRadBin;
        centfnctr.centShellW = params.centShellW;
        centfnctr.img = &imgCent;

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
        if (params.verbose) cout << "\t" << center[0] << " " << center[1] << endl;

	// Show image center
        if (params.pltCent) {
          TH2F* cImg = plt.plotRC(imgCent, "centeredImage"+to_string(icnt), maximum, "1500");
          for (int ir=(center[0])-5; ir<=center[0]+5; ir++) {
            for (int ic=(center[1])-5; ic<=center[1]+5; ic++) {
              if (sqrt(pow(ir-center[0],2)+pow(ic-center[1],2)) < 5) cImg->SetBinContent(ic, ir, 0);
            }
          }
          plt.print2d(cImg, "centeredImage"+to_string(icnt));
          delete cImg;
        }

      }
      centerR /= Ncents;
      centerC /= Ncents;
      if (params.verbose) cout << "INFO: Found center " << centerR 
			<< " "<< centerC << "!\n\n";


      /////////////////////////////////////////////////////
      /////  Initializing reference images variables  /////
      /////////////////////////////////////////////////////

      /*
      ///  Finding reference images and summing them  ///
      cout << "    Making reference image!\n";
      if (hasRef) {
        string refAddr;
        cv::Mat refMat;
        std::vector< std::vector<double> > imgVec;
        for (uint ir=0; ir<imgOrig.size(); ir++) {
          std::fill(imgOrig[ir].begin(), imgOrig[ir].end(), 0.0);
        }

        for (uint i=0; i<Nrefs; i++) {
          // Retrieve image
          refMat = cv::imread(imgINFO[ifl+i].path + imgINFO[ifl+i].fileName,
              CV_LOAD_IMAGE_ANYDEPTH); 
          imgProc::threshold(refMat, params.hotPixel);

          if (params.imgMatType.compare("uint16") == 0) {
            imgVec = imgProc::getImgVector<uint16_t>(refMat, 
                imgSize, centerR, centerC, false);
          }
          else if (params.imgMatType.compare("uint32") == 0) {
            imgVec = imgProc::getImgVector<uint32_t>(refMat, 
                imgSize, centerR, centerC, false);
          }
          else {
            cerr << "ERROR: Do not recognize imgMatType = " 
                + params.imgMatType + "!!!\n";
            exit(0);
          }

          // Add images
          for (uint ir=0; ir<imgOrig.size(); ir++) {
            std::transform(imgVec[ir].begin(), imgVec[ir].end(), 
                imgOrig[ir].begin(), imgOrig[ir].begin(),
                std::plus<double>());
          }
        }
        
        for (uint ir=0; ir<imgOrig.size(); ir++) {
          tools::vScale<double>(imgOrig[ir], 1./((double)Nrefs));
        }
      }

      stagePos = 0;		
      imgIsRef = 1;
      imgNum = -1;

      /// Calculate image norm ///
      imgNorm = imgProc::imageNorm(imgOrig);

      /// Finish and save reference image ///
      tree->Fill();
      */
    } // Initial things to do for new file






    ///////////////////////////////////////
    /////  Processing images in scan  /////
    ///////////////////////////////////////

    // Skipping reference images, already processed above
    //if (hasRef && Nrefs) {
    //  Nrefs--;
    //  continue;
    //}
    
    // Looping through non reference images
    cout << "\tStage position: " 
      + to_string(imgINFO[ifl].stagePos) << endl;

    /////  Reference image done at beginning  /////
    imgIsRef = 0;
    if (params.hasRef) {
      if (imgINFO[ifl].stagePos < params.refStagePos) {
        imgIsRef = 1;
      }
    }

    /////  Image number (ordered)  /////
    imgNum = imgINFO[ifl].imgNum;

    /////  Filling various variable  /////
    stagePos = imgINFO[ifl].stagePos;
 
    ///////  Load image  ///////
    string imgAddr = imgINFO[ifl].path + imgINFO[ifl].fileName;
    if (params.verbose) cout << "INFO: Trying to open " << imgAddr << "\t .....";
    cv::Mat imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
    imgProc::threshold(imgMat, params.hotPixel);
    if (params.verbose) cout << "\tpassed!\n\n";

    if (params.imgMatType.compare("uint16") == 0) {
      imgOrig = imgProc::getImgVector<uint16_t>(imgMat, 
          imgSize, centerR, centerC, 
          params.holeRad, params.holeR, 
          params.holeC, params.nanVal, false);
    }
    else if (params.imgMatType.compare("uint32") == 0) {
      imgOrig = imgProc::getImgVector<uint32_t>(imgMat, 
          imgSize, centerR, centerC, 
          params.holeRad, params.holeR, 
          params.holeC, params.nanVal, false);
    }
    else {
      cerr << "ERROR: Do not recognize imgMatType = " 
          + params.imgMatType + "!!!\n";
      exit(0);
    }


    ///// Remove readout noise  /////
    if (params.verbose) std::cout << "INFO: Readout noise subtraction.\n";
    if (params.hasLaserBkg) {
      imgProc::removeAvgReadOutNoise(imgOrig, imgSize/2, imgSize/2, 
                  0.9*imgSize/2., 0.99*imgSize/2., params.nanVal, &params.nanMap);
    }
    else {
      imgProc::removeAvgReadOutNoise(imgOrig, imgSize/2, imgSize/2, 
                  0.9*imgSize/2., 0.99*imgSize/2., params.nanVal);
    }

    for (int ir=0; ir<imgSize; ir++) {
      for (int ic=0; ic<imgSize; ic++) {
        imgSubBkg[ir][ic] = imgOrig[ir][ic];
      }
    }

    //////////////////////////////////////////////////////
    /////  Finding and subtracting laser background  /////
    //////////////////////////////////////////////////////
    // We cluster background spots based on the assymetry of 
    //    the image. 

    if (params.verbose) std::cout << "INFO: Laser background removal.\n";
    if (params.hasLaserBkg) {
      // Fill in hole based on symmetry
      std::vector< std::vector<double> > imgLaser = imgOrig;
      int hCentR = imgSize/2;
      int hCentC = imgSize/2;
      /*
      for (int ir=hCentR-params.holeRad-3; 
          ir<=hCentR+params.holeRad+3; ir++) {
        for (int ic=hCentC-params.holeRad-3; 
            ic<=hCentC+params.holeRad+3; ic++) {
          if (std::pow(ir-hCentR,2) + std::pow(ic-hCentC,2) 
              < std::pow(params.holeRad, 2)) {
            imgOrig[ir][ic] = imgOrig[(imgSize-1)-ir][(imgSize-1)-ic];
          }
        }
      }
      */

      for (int ir=hCentR-params.holeRad-3; 
          ir<=hCentR+params.holeRad+3; ir++) {
        for (int ic=hCentC-params.holeRad-3; 
            ic<=hCentC+params.holeRad+3; ic++) {
          if (std::pow(ir-hCentR,2) + std::pow(ic-hCentC,2) 
              < std::pow(params.holeRad, 2)) {
            imgLaser[ir][ic] = imgLaser[(imgSize-1)-ir][(imgSize-1)-ic];
          }
        }
      }


      ///  Finding asymmetric parts of the image  ///
      std::vector< std::vector<double> > oddImgReal = imgProc::asymmetrize(imgLaser, 
          imgSize/2, imgSize/2, imgSize, imgSize, 
          oddImgImgn, fftFref, fftIn, fftBref, fftOut);

      // Building map of ratio of noise/"signal" = asymmetric/symmetric
      for (int ir=0; ir<imgSize; ir++) { 
        for (int ic=0; ic<imgSize; ic++) {
          symImg[ir][ic] = imgLaser[ir][ic] - oddImgReal[ir][ic];
          stdRatioImg[ir][ic] = oddImgReal[ir][ic]/std::max(std::abs(symImg[ir][ic]),0.01)
                                  /sqrt(pow(ir - imgSize/2, 2) + pow(ic - imgSize/2, 2));
        }
      }

      ///  Find large clusters corresponding to laser spots  ///
      std::vector< std::pair<uint, uint> > removePairs = imgProc::findClusters(
              stdRatioImg, params.holeRad*1.2,
              params.coreValThresh, params.coreRad, 
              params.minClusterSize, params.minPixelSize, params.minDensity, params.clusterRad,
              params.borderValThresh, params.borderRad, params.padRad, imgBkg);


      ///  Remove laser background from image by setting pixels = params.nanVal  ///
      for (uint ip=0; ip<removePairs.size(); ip++) {
        imgSubBkg[removePairs[ip].first][removePairs[ip].second] = params.nanVal;
      }
      /*
      ///  Remove laser background from image by setting pixels = params.nanVal  ///
      for (ir=0; ir<imgSize; ir++) {
        for (ic=0; ic<imgSize; ic++) {
          imgSubBkg[ir][ic] = imgOrig[ir][ic];
        }
      }

      for (ir=imgSize/2-40; ir<imgSize/2+40; ir++) {
        for (ic=imgSize/2-120; ic<imgSize/2; ic++) {
          imgSubBkg[ir][ic] = params.nanVal;
        }
      }
      for (ir=imgSize/2-200; ir<imgSize/2+50; ir++) {
        for (ic=imgSize-225; ic<imgSize-115; ic++) {
          imgSubBkg[ir][ic] = params.nanVal;
        }
      }
      for (ir=imgSize/2-60; ir<imgSize/2+35; ir++) {
        for (ic=imgSize-75; ic<imgSize; ic++) {
          imgSubBkg[ir][ic] = params.nanVal;
        }
      }
      */
      
      if (params.pltVerbose) {
        
        std::vector< std::vector<double> > pltStdRat(imgSize);
        for (int ir=0; ir<imgSize; ir++) {
          pltStdRat[ir].resize(imgSize);
          for (int ic=0; ic<imgSize; ic++) {
            if (std::pow(ir-hCentR,2) + std::pow(ic-hCentC,2) 
                < std::pow(params.holeRad, 2)) {
              pltStdRat[ir][ic] = 0;
            }
            else {
              pltStdRat[ir][ic] = stdRatioImg[ir][ic];
            }
          }
        }
 
        delete plt.printRC(oddImgReal, 
            "oddImgReal_"+curRun+"_"+to_string(curScan)+"_"+to_string(stagePos));//, 
            //pltOpts, pltVals);

        pltOpts[0] = minimum;	pltVals[0] = "0";
        pltOpts[1] = maximum;	pltVals[1] = "5e6";
        delete plt.printRC(imgBkg, 
            "imgBkg_"+curRun+"_"+to_string(curScan)+"_"+to_string(stagePos));

        pltOpts[0] = minimum;	pltVals[0] = "0"; //"1e4";
        pltOpts[1] = maximum;	pltVals[1] = to_string(params.coreValThresh); //"1e6";
        //pltOpts.push_back(logz);  pltVals.push_back("");
        delete plt.printRC(pltStdRat, 
            "imgStdRat_"+curRun+"_"+to_string(curScan)+"_"+to_string(stagePos),
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
            "imgStdRatRemove_"+curRun+"_"+to_string(curScan)+"_"+to_string(stagePos),
            pltOpts, pltVals);
      }
    }


    ////////////////////////////////////////
    /////  Plot results for debugging  /////
    ////////////////////////////////////////

    if (params.pltVerbose) {
      int hCentR = imgSize/2;
      int hCentC = imgSize/2;
      std::vector< std::vector<double> > pltOrig(imgSize);
      std::vector< std::vector<double> > pltSubBkg(imgSize);
      for (int ir=0; ir<imgSize; ir++) {
        pltOrig[ir].resize(imgSize);
        pltSubBkg[ir].resize(imgSize);
        for (int ic=0; ic<imgSize; ic++) {
          if (std::pow(ir-hCentR,2) + std::pow(ic-hCentC,2) 
              < std::pow(params.holeRad, 2)) {
            pltOrig[ir][ic] = 0;
            pltSubBkg[ir][ic] = 0;
          }
          else {
            pltOrig[ir][ic] = imgOrig[ir][ic];
            pltSubBkg[ir][ic] = imgSubBkg[ir][ic];
          }
        }
      }
      delete plt.printRC(pltOrig, 
          "imgOrig_"+curRun+"_"+to_string(curScan)+"_"+to_string(stagePos));
      delete plt.printRC(pltSubBkg, 
          "imgSub_"+curRun+"_"+to_string(curScan)+"_"+to_string(stagePos));
    }


    //////////////////////////////////////////////////////////
    /////  Filling image variables and filling the tree  /////
    //////////////////////////////////////////////////////////

    /////  Background subtraction and image norm  /////
    if (params.verbose) cout << "INFO: Subtracting background and calculating norm  .....  ";

    imgNorm = imgProc::imageNorm(imgOrig);

    if (params.verbose) cout << "passed!\n\n";

    //////  Legendre Fit  //////
    assert(imgSize%5 == 0);
    assert(imgSize%5 == 0);
    const int lgFit_Rows = imgSize/5;
    const int lgFit_Cols = imgSize/5;
    // Check if g matrix already exists, else make new one
    string matrix_folder = "/reg/neh/home/khegazy/analysis/legendreFitMatrices/";
    string matrix_fileName = "gMatrix_row-" + to_string(lgFit_Rows)
        + "_col-" + to_string(lgFit_Cols) + "_Nrad-" + to_string(NradBins)
        + "_Nlg-" + to_string(Nlegendres) + ".dat";
    if (params.verbose) std::cout << "INFO: Looking for " + matrix_fileName << endl;
    if (access((matrix_folder + matrix_fileName).c_str(), F_OK) == -1) {
      cout << "INFO: Making new g matrix\n";
      system(("python " + matrix_folder + "makeLgMatrix.py --NradBins="
            + to_string(NradBins) + " --Ncols=" + to_string(lgFit_Cols)
            + " --Nrows=" + to_string(lgFit_Rows)
            + " --Nlg=" + to_string(Nlegendres)).c_str());
    }

    // Import the g matrix
    const int NgMat = lgFit_Rows*lgFit_Cols*NradBins*Nlegendres;
    const int Npix = lgFit_Rows*lgFit_Cols;
    double* gInp = new double[NgMat];
    FILE* inpFile = fopen((matrix_folder + matrix_fileName).c_str(), "rb");
    fread(gInp, sizeof(double), NgMat, inpFile);
    Eigen::Map< Eigen::Matrix<double, Npix, 
        Nlegendres*NradBins, Eigen::RowMajor> > gOrig(gInp); 

    clock_t begin = clock();
    //legCoeffs = imgProc::legendreFit(imgSubBkg, 5, Nlegendres, 
    //    NradBins, lgFit_Rows, lgFit_Cols, params.nanVal, gOrig);
    
    std::fill(legCoeffs.begin(), legCoeffs.end(), 0);
    std::fill(legCounts.begin(), legCounts.end(), 0);
    for (int ir=0; ir<imgSubBkg.size(); ir++) {
      for (int ic=0; ic<imgSubBkg[ir].size(); ic++) {
        if (imgSubBkg[ir][ic] != params.nanVal) {
          double rad = sqrt(pow(ir-params.imgSize/2,2) + pow(ic-params.imgSize/2,2));
          int radInd = (int)(params.NradBins*rad/(params.imgSize/2));
          if (radInd >= params.NradBins) continue;
          legCoeffs[radInd] += imgSubBkg[ir][ic];
          legCounts[radInd] += 1;
        }
        else {
          //cout<<"skipping nan"<<endl;
        }
      }
    }

    for (int ir=0; ir<params.NradBins; ir++) {
      if (legCounts[ir] != 0) {
        legCoeffs[ir] /= legCounts[ir];
      }
      else legCoeffs[ir] = 0;
    }


    clock_t end = clock();
    cout<<"TIME: "<<double(end - begin) / CLOCKS_PER_SEC<<endl;

    if (params.pltVerbose) {
      std::vector<double> test(NradBins);
      for (int i=0; i<params.Nlegendres; i++) {
        for (int j=0; j<NradBins; j++) {
          test[j] = legCoeffs[i*NradBins + j];
        }
        plt.print1d(test, "testLeg" + to_string(i) + "_" + to_string(imgINFO[ifl].imgNum));
      }
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

                                        
