#include "preProcessing.h"
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
          //N2 return imgProc::centerOddLeg(img, vect[0], vect[1], 1, 70, 40, true, nanVal);
        }
};
         
/*
class centerfnctr {

    public:
        std::vector< std::vector<double> >* img;
        int count;
        int minRadBin, centShellW;
 	double nanVal;
        //centerfnctr(std::vector< std::vector<double> > &img_i) {img = &img_i;}
        //centerfnctr(); 
        double operator() (std::vector<double> vect) {
          count++;
          //cout<<"calling fxn          ";
          return imgProc::centerSymXsqr(img, vect[0], vect[1], 8, 1, minRadBin, centShellW, true, nanVal);
          //return fabs(1-imgProc::centerSymXsqr(img, vect[0], vect[1], 8, 3, 100, 20));
        }  
};
*/ 


int main(int argc, char* argv[]) {

  experimentParameters  expP = experimentParameters();
  std::string outputDir = "/reg/ued/ana/scratch/nitroBenzene/rootFilesNew/";

  // Image Parameters
  const int roiw = 835;
  const int roih = 835;
  const int hotPixel = 50000;

  // Number of Legendres
  const int Nlg = 6;
  const int NradBins = 50;

  // Laser Background Removal Parameters
  const double  decayConst = -1.0/700.0;
  const double  coreValThresh = 60; //65; //5e5;
  const int     coreRad = 4;
  const int     minClusterSize = 450;
  const int     minPixelSize = 200;
  const double  minDensity = 0;//0.15;//0.2;
  const int     clusterRad = 1;
  const double  borderValThresh = 26; //1e5;
  const int     borderRad = 5;
  const int     padRad = 10;

  // Center Finding Parameters
  int symVLeg = 1;
  int minRadBin = 60;
  int centShellW = 70;
  int holeR = 547;
  int holeC = 492;
  int holeRad = 60;
  int uv1R = 580;
  int uv1C = 415;
  int uv1Rad = 70;
  double nanVal = 1.23456789e-12;

  // Debugging
  bool pltCent = false;
  bool verbose = false;
  bool pltVerbose = false;

  // Make runLists
  bool doRunLists = false;


  if (argc != 2) {
    cerr << "ERROR: Must run program a list of files" 
      << " ./preProcessing.exe runList.txt!!!" << endl;
    exit(0);
  }

  string runListName(argv[1]);

  // Make sure the image has odd number of bins
  if (!(roih%2) || !(roiw%2)) {
    cerr << "ERROR: roih and roiw must be an odd number!!!" << endl;
    exit(0);
  }
 

  PLOTclass plt;

  // Indices
  int ir, ic;
  int irow, icol, index;
  uint k, ifl, iffl, runScanStartItr;
  std::vector<PLOToptions> pltOpts(2);
  std::vector<string> pltVals(2);


  std::vector<double> legCoeffs;

  std::vector< std::vector<double> > imgOrig(roih);
  std::vector< std::vector<double> > imgSubBkg(roih);
  std::vector< std::vector<double> > imgTemp(roih);
  std::vector< std::vector<double> > imgBkg(roih);
  std::vector< std::vector<double> > imgCent;
  float imgNorm;
  for (ir=0; ir<roih; ir++) {
    imgOrig[ir].resize(roiw, 0);
    imgSubBkg[ir].resize(roiw, 0);
    imgTemp[ir].resize(roiw, 0);
    imgBkg[ir].resize(roiw, 0);
  }

  std::vector< std::vector<double> > oddImgImgn;
  std::vector< std::vector<double> > oddImgReal;
  std::vector< std::vector<double> > symImg(roih);
  std::vector< std::vector<double> > stdRatioImg(roih);
  for (ir=0; ir<roih; ir++) { 
    symImg[ir].resize(roiw, 0);
    stdRatioImg[ir].resize(roiw, 0);
  }

  fftw_complex* fftIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*roih*roiw);
  fftw_complex* fftOut = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*roih*roiw);
  fftw_plan fftFref = fftw_plan_dft_2d(roih, roiw, fftIn, fftOut, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan fftBref = fftw_plan_dft_2d(roih, roiw, fftIn, fftOut, FFTW_BACKWARD, FFTW_MEASURE);


  std::vector< std::vector<double> > radBins;

  centerfnctr centfnctr;
  std::vector<double> center(2);
  int centerC, centerR, Ncents;

  int Nrows, Ncols;
  double count;

  int imgIsRef;
  int Nrefs;
  bool hasRef;

  int imgNum;
  int curRun;
  string fileName;
  string line;
  size_t ipos, fpos;
  string date, scan, curScan, curDate, rFileName;
  float stagePos;
  float t0SP, t0Time;
  std::vector<imgInfoStruct> imgINFO;

  TFile* file=NULL;
  TTree* tree=NULL;

  ////////////////////////////////
  ////  Retrieving File Info  ////
  ////////////////////////////////

  ppFunct::getScanRunInfo(imgINFO, runListName, verbose);


  //////////////////////////////////////////////////////
  /////  Creating smaller runList files if needed  /////
  //////////////////////////////////////////////////////

  if (doRunLists) {
    ppFunct::makeRunLists(imgINFO);
  }


  ////////////////////////////////////////////////////////////
  /////  Create root file for each run-scan combination  /////
  ////////////////////////////////////////////////////////////

  curScan=""; curDate=""; curRun=-1;
  file=NULL; tree=NULL;
  for (ifl=0; ifl<imgINFO.size(); ifl++) {

    
    ////////////////////////////////////////
    // For a new run/scan do setup        //
    //    -Run/scan numbers               //
    //    -Number of images in scan       //
    //    -Close previous root file       //
    //    -Create new root file           //
    //    -Find image center              //
    //    -Check for reference images     //
    //    -Make reference image           //
    //    -Save reference image           //
    ////////////////////////////////////////
    
    if ((curScan != imgINFO[ifl].scan) || (curRun != imgINFO[ifl].run) 
        || (curDate != imgINFO[ifl].date)) {

      ///////  Get information about new run-scan combination  ///////
      runScanStartItr = ifl;
      curScan = imgINFO[ifl].scan;
      curRun = imgINFO[ifl].run;
      curDate = imgINFO[ifl].date; 

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

      rFileName =  outputDir + "/" + subFolder + "/"
      //rFileName =  "testing/" + subFolder + "/"
            + "Date-" + curDate + "_" 
            + "Scan-" + curScan + "_"
            + "Run-" + to_string(curRun) + ".root";

      cout << "  Making file " << rFileName << endl;

      file = TFile::Open(rFileName.c_str(), "RECREATE");
      tree = new TTree("physics","physics");

      tree->Branch("date", 	&curDate);
      tree->Branch("scan", 	&curScan);
      tree->Branch("runNum", 	&curRun, 	"runNum/I");
      tree->Branch("imgNum", 	&imgNum, 	"imgNum/I");
      tree->Branch("imgIsRef", 	&imgIsRef, 	"imgIsRef/I");
      tree->Branch("stagePos", 	&stagePos, 	"stagePos/F");
      tree->Branch("t0StagePos",&t0SP,		"t0StagePos/F");
      tree->Branch("t0Time",	&t0Time,	"t0Time/F");
      tree->Branch("centerC", 	&centerC, 	"centerC/I");
      tree->Branch("centerR", 	&centerR, 	"centerR/I");
      tree->Branch("imgNorm", 	&imgNorm, 	"imgNorm/F");
      tree->Branch("nanVal", 	&nanVal,   "nanVal/F");
      tree->Branch("imgOrig", 	&imgOrig);
      tree->Branch("imgBkg", 	&imgBkg);
      tree->Branch("imgSubBkg", &imgSubBkg);
      tree->Branch("legCoeffs", &legCoeffs);
      if (verbose) cout << "INFO: Tree and file are setup!\n\n";


      ///// Finding image center /////
      if (verbose) cout << "INFO: Starting center Finding!\n\n";
      centerR = centerC = Ncents = 0;
      for (int icnt=ifl; icnt<std::min(ifl+3, (uint)imgINFO.size()); icnt++) {

        /// Filling image std::vector ///
        string imgAddr = imgINFO[icnt].path + imgINFO[icnt].fileName;
        cv::Mat imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
        cout<<"image addr: "<<imgAddr<<endl;
        imgProc::threshold(imgMat, hotPixel);

        bool hasLaserBkg = false;
        if (expP.laserBkg.find(imgINFO[ifl].date) != expP.laserBkg.end()) {
          if (expP.laserBkg[imgINFO[ifl].date].find(imgINFO[ifl].scan)
                 != expP.laserBkg[imgINFO[ifl].date].end()) {
            hasLaserBkg = true;
          }
        }


        Nrows = imgMat.rows;
        Ncols = imgMat.cols;
        imgCent.resize(Nrows);
        uint32_t* pv = imgMat.ptr<uint32_t>(0);
        for (ir=0; ir<Nrows; ir++) {
          imgCent[ir].resize(Ncols);
          for (ic=0; ic<Ncols; ic++) {
            index = ir*Ncols + ic;
            if (sqrt(pow(ir-holeR,2) + pow(ic-holeC,2)) < holeRad) {
              imgCent[ir][ic] = nanVal;
            }
            else if (sqrt(pow(ir-uv1R,2) + pow(ic-uv1C,2)) < uv1Rad) {
              imgCent[ir][ic] = nanVal;
            }
            else {
              imgCent[ir][ic] = pv[index];
            }
          }
        }

        // Remove readout noise
        imgProc::removeReadOutNoise(imgCent);

        center[0] = 575;
        center[1] = 495;
        centfnctr.symVleg = symVLeg;
        centfnctr.nanVal = nanVal;
        centfnctr.minRadBin = minRadBin;
        centfnctr.centShellW = centShellW;
        centfnctr.img = &imgCent;

        tools::powellMin<centerfnctr> (centfnctr, center, 10, 1, 1, 0.01);

        if ((curDate == "20131102") && (curScan == "LongScan1")) {
          center[0] = 569;
          center[1] = 490;
        }
        centerR += center[0];
        centerC += center[1];
	Ncents += 1;
        if (verbose) cout << "\t" << center[0] << " " << center[1] << endl;

	// Show image center
        if (pltCent) {
          TH2F* cImg = plt.plotRC(imgCent, "centeredImage"+to_string(icnt), maximum, "1500");
          for (ir=(center[0])-5; ir<=center[0]+5; ir++) {
            for (ic=(center[1])-5; ic<=center[1]+5; ic++) {
              if (sqrt(pow(ir-center[0],2)+pow(ic-center[1],2)) < 5) cImg->SetBinContent(ic, ir, 0);
            }
          }
          plt.print2d(cImg, "centeredImage"+to_string(icnt));
          delete cImg;
        }

      }
      centerR /= Ncents;
      centerC /= Ncents;
      if (verbose) cout << "INFO: Found center " << centerR 
			<< " "<< centerC << "!\n\n";


      ///////  Initializing reference images variables  ///////
      Nrefs = 0;
      hasRef = false;	
      // Does scan on day 'date' have a reference
      if (expP.refStagePos.find(imgINFO[ifl].date) != expP.refStagePos.end()) {
	if (expP.refStagePos[imgINFO[ifl].date].find(imgINFO[ifl].scan) 
		!= expP.refStagePos[imgINFO[ifl].date].end()) {
	  hasRef = true;
	}
      }

      // Find t0 from experiment log
      t0SP = t0Time = 0;
      if (expP.t0StagePos.find(imgINFO[ifl].date) != expP.t0StagePos.end()) {
	if (expP.t0StagePos[imgINFO[ifl].date].find(imgINFO[ifl].scan) 
		!= expP.t0StagePos[imgINFO[ifl].date].end()) {
	  t0SP = expP.t0StagePos[imgINFO[ifl].date][imgINFO[ifl].scan];
	  t0Time = (t0SP/C_SI)*1e15;
	}
      }


      ///  Finding reference images and summing them  ///
      cout << "    Making reference image!\n";
      string refAddr;
      cv::Mat refMat;
      imgNorm = 0;
      if (hasRef) {
        iffl = ifl;
        std::vector< std::vector<double> > inpImg;
        while ((curScan == imgINFO[iffl].scan) 
	    && (curRun == imgINFO[iffl].run) && (curDate == imgINFO[iffl].date)
            && (iffl < (int)imgINFO.size())) {
          if (imgINFO[iffl].stagePos < expP.refStagePos[imgINFO[iffl].date][imgINFO[iffl].scan]) {
            refAddr = imgINFO[iffl].path+imgINFO[iffl].fileName;
            refMat = cv::imread(refAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
            imgProc::threshold(refMat, hotPixel);

            inpImg = imgProc::getImgVector(refMat, roih, roiw, centerR, centerC);
            for (ir=0; ir<roih; ir++) {
              for (ic=0; ic<roiw; ic++) {
                imgOrig[ir][ic] += inpImg[ir][ic];
              }
            }

            Nrefs++;
          }
          else {
            break;
          }
          iffl++;
        }

        for (ir=0; ir<roih; ir++) {
          for (ic=0; ic<roiw; ic++) {
            imgOrig[ir][ic] /= Nrefs;
          }
        }
      }
      // If no reference take first image as reference
      else {
        refAddr = imgINFO[ifl].path + imgINFO[ifl].fileName;
        refMat = cv::imread(refAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
        imgProc::threshold(refMat, hotPixel);
        uint32_t* pv = refMat.ptr<uint32_t>(0);
        for (ir=0; ir<roih; ir++) {
          irow = ir - roih/2 + centerR;
          for (ic=0; ic<roiw; ic++) {
            icol = ic - roiw/2 + centerC;
            index = irow*Ncols + icol;
            imgOrig[ir][ic] = pv[index];
          }
        }
      }

      stagePos = 0;		
      imgIsRef = 1;
      imgNum = -1;

      /// Calculate image norm ///
      imgNorm = imgProc::imageNorm(imgOrig);

      /// Finish and save reference image ///
      tree->Fill();
    } // Initial things to do for new file






    ///////////////////////////////////////
    /////  Processing images in scan  /////
    ///////////////////////////////////////

    // Skipping reference images, already processed above
    if (hasRef && Nrefs) {
      Nrefs--;
      continue;
    }
    
    // Looping through non reference images
    cout << "\tStage position: " 
      + to_string(imgINFO[ifl].stagePos) << endl;

    /////  Reference image done at beginning  /////
    imgIsRef = 0;

    /////  Image number (ordered)  /////
    imgNum = k;

    /////  Filling various variable  /////
    stagePos = imgINFO[ifl].stagePos;
 
    ///////  Load image  ///////
    string imgAddr = imgINFO[ifl].path + imgINFO[ifl].fileName;
    if (verbose) cout << "INFO: Trying to open " << imgAddr << "\t .....";
    cv::Mat imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
    imgProc::threshold(imgMat, hotPixel);
    if (verbose) cout << "\tpassed!\n\n";

    imgOrig = imgProc::getImgVector(imgMat, roih, roiw, centerR, centerC);

    // Fill in hole based on symmetry
    int hCentR = roih/2;
    int hCentC = roiw/2;
    for (ir=hCentR-holeRad-3; ir<=hCentR+holeRad+3; ir++) {
      for (ic=hCentC-holeRad-3; ic<=hCentC+holeRad+3; ic++) {
        if (pow(ir-hCentR,2) + pow(ic-hCentC,2) < holeRad*holeRad) {
          imgOrig[ir][ic] = imgOrig[(roih-1)-ir][(roiw-1)-ic];
        }
      }
    }


    ///////  Finding and subtracting background  ///////

    ///  Finding asymmetric parts of the image  ///
    oddImgReal = imgProc::asymmetrize(imgOrig, roih/2, roiw/2, roih, roiw, 
        oddImgImgn, fftFref, fftIn, fftBref, fftOut);

    // Building map of ratio of noise/"signal" = asymmetric/symmetric
    for (ir=0; ir<roih; ir++) { 
      for (ic=0; ic<roiw; ic++) {
        symImg[ir][ic] = imgOrig[ir][ic] - oddImgReal[ir][ic];
        stdRatioImg[ir][ic] = std::exp(decayConst*sqrt(pow(ir - roih/2, 2) + pow(ic - roiw/2, 2)))
                                *oddImgReal[ir][ic]/std::max(std::abs(symImg[ir][ic]),0.01);
        //stdRatioImg[ir][ic] = //std::sqrt(pow(ir - roih/2, 2) + pow(ic - roiw/2, 2))
          //                          oddImgReal[ir][ic]
            //                        /std::sqrt(std::sqrt(pow(ir - roih/2, 2) + pow(ic - roiw/2, 2)))
              //                      /((std::max(symImg[ir][ic],0.0)) + 0.01); //sqrt(pow(ir - roih/2, 2) + pow(ic - roiw/2, 2));
        //cout<<"stdImg "<<ir<<" "<<ic<<"   "<<stdRatioImg[ir][ic]<<"  "<<(pow(ir - roih/2, 2) + pow(ic - roiw/2, 2))<<"  "<<oddImgReal[ir][ic]<<"  "<<(sqrt(std::abs(symImg[ir][ic])) + 0.5)<<"  "<<symImg[ir][ic]<<endl;
      }
    }

    ///  Find large clusters corresponding to laser spots  ///
    std::vector< std::pair<uint, uint> > removePairs = imgProc::findClusters(
            stdRatioImg, holeRad*1.2,
            coreValThresh, coreRad, 
            minClusterSize, minPixelSize, minDensity, clusterRad,
            borderValThresh, borderRad, padRad, imgBkg);


    ///  Remove laser background from image by setting pixels = nanVal  ///
    for (ir=0; ir<roih; ir++) {
      for (ic=0; ic<roiw; ic++) {
        imgSubBkg[ir][ic] = imgOrig[ir][ic];
      }
    }
    for (uint ip=0; ip<removePairs.size(); ip++) {
      imgSubBkg[removePairs[ip].first][removePairs[ip].second] = nanVal;
    }

    if (pltVerbose) {
      std::vector< std::vector<double> > pltOrig(roih);
      std::vector< std::vector<double> > pltSubBkg(roih);
      std::vector< std::vector<double> > pltStdRat(roih);
      for (ir=0; ir<roih; ir++) {
        pltOrig[ir].resize(roiw);
        pltSubBkg[ir].resize(roiw);
        pltStdRat[ir].resize(roiw);
        for (ic=0; ic<roiw; ic++) {
          if (pow(ir-hCentR,2) + pow(ic-hCentC,2) < holeRad*holeRad) {
            pltOrig[ir][ic] = 0;
            pltSubBkg[ir][ic] = 0;
            pltStdRat[ir][ic] = 0;
          }
          else {
            pltOrig[ir][ic] = imgOrig[ir][ic];
            pltSubBkg[ir][ic] = imgSubBkg[ir][ic];
            pltStdRat[ir][ic] = stdRatioImg[ir][ic];
          }
        }
      }
      delete plt.printRC(pltOrig, 
          "imgOrig_"+curScan+"_"+to_string(curRun)+"_"+to_string(stagePos));
      delete plt.printRC(oddImgReal, 
          "oddImgReal_"+curScan+"_"+to_string(curRun)+"_"+to_string(stagePos));//, 
          //pltOpts, pltVals);
      delete plt.printRC(pltSubBkg, 
          "imgSub_"+curScan+"_"+to_string(curRun)+"_"+to_string(stagePos));

      pltOpts[0] = minimum;	pltVals[0] = "0";
      pltOpts[1] = maximum;	pltVals[1] = "5e6";
      delete plt.printRC(imgBkg, 
          "imgBkg_"+curScan+"_"+to_string(curRun)+"_"+to_string(stagePos));

      pltOpts[0] = minimum;	pltVals[0] = "1e4";
      pltOpts[1] = maximum;	pltVals[1] = "1e6";
      //pltOpts.push_back(logz);  pltVals.push_back("");
      delete plt.printRC(pltStdRat, 
          "imgStdRat_"+curScan+"_"+to_string(curRun)+"_"+to_string(stagePos));//,
          //pltOpts, pltVals);
    }


    ////////////////////////////////////////////////////
    //  Filling image variables and filling the tree  //
    ////////////////////////////////////////////////////

    /////  Background subtraction and image norm  /////
    if (verbose) cout << "INFO: Subtracting background and calculating norm  .....  ";

    imgNorm = imgProc::imageNorm(imgOrig);

    if (verbose) cout << "passed!\n\n";

    //////  Legendre Fit  //////
    assert(roih%5 == 0);
    assert(roiw%5 == 0);
    const int lgFit_Rows = roih/5;
    const int lgFit_Cols = roiw/5;
    // Check if g matrix already exists, else make new one
    string matrix_folder = "/reg/neh/home/khegazy/analysis/legendreFitMatrices/";
    string matrix_fileName = "gMatrix_row-" + to_string(lgFit_Rows)
        + "_col-" + to_string(lgFit_Cols) + "_Nrad-" + to_string(NradBins)
        + "_Nlg-" + to_string(Nlg) + ".dat";
    if (access((matrix_folder + matrix_fileName).c_str(), F_OK) == -1) {
      cout << "INFO: Making new g matrix\n";
      system(("python " + matrix_folder + "makeLgMatrix.py --NradBins="
            + to_string(NradBins) + " --Ncols=" + to_string(lgFit_Cols)
            + " --Nrows=" + to_string(lgFit_Rows)
            + " --Nlg=" + to_string(Nlg)).c_str());
    }

    // Import the g matrix
    const int NgMat = lgFit_Rows*lgFit_Cols*NradBins*Nlg;
    const int Npix = lgFit_Rows*lgFit_Cols;
    double* gInp = new double[NgMat];
    FILE* inpFile = fopen((matrix_folder + matrix_fileName).c_str(), "rb");
    fread(gInp, sizeof(double), NgMat, inpFile);
    Eigen::Map< Eigen::Matrix<double, Npix, Nlg*NradBins, Eigen::RowMajor> > gOrig(gInp); 

    clock_t begin = clock();
    legCoeffs = imgProc::legendreFit(imgSubBkg, 5, Nlg, NradBins, lgFit_Rows, lgFit_Cols, nanVal, gOrig);
    clock_t end = clock();
    cout<<"TIME: "<<double(end - begin) / CLOCKS_PER_SEC<<endl;

    if (pltVerbose) {
      std::vector<double> test(NradBins);
      for (int i=0; i<Nlg; i++) {
        for (int j=0; j<NradBins; j++) {
          test[j] = legCoeffs[i*NradBins + j];
        }
        plt.print1d(test, "testLeg" + to_string(i));
      }
      cout<<"done plotting"<<endl;
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

                                        
