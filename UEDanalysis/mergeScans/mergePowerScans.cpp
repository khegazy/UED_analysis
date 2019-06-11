#include "../analysis.h"
#include "mergeTools/mergeClass.h"
#include <TLegend.h>

using namespace std;


int main(int argc, char* argv[]) {

  if (argc<2) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe 'fileList.txt' !!!"<<endl;
    cerr<<"         Can also run ./analysis.exe 'fileList.txt' 'treeName'"<<endl;
    exit(0);
  }


  ///////////////////////////////////////////////////////////
  /////  Load environment and get the number of events  /////
  ///////////////////////////////////////////////////////////
  
  string fileList(argv[1]);
  string treeName("physics");
  if (argc==3) string treeName(argv[2]);
  analysisClass analysis(fileList, treeName);  // Use this when specifying treeName

  uint64_t Nentries;
  Nentries = analysis.setupEnvironment();  // Alter this function as needed for specific setup

  cout.setf(ios::scientific);
  
  //////////////////////////////////
  /////  Setting up variables  /////
  //////////////////////////////////

  auto iPos = fileList.find("run-");
  if (iPos == string::npos) {
    cerr << "ERROR: Cannot find 'power' in runName!\n";
    exit(0);
  }
  iPos += 4;
  auto fPos = fileList.find(".txt");
  std::string runName = fileList.substr(iPos, 8); //fPos - iPos);
  cout<<"RunName: "<<runName<<endl;
  mergeClass merge(runName);
  runName += "_PowerScan";

  std::string prefix = "";
  bool scanSearch = false;
  /////  Importing variables from command line  /////
  for (int iarg=2; iarg<argc; iarg+=2) {
    if (strcmp(argv[iarg],"-Odir") == 0) {
      string str(argv[iarg+1]); 
      merge.mergeScansOutputDir = str;
    }
    else if (strcmp(argv[iarg],"-ScanSearch") == 0) {
      string str(argv[iarg+1]);
      prefix = str + "-";
      merge.mergeScansOutputDir = merge.scanSearchOutputDir;
      scanSearch = true;
    }
    else if (strcmp(argv[iarg], "-OdirSuffix") == 0) {
      string str(argv[iarg+1]);
      merge.mergeScansOutputDir += (str + "/");

      if (!tools::fileExists(merge.mergeScansOutputDir)) {
        system(("mkdir " + merge.mergeScansOutputDir).c_str());
      }
    }
    else if (strcmp(argv[iarg], "-CompareRef") == 0) {
      std::string str(argv[iarg+1]);
      merge.compareReference(str);
    }
    else {
      cerr<<"ERROR!!! Option "<<argv[iarg]<<" does not exist!"<<endl;
      exit(0);
    }
  }


  /////  Plotting variables  /////
  std::vector<PLOToptions> opts(6);
  std::vector<std::string> vals(6);
  std::vector<PLOToptions> oppts(3);
  std::vector<std::string> vaals(3);
  oppts[0] = yLabel;   vaals[0] = "ROI Ratio";
  oppts[1] = xLabel;   vaals[1] = "Throttle [arb units]";
  oppts[2] = draw;     vaals[2] = "COLZ";
  opts[1] = xLabel;   vals[1] = "Q";
  opts[2] = draw;     vals[2] = "COLZ";
  opts[3] = minimum;  vals[3] = "-0.2";
  opts[4] = maximum;  vals[4] = "0.05";
  opts[5] = xSpan;    vals[5] = "0,"+to_string(merge.maxQazm);

  std::map< string, std::vector<double> > refAutCor, refLineOut;

  int curScan;
  
  ///////////////////////////
  /////  Merging Scans  /////
  ///////////////////////////

  if (merge.verbose)  {
    std::cout << "Begin to merge runs" << endl;
  }

  std::vector<string> runInds;
  std::map< string, std::map< double, double* > >  diffP_arrays;
  std::map< string, int > runShifts;

  ///// Loop through events in the file and saving to maps  /////
  curScan = -1;
  int pInd;
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);


    // Adding image parameters
    merge.addLabTimeParameter(timeStamp, "scan", scan);
    merge.addLabTimeParameter(timeStamp, "imgNorm", imgNorm);
    merge.addLabTimeParameter(timeStamp, "centerRstdRatio", centerRstdRatio);
    merge.addLabTimeParameter(timeStamp, "centerCstdRatio", centerCstdRatio);
    for (uint i=0; i<imgRadSTD->size(); i++) {
      merge.addLabTimeParameter(
          timeStamp, 
          "imgRadSTD" + to_string(i), 
          (*imgRadSTD)[i]);
    }
    // Adding lab/setup parameters
    for (auto const & pvInd : merge.pvMap) {
      if (pvInd.first.compare("pressure") == 0) {
        merge.addLabTimeParameter(timeStamp, "pressure", pressure);
      }
      else if (pvInd.first.compare("UVcounts") == 0) {
        merge.addLabTimeParameter(timeStamp, "UVcounts", UVcounts);
      }
      else if (pvInd.first.compare("bunkerTemp") == 0) {
        merge.addLabTimeParameter(timeStamp, "bunkerTemp", bunkerTemp);
      }
      else if (pvInd.first.compare("highBayTemp") == 0) {
        merge.addLabTimeParameter(timeStamp, "highBayTemp", bunkerTemp);
      }
      else {
        std::cerr << "ERROR: Do not know how to add PV " 
          << pvInd.first << "!!!\n";
        exit(1);
      }
    }


    // Ignore reference images taken before scan
    if (imgIsRef) {
      if (merge.verbose) std::cout << "INFO: Adding reference image.\n";
      
    merge.addLabTimeParameter(timeStamp, "stagePos", stagePos);
    //cout<<"ref add: "<<timeStamp<<" "<<stagePos<<endl;
      merge.addReference(
          scan, stagePos, 
          timeStamp, azmAvg, 
          legCoeffs, imgNorm);
      continue;
    }

    ///  Insert entries to respective maps  ///
    merge.addLabTimeParameter(timeStamp, "stagePos", throttle*1e8+stagePos);
    //cout<<"td add: "<<timeStamp<<" "<<throttle*1e8+stagePos<<endl;
    merge.addEntry(
        scan, throttle*1e8+stagePos, 
        timeStamp, azmAvg, 
        legCoeffs, imgNorm);

  }

  ///  Cut on lab time parameters  ///
  
  if (merge.verbose) {
    std::cout << "INFO: Cutting on lab time parameters.\n";
  }

  /*
  for (auto itrt = merge.stagePosInds.begin();
      itrt != merge.stagePosInds.end(); itrt++) {
    cout<<"itrt: "<<itrt->first<<" "<<itrt->second<<endl;
  }
  for (auto & itr : merge.labTimeMap) {
    cout<<"itr: "<<itr.first<<" "<<itr.second.first<<" "<<itr.second.second<<endl;
  }
  */
  /*
  merge.basicGreaterThanCut("centerRstdRatio", 3);
  merge.basicGreaterThanCut("centerCstdRatio", 3);
  for (uint i=0; i<imgRadSTD->size(); i++) {
    merge.stdParamCut("imgRadSTD" + to_string(i), 3);
  }
  */

  ///  Subtract low order polynomials  ///
 
  if (merge.verbose) {
    std::cout << "INFO: Removing low order polynomials.\n";
  }

  merge.removeImgNormOutliers();
  merge.removeLowPolynomials();


  ///  Prune data for outliers and sparse time steps  ///

  if (merge.verbose) {
    std::cout << "INFO: Removing outliers.\n";
  }


  merge.scaleByFit();

  merge.removeOutliers();

  
  ///  Merging Scans  ///

  if (merge.verbose) {
    std::cout << "INFO: Merging Scans.\n";
  }

  merge.mergeScans();

  // Normalize line outs
  //merge.normalizeScansResults();

  // Get Mean, STD, and SEM
  merge.getRunMeanSTDSEM();


  /////  Saving  /////

  // Clean up NANVAL for saving and plotting
  for (int ir=0; ir<merge.azmReference.size(); ir++) {
    if (merge.azmReference[ir] == NANVAL) {
      merge.azmReference[ir] = 0;
    }
  }
  for (int ilg=0; ilg<merge.legReference.size(); ilg++) {
    for (int ir=0; ir<merge.legReference[ilg].size(); ir++) {
      if (merge.legReference[ilg][ir] == NANVAL) {
        merge.legReference[ilg][ir] = 0;
      }
    }
  }


  /////  Time Domain Changes  /////

  // Subtract T0
  if (merge.verbose) std::cout << "INFO: Subtracting T0.\n";
  merge.subtractT0();

  /////  Normalize and get statistics  /////

  /*
  // Normalize to get sM(s)
  if (merge.verbose) std::cout << "INFO: sMs Normalization.\n";
  merge.sMsNormalize();

  // Gaussian smooth
  if (merge.verbose) std::cout << "INFO: Gaussian Smearing Q.\n";
  //merge.gaussianFilterQ();

  // Make pair correlations
  if (merge.verbose) std::cout << "INFO: Making pair correlations.\n";
  merge.makePairCorrs();

  // Get Mean, STD, and SEM
  if (merge.verbose) std::cout << "INFO: Calculating mean, STD, SEM.\n";
  merge.getRunMeanSTDSEM();
  */


  // Calculating Differences
  double range1, range2;
  int imgNormBinMin = (int)(merge.NradAzmBins*3.5/merge.maxQazm);
  int imgNormBinMax = (int)(merge.NradAzmBins*4.9/merge.maxQazm);
  int range1Begin   = (int)(merge.NradAzmBins*merge.range1Qbegin/merge.maxQazm);
  int range1End     = (int)(merge.NradAzmBins*merge.range1Qend/merge.maxQazm);
  int range2Begin   = (int)(merge.NradAzmBins*merge.range2Qbegin/merge.maxQazm);
  int range2End     = (int)(merge.NradAzmBins*merge.range2Qend/merge.maxQazm);

  std::vector<double> throttles;
  std::vector<double> smooth;
  std::vector<double> curRatio(merge.NradAzmBins);
  std::vector<double> powerScanRatio, powerScanRatioSEM;

  for (auto itrt = merge.stagePosInds.begin();
      itrt != merge.stagePosInds.end(); itrt++) {

    smooth = imgProc::gaussianSmooth1d(merge.azimuthalAvg[itrt->second], 7, 35);
    for (int iq=0; iq<merge.NradAzmBins; iq++) {
      merge.azimuthalAvg[itrt->second][iq] = smooth[iq];
    }

    smooth = imgProc::gaussianSmooth1d(merge.azimuthalsMs[itrt->second], 7, 35);
    for (int iq=0; iq<merge.NradAzmBins; iq++) {
      merge.azimuthalsMs[itrt->second][iq] = smooth[iq];
    }

    for (int iq=0; iq<merge.NradAzmBins; iq++) {
      if ((merge.azimuthalAvg[itrt->second][iq] == 0) ||
          (merge.azimuthalAvg[itrt->second][iq] == NANVAL) ||
          (merge.azmReference[iq] == 0) ||
          (merge.azmReference[iq] == NANVAL)) {
        curRatio[iq] = NANVAL;
      }
      else{
        curRatio[iq] = merge.azimuthalAvg[itrt->second][iq]
                        /merge.azmReference[iq];
      }
    }
    
    // Summing Transient Signal
    double mean   = 0;
    double stdev  = 0;
    for (int iq=imgNormBinMin; iq<imgNormBinMax; iq++) {
      //powerScanRatio[throttleCount] += fabs(powerScanDiff[iq])*merge.sMsAzmNorm[iq];
      //mean += (merge.azimuthalAvg[itrt->second][iq])/merge.atmAzmDiff[iq];
      mean += curRatio[iq];
    }
    mean /= imgNormBinMax - imgNormBinMin;
    for (int iq=imgNormBinMin; iq<imgNormBinMax; iq++) {
      stdev += std::pow(mean 
                - curRatio[iq],
                //- (merge.azimuthalAvg[itrt->second][iq])/merge.atmAzmDiff[iq],
                2);
    }
    stdev = std::sqrt(stdev/(imgNormBinMax - imgNormBinMin))
              /std::sqrt(double(imgNormBinMax - imgNormBinMin));

    powerScanRatio.push_back(mean);
    powerScanRatioSEM.push_back(stdev);


    // Saving and plotting 
    throttle = (int)(itrt->first*1e-8);
    throttles.push_back(throttle);
    opts[3] = minimum;  vals[3] = "-0.002";
    opts[4] = maximum;  vals[4] = "0.002";
    prefix = "throttle-" + to_string(throttle) + "_";
    save::saveDat<double>(merge.azimuthalAvg[itrt->second],
        merge.mergeScansOutputDir + "data-"
        + runName + "-" + prefix + "azmAvgDiff["
        + to_string(merge.NradAzmBins) + "].dat");
    plt.print1d(merge.azimuthalAvg[itrt->second],
        "./plots/data-"
        + runName + "-" + prefix + "azmAvgDiff", opts, vals);


    save::saveDat<double>(merge.azimuthalsMs[itrt->second],
        merge.mergeScansOutputDir + "data-"
        + runName + "-" + prefix + "sMsAzmAvgDiff["
        + to_string(merge.NradAzmBins) + "].dat");
    opts[3] = minimum;  vals[3] = "-0.05";
    opts[4] = maximum;  vals[4] = "0.05";
    plt.print1d(merge.azimuthalsMs[itrt->second],
        "./plots/data-"
        + runName + "-" + prefix + "sMsAzmAvgDiff", opts, vals);

  }


  for (int v=0; v<powerScanRatio.size(); v++) {
    cout<<"vals: "<<v<<" "<<powerScanRatio[v]<<endl;
  }
  save::saveDat<double>(powerScanRatio,
      "./results/data-"
      + runName + "_powerScanLinearity["
      + to_string(powerScanRatio.size()) + "].dat");
  save::saveDat<double>(powerScanRatioSEM,
      "./results/data-"
      + runName + "_powerScanLinearitySEM["
      + to_string(powerScanRatioSEM.size()) + "].dat");
  save::saveDat<double>(throttles,
      "./results/data-"
      + runName + "_throttles["
      + to_string(powerScanRatio.size()) + "].dat");

  plt.print1d(powerScanRatio,
      "./plots/data-"
      + runName + "_powerScanLinearity", oppts, vaals);



  /////////////////////////
  /////  Cleaning up  /////
  /////////////////////////

  if (merge.verbose) {
    std::cout << "Cleaning up" << endl;
  }

  return 1;
}
