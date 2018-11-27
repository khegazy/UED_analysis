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
    iPos = fileList.find("RUN-");
  }
  iPos += 4;
  auto fPos = fileList.find(".txt");
  std::string runName = fileList.substr(iPos, 13); //fPos - iPos);
  cout<<"RunName: "<<runName<<endl;
  mergeClass merge(runName);

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
    else {
      cerr<<"ERROR!!! Option "<<argv[iarg]<<" does not exist!"<<endl;
      exit(0);
    }
  }

  /////  Plotting variables  /////
  std::vector<PLOToptions> opts(5);
  std::vector<std::string> vals(5);
  std::vector<PLOToptions> oppts(3);
  std::vector<std::string> vaals(3);
  oppts[0] = yLabel;   vaals[0] = "Time [ps]";
  oppts[1] = xLabel;   vaals[1] = "Scattering Q [arb units]";
  oppts[2] = draw;     vaals[2] = "COLZ";
  opts[0] = yLabel;   vals[0] = "Time [ps]";
  opts[1] = xLabel;   vals[1] = "Scattering Q [arb units]";
  opts[2] = draw;     vals[2] = "COLZ";
  opts[3] = minimum;  vals[3] = "-0.2";
  opts[4] = maximum;  vals[4] = "0.2";

  std::map<string, double > cbar;
  cbar["20161102_LongScan1_Leg0"] = 0.2;
  cbar["20161102_LongScan1_Leg2"] = 0.2;

  std::map< string, std::vector<double> > refAutCor, refLineOut;

  int curScan;
  
  ////////////////////////////////////////
  /////  Compare simulation results  /////
  ////////////////////////////////////////
  
  //if (merge.compareSims) {
  //  merge.compareSimulations(merge.molName);
  //}

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

    if (std::find(merge.badScans.begin(), merge.badScans.end(), scan) 
          != merge.badScans.end()) continue;
    auto badImgItr = merge.badImages.find(scan);
    if (badImgItr != merge.badImages.end()) {
      if (std::find(badImgItr->second.begin(), badImgItr->second.end(), stagePos) 
            != badImgItr->second.end()) {
        continue;
      }
    }

    merge.addLabTimeParameter(timeStamp, scan, stagePos, imgNorm);

    // Ignore reference images taken before scan
    if (imgIsRef) {
      if (merge.verbose) std::cout << "INFO: Adding reference image.\n";

      merge.addReference(scan, stagePos, filtAzmAvg, legCoeffs, 1);
      //merge.addReference(scan, stagePos, azmAvg, legCoeffs, imgNorm);
      continue;
    }

    ///  Insert entries to respective maps  ///
   
    merge.addEntry(scan, stagePos, filtAzmAvg, legCoeffs, 1);
    //merge.addEntry(scan, stagePos, azmAvg, legCoeffs, imgNorm);

  }


  ////////////////////////////////////////////
  /////  Subtract low order polynomials  /////
  ////////////////////////////////////////////
 
  if (merge.verbose) {
    std::cout << "INFO: Removing low order polynomials.\n";
  }
  merge.removeLowPolynomials();


  ///////////////////////////////////////////////////////////
  /////  Prune data for outliers and sparse time steps  /////
  ///////////////////////////////////////////////////////////


  //cout<<"removing outliers"<<endl;
  merge.removeOutliers();

  merge.mergeScans();

  // Get Mean and STD
  merge.getMeanSTD();


  /////  Saving  /////

  // Clean up NANVAL for saving and plotting
  for (int ir=0; ir<merge.azmReference.size(); ir++) {
    if (merge.azmReference[ir] == NANVAL) merge.azmReference[ir] = 0;
  }
  for (int ilg=0; ilg<merge.legReference.size(); ilg++) {
    for (int ir=0; ir<merge.legReference[ilg].size(); ir++) {
      if (merge.legReference[ilg][ir] == NANVAL) merge.legReference[ilg][ir] = 0;
    }
  }

  if (merge.verbose) 
    std::cout << "INFO: saving merged references and before t0 subtraction.\n";
  // References
  plt.print1d(merge.azmReference, "testingRef");
  save::saveDat<double>(merge.azmReference,
      merge.mergeScansOutputDir + 
      "/data-" + runName + "-" + prefix+ 
      "referenceAzm[" +
      to_string(merge.NradAzmBins) + "].dat");
  save::saveDat<double>(merge.legReference,
      merge.mergeScansOutputDir + 
      "/data-" + runName + "-" + prefix+ 
      "referenceLeg[" + 
      to_string(merge.Nlegendres) + 
      "," + to_string(merge.NradLegBins) + "].dat");
  save::saveDat<double>(merge.runAzmRefMeans,
      merge.mergeScansOutputDir + 
      "/data-" + runName + "-" + prefix + 
      "referenceAzmMean[" +
      to_string(merge.NradAzmBins) + "].dat");
  save::saveDat<double>(merge.runAzmRefSTD,
      merge.mergeScansOutputDir + 
      "/data-" + runName + "-" + prefix +
      "referenceAzmStandardDev[" +
      to_string(merge.NradAzmBins) + "].dat");
  save::saveDat<double>(merge.runAzmMeans,
      merge.mergeScansOutputDir + "data-"
      + runName + "-" + prefix + "mean[" 
      + to_string(merge.runAzmMeans.size()) + ","
      + to_string(merge.runAzmMeans[0].size()) + "].dat");
  save::saveDat<double>(merge.runAzmSTD,
      merge.mergeScansOutputDir + "data-"
      + runName + "-" + prefix + "standardDev[" 
      + to_string(merge.runAzmSTD.size()) + ","
      + to_string(merge.runAzmSTD[0].size()) + "].dat");

 

  merge.subtractT0();
  merge.normalize();

  // Get Mean and STD
  merge.getMeanSTD();

  //merge.smearTime();


  // Clean up NANVAL for saving and plotting
  for (int itm=0; itm<merge.azimuthalAvg.size(); itm++) {
    for (int ir=0; ir<merge.azimuthalAvg[itm].size(); ir++) {
      if (merge.azimuthalAvg[itm][ir] == NANVAL) merge.azimuthalAvg[itm][ir] = 0;
    }
  }


  /////  Saving and plotting  /////
  if (merge.verbose) 
    std::cout << "INFO: saving merged references and before t0 subtraction.\n";

  save::saveDat<double>(merge.azimuthalAvg, 
      merge.mergeScansOutputDir + "data-"
      + runName + "-" + prefix + "azmAvgDiff["
      + to_string(merge.azimuthalAvg.size()) + ","
      + to_string(merge.azimuthalAvg[0].size()) + "].dat");
  plt.printRC(merge.azimuthalAvg, 
      "./plots/data-"
      + runName + "-" + prefix + "azmAvgDiff", opts, vals);

  save::saveDat<double>(merge.azimuthalsMs, 
      merge.mergeScansOutputDir + "data-"
      + runName + "-" + prefix + "sMsAzmAvgDiff["
      + to_string(merge.azimuthalAvg.size()) + ","
      + to_string(merge.azimuthalAvg[0].size()) + "].dat");
  plt.printRC(merge.azimuthalAvg, 
      "./plots/data-"
      + runName + "-" + prefix + "sMsAzmAvgDiff", opts, vals);

  save::saveDat<double>(merge.runAzmMeans,
      merge.mergeScansOutputDir + "data-"
      + runName + "-" + prefix + "sMsMean[" 
      + to_string(merge.runAzmMeans.size()) + ","
      + to_string(merge.runAzmMeans[0].size()) + "].dat");
  save::saveDat<double>(merge.runAzmSTD,
      merge.mergeScansOutputDir + "data-"
      + runName + "-" + prefix + "sMsStandardDev[" 
      + to_string(merge.runAzmSTD.size()) + ","
      + to_string(merge.runAzmSTD[0].size()) + "].dat");

  save::saveDat<double>(merge.runAzmRefMeans,
      merge.mergeScansOutputDir + "data-"
      + runName + "-" + prefix
      + "referenceAzmsMsMean[" 
      + to_string(merge.NradAzmBins) + "].dat");
  save::saveDat<double>(merge.runAzmRefSTD,
      merge.mergeScansOutputDir + "data-"
      + runName + "-" + prefix
      + "referenceAzmsMsStandardDev[" 
      + to_string(merge.NradAzmBins) + "].dat");
 


  if (scanSearch) {
    cout<<"plotting to "<<merge.mergeScansOutputDir + "/plots/data-"+ runName + "-" + prefix + "azmAvgDiff";
    plt.printRC(merge.azimuthalAvg, 
        merge.mergeScansOutputDir + "/plots/data-"
        + runName + "-" + prefix + "azmAvgDiff", opts, vals);
  }


  /////////////////////////
  /////  Cleaning up  /////
  /////////////////////////

  if (merge.verbose) {
    std::cout << "Cleaning up" << endl;
  }

  return 1;
}
