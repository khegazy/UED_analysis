#include "../analysis.h"
#include "mergeTools/mergeClass.h"
#include <TLegend.h>

using namespace std;


int main(int argc, char* argv[]) {

  if (argc<2) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe 'fileList.txt' !!!"<<endl;
    cerr<<"         Can also run ./analysis.exe 'fileList.txt' 'treeName'"<<endl;
    exit(1);
  }


  ///////////////////////////////////////////////////////////
  /////  Load environment and get the number of events  /////
  ///////////////////////////////////////////////////////////
 
  std::string input_run(argv[1]);
  analysisClass analysis(input_run);
  std::vector<std::string> variables{
      "imgNum", "stagePos", "timeStamp", "imgNorm", "I0norm", "imgIsRef",
      "centerR", "centerC", "centerRstdRatio", "centerCstdRatio",
      "I0centerR", "I0centerC", "imgRadSTD", "pressure", "UVcounts", 
      "legCoeffs", "legCoeffs_nanMap", "azmAvg", "azmAvg_nanMap"};
  cout<<"scan size "<<scans.size();
  analysis.init_preProcesses_data(variables);
  cout<<"  "<<imgNums->size()<<endl;

  uint64_t Nentries = scans.size();
  cout.setf(ios::scientific);
  
  //////////////////////////////////
  /////  Setting up variables  /////
  //////////////////////////////////

  cout<<"RunName: "<<analysis.run<<endl;
  mergeClass merge(analysis.run);

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
      exit(1);
    }
  }

  std::map< string, std::vector<double> > refAutCor, refLineOut;
  std::vector<long int> labTimeStamps;

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
  opts[3] = minimum;  vals[3] = "-0.1";
  opts[4] = maximum;  vals[4] = "0.1";

  std::map<string, double > cbar;
  cbar["20161102_LongScan1_Leg0"] = 0.2;
  cbar["20161102_LongScan1_Leg2"] = 0.2;


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
  std::vector<int> remove_inds;
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {

    if (std::find(merge.badScans.begin(), merge.badScans.end(), scans[ievt]) 
          != merge.badScans.end()) {
      remove_inds.push_back(ievt);
      continue;
    }
    auto badImgItr = merge.badImages.find(scans[ievt]);
    if (badImgItr != merge.badImages.end()) {
      if (std::find(badImgItr->second.begin(), badImgItr->second.end(),
              (*stagePos)[ievt]) 
            != badImgItr->second.end()) {
        remove_inds.push_back(ievt);
      }
    }
  }

  sort(remove_inds.begin(), remove_inds.end(), greater<int>());
  for (auto ind : remove_inds) {
    scans.erase(scans.begin()+ind);
    imgNums->erase(imgNums->begin()+ind);
    stagePos->erase(stagePos->begin()+ind);
    timeStamps->erase(timeStamps->begin()+ind);
    imgNorms->erase(imgNorms->begin()+ind);
    I0norms->erase(I0norms->begin()+ind);
    imgIsRefs->erase(imgIsRefs->begin()+ind);
    centerRs->erase(centerRs->begin()+ind);
    centerCs->erase(centerCs->begin()+ind);
    centerRstdRatios->erase(centerRstdRatios->begin()+ind);
    centerCstdRatios->erase(centerCstdRatios->begin()+ind);
    I0centerRs->erase(I0centerRs->begin()+ind);
    I0centerCs->erase(I0centerCs->begin()+ind);
    imgRadSTDs->erase(imgRadSTDs->begin()+ind);
    legCoeffs->erase(legCoeffs->begin()+ind);
    legCoeffs_nanMaps->erase(legCoeffs_nanMaps->begin()+ind);
    azmAvgs->erase(azmAvgs->begin()+ind);
    azmAvg_nanMaps->erase(azmAvg_nanMaps->begin()+ind);
    pressures->erase(pressures->begin()+ind);
    UVcounts->erase(UVcounts->begin()+ind);
  }
  std::cerr << "WARNING: RESHAPING LEGENDRES (HACK) SINCE PREPROC DOESN'T CALC THEM\n";
  for (uint i=0; i<legCoeffs->size(); i++) {
    (*legCoeffs)[i].resize(merge.NradLegBins, 0);
    (*legCoeffs_nanMaps)[i].resize(merge.NradLegBins, 0);
  }

  Nentries = scans.size();
  std::vector< std::vector<double> > rpVals(merge.NradAzmBins/5);
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    // Remove Q regions of certain stage positions
    merge.removeBadRegions(&(*azmAvg_nanMaps)[ievt], (*stagePos)[ievt]);

    bool fillLabParams = false;

    // Adding Reference Images
    if (std::find(merge.refSubtractStagePos.begin(), 
            merge.refSubtractStagePos.end(), (*stagePos)[ievt]) 
          != merge.refSubtractStagePos.end()) {
      if (merge.verbose) std::cout << "INFO: Adding reference image.\n";
   
      merge.addReference(
          scans[ievt], (*stagePos)[ievt], (*timeStamps)[ievt],
          &(*azmAvgs)[ievt], &(*azmAvg_nanMaps)[ievt], 
          &(*legCoeffs)[ievt], &(*legCoeffs_nanMaps)[ievt], (*imgNorms)[ievt]);
      fillLabParams = true;
    }

    // Adding Time Dependent Images
    if (!(*imgIsRefs)[ievt]) {
      merge.addEntry(
          scans[ievt], (*stagePos)[ievt], (*timeStamps)[ievt],
          &(*azmAvgs)[ievt], &(*azmAvg_nanMaps)[ievt],
          &(*legCoeffs)[ievt], &(*legCoeffs_nanMaps)[ievt], (*imgNorms)[ievt]);
      //if (stagePos == 1542750) {
      //  int k = 0;
      //  for (int i=4; i<merge.NradAzmBins; i+=5) {
      //    rpVals[k].insert(rpVals[k].end(), 
      //        (*radPixDist)[i].begin(), 
      //        (*radPixDist)[i].end());
      //    k+=1;
      //  }
      //  if (k != rpVals.size()) {
      //    cout<<"SIZES: "<<k<<" "<<rpVals.size()<<endl;
      //  }
      //}
      fillLabParams = true;
    }
    
    //if (stagePos == 1542750) {
    //  plt.print1d((*azmAvg), "ttestLO_1542750_" + to_string(scan));
    //}

    // Adding image parameters
    merge.addLabTimeParameter((*timeStamps)[ievt], "scan", scans[ievt]);
    merge.addLabTimeParameter((*timeStamps)[ievt], "stagePos", (*stagePos)[ievt]);
    merge.addLabTimeParameter((*timeStamps)[ievt], "imgNorm", (*imgNorms)[ievt]);
    merge.addLabTimeParameter((*timeStamps)[ievt], "I0norm", (*I0norms)[ievt]);
    merge.addLabTimeParameter((*timeStamps)[ievt], "centerR", (*centerRs)[ievt]);
    merge.addLabTimeParameter((*timeStamps)[ievt], "centerC", (*centerCs)[ievt]);
    merge.addLabTimeParameter((*timeStamps)[ievt], "I0centerR", (*I0centerRs)[ievt]);
    merge.addLabTimeParameter((*timeStamps)[ievt], "I0centerC", (*I0centerCs)[ievt]);
    merge.addLabTimeParameter((*timeStamps)[ievt], "centerRstdRatio", (*centerRstdRatios)[ievt]);
    merge.addLabTimeParameter((*timeStamps)[ievt], "centerCstdRatio", (*centerCstdRatios)[ievt]);
    // FIX preproc calculation then uncomment
    for (uint i=0; i<(*imgRadSTDs)[ievt].size(); i++) {
      merge.addLabTimeParameter(
          (*timeStamps)[ievt], 
          "imgRadSTD" + to_string(i), 
          (*imgRadSTDs)[ievt][i]);
    }
    // Adding lab/setup parameters
    for (auto const & pvInd : merge.pvMap) {
      if (pvInd.first.compare("pressure") == 0) {
        merge.addLabTimeParameter((*timeStamps)[ievt], "pressure", (*pressures)[ievt]);
      }
      else if (pvInd.first.compare("UVcounts") == 0) {
        merge.addLabTimeParameter((*timeStamps)[ievt], "UVcounts", (*UVcounts)[ievt]);
      }
      else if (pvInd.first.compare("bunkerTemp") == 0) {
        merge.addLabTimeParameter((*timeStamps)[ievt], "bunkerTemp", (*bunkerTemps)[ievt]);
      }
      else if (pvInd.first.compare("highBayTemp") == 0) {
        merge.addLabTimeParameter((*timeStamps)[ievt], "highBayTemp", (*bunkerTemps)[ievt]);
      }
      else {
        std::cerr << "ERROR: Do not know how to add PV " 
          << pvInd.first << "!!!\n";
        exit(1);
      }
    }

    if (merge.saveMergeIntermediates) {
      save::saveDat<double>(
          (*imgOrigs)[ievt],
          merge.saveMergeInterFolder 
            + "/imgOrig-" + merge.run + "_scan-" + to_string(scans[ievt])
            + "_stagePos-" + to_string((*stagePos)[ievt]) + ".dat");
      save::saveDat<double>(
          (*imgSubBkgs)[ievt],
          merge.saveMergeInterFolder 
            + "/imgSubBkg-" + merge.run + "_scan-" + to_string(scans[ievt])
            + "_stagePos-" + to_string((*stagePos)[ievt]) + ".dat");
    }

  }

  ///  Check if I0 moved within this time frame  ///
  for (auto& itr : merge.labTimeParams) {
    labTimeStamps.push_back(itr.first);
  }

  long int tm, tm1;
  double centShift, intsChange;
  std::vector<double> tempGasShift(labTimeStamps.size());
  std::vector< std::vector<double> > plotI0Shift(2);
  plotI0Shift[0].resize(labTimeStamps.size());
  plotI0Shift[1].resize(labTimeStamps.size());
  std::vector< std::vector<double> > plotShift(2);
  plotShift[0].resize(labTimeStamps.size());
  plotShift[1].resize(labTimeStamps.size());
  for (int itm=0; itm<(int)labTimeStamps.size() - 1; itm++) {
    tm  = labTimeStamps[itm];
    tm1 = labTimeStamps[itm+1];
    centShift = sqrt(
        std::pow(
            merge.labTimeParams[tm]["I0centerR"]
            - merge.labTimeParams[tm1]["I0centerR"],
        2)
        + std::pow(
            merge.labTimeParams[tm]["I0centerC"]
            - merge.labTimeParams[tm1]["I0centerC"],
        2));

    plotI0Shift[0][itm] = merge.labTimeParams[tm]["I0centerR"];
    plotI0Shift[0][itm+1] = merge.labTimeParams[tm1]["I0centerR"];
    plotI0Shift[1][itm] = merge.labTimeParams[tm]["I0centerC"];
    plotI0Shift[1][itm+1] = merge.labTimeParams[tm1]["I0centerC"];
    plotShift[0][itm] = merge.labTimeParams[tm]["centerR"];
    plotShift[0][itm+1] = merge.labTimeParams[tm1]["centerR"];
    plotShift[1][itm] = merge.labTimeParams[tm]["centerC"];
    plotShift[1][itm+1] = merge.labTimeParams[tm1]["centerC"];
    intsChange = std::abs(
        merge.labTimeParams[tm]["imgNorm"]
        - merge.labTimeParams[tm1]["imgNorm"])
        /merge.labTimeParams[tm]["imgNorm"];
    
    tempGasShift[itm] = intsChange*centShift;
  }
  tempGasShift[labTimeStamps.size()-1] = 1;
  std::vector<TH1*> hists(2);
  hists[0] = plt.print1d(plotI0Shift[0], "plots/data-" + merge.run + "_I0centerR");
  hists[1] = plt.print1d(plotI0Shift[1], "plots/data-" + merge.run + "_I0centerC");
  hists[0] = plt.print1d(plotShift[0], "plots/data-" + merge.run + "_centerR");
  hists[1] = plt.print1d(plotShift[1], "plots/data-" + merge.run + "_centerC");

  // Want to remove both images before and after nozzle changed
  //for (int itm=0; itm<(int)labTimeStamps.size() - 1; itm++) {
  //  tm  = labTimeStamps[itm];
  //  merge.labTimeParams[tm]["gasShift"] = std::max(
  //      tempGasShift[itm],
  //      tempGasShift[itm+1]);
  //}
  //save::saveDat<double>(tempGasShift,
  //    merge.saveMergeInterFolder
  //    + "/gasShift_run-" + run + "_distribution["
  //    + to_string(tempGasShift.size()) + "].dat");

  ///  Cut on lab time parameters  ///
  
  if (merge.verbose) {
    std::cout << "INFO: Cutting on lab time parameters.\n";
  }

  merge.basicGreaterThanCut("centerRstdRatio", 3);
  merge.basicGreaterThanCut("centerCstdRatio", 3);
  //merge.basicGreaterThanCut("gasShift", merge.gasShiftCut);
  merge.stdParamCut("imgNorm", 3);
  //merge.stdParamCut("imgNorm", 3);
  //merge.stdParamCut("imgNorm", 3);
  merge.stdParamCut("I0norm", 3);
  for (uint i=0; i<imgRadSTDs->size(); i++) {
    merge.stdParamCut("imgRadSTD" + to_string(i), 3);
  }



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
  int Nnans = 0;
  for (int i=0; i<merge.azimuthalAvg.size(); i++) {
    for (int j=0; j<merge.azimuthalAvg[i].size(); j++) {
      if (std::isnan(merge.azimuthalAvg[i][j])) Nnans++;
    }
  }


  // Normalize line outs
  //merge.normalizeScansResults();

  // Get Mean, STD, and SEM
  merge.getRunMeanSTDSEM();


  /////  Saving  /////

  // Clean up NANVAL for saving and plotting
  for (int ir=0; ir<(int)merge.azmReference.size(); ir++) {
    if (merge.azmReference[ir] == NANVAL) {
      merge.azmReference[ir] = 0;
    }
  }
  for (int ilg=0; ilg<(int)merge.legReference.size(); ilg++) {
    for (int ir=0; ir<(int)merge.legReference[ilg].size(); ir++) {
      if (merge.legReference[ilg][ir] == NANVAL) {
        merge.legReference[ilg][ir] = 0;
      }
    }
  }

  if (merge.verbose) 
    std::cout << "INFO: saving merged references and before t0 subtraction.\n";
  // References
  save::saveDat<double>(merge.azmReference,
      merge.mergeScansOutputDir + 
      "/data-" + merge.run + "-" + prefix+ 
      "referenceAzm[" +
      to_string(merge.NradAzmBins) + "].dat");
  save::saveDat<double>(merge.legReference,
      merge.mergeScansOutputDir + 
      "/data-" + merge.run + "-" + prefix+ 
      "referenceLeg[" + 
      to_string(merge.Nlegendres) + 
      "," + to_string(merge.NradLegBins) + "].dat");

  //for (int i=0; i<rpVals.size(); i++) {
  //  double mean = std::accumulate(rpVals[i].begin(), rpVals[i].end(), 0);
  //  mean /= rpVals[i].size();

  //  double var = 0;
  //  for (int k=0; k<rpVals[i].size(); k++) {
  //    var += std::pow(mean - rpVals[i][k], 2);
  //  }
  //  var /= rpVals[i].size();
  //  var = std::sqrt(var);
  //  var /= std::sqrt((float)rpVals[i].size());

  //  cout<<"SEM check radius "<<i
  //    <<" SEM / STD/sqrt(N) / ratio: "<<merge.runAzmSEM[tmInd][4+i*5]
  //    <<"  "<<var<<"  "<<merge.runAzmSEM[tmInd][4+i*5]/var<<"  "<<rpVals[i].size()<<endl;
  //}


  /////  Time Domain Changes  /////

  // Subtract T0
  if (merge.verbose) std::cout << "INFO: Subtracting T0.\n";
  merge.subtractT0();

  /////  Normalize and get statistics  /////

  // Gaussian smooth
  if (merge.verbose) std::cout << "INFO: Gaussian Smearing Q.\n";
  //merge.gaussianFilterQ();

  // Make pair correlations
  if (merge.verbose) std::cout << "INFO: Making pair correlations.\n";
  //merge.makePairCorrs();

  // Get Mean, STD, and SEM
  if (merge.verbose) std::cout << "INFO: Calculating mean, STD, SEM.\n";
  if (merge.useBootstrapSEM) {
    if (merge.computeBootstrapSEM) {
      merge.bootstrapSEM();

      // Testing the number of bootstrap loops are needed
      if (merge.testMergeNbootStrap) {
        merge.testSEMbootstrap();
      }
    }
  }
  else {
    merge.getRunMeanSTDSEM();
  }

  // Normalize to get sM(s)
  if (merge.verbose) std::cout << "INFO: sMs Normalization.\n";
  merge.sMsNormalize();

  std::vector<double> smooth, smoothSMS;
  for (int it=0; it<(int)merge.azimuthalAvg.size(); it++) {
    smooth    = imgProc::gaussianSmooth1d(
                  merge.azimuthalAvg[it], 
                  merge.mergeGSmoothSTD, 
                  7*merge.mergeGSmoothSTD);
    smoothSMS = imgProc::gaussianSmooth1d(
                  merge.azimuthalsMs[it],
                  merge.mergeGSmoothSTD, 
                  7*merge.mergeGSmoothSTD);
    for (int iq=0; iq<merge.NradAzmBins; iq++) {
      merge.azimuthalAvg[it][iq] = smooth[iq];
      merge.azimuthalsMs[it][iq] = smoothSMS[iq];
    }
  }


  // Clean up NANVAL for saving and plotting
  for (int itm=0; itm<(int)merge.azimuthalAvg.size(); itm++) {
    for (int ir=0; ir<(int)merge.azimuthalAvg[itm].size(); ir++) {
      if (merge.azimuthalAvg[itm][ir] == NANVAL) merge.azimuthalAvg[itm][ir] = 0;
    }
  }

  // Smooth in time
  merge.smearTimeGaussian();
  //merge.smearTimeFFT();


  /////  Saving and plotting  /////
  if (merge.verbose) 
    std::cout << "INFO: saving merged references after t0 subtraction.\n";

  if (!(merge.computeBootstrapSEM && merge.SEMisBootstrap) 
      && !(!merge.useBootstrapSEM && !merge.SEMisBootstrap)) {
    std::cerr << "WARNING: Cannot save SEM because method and "
      << "SEMisBootstrap do not align!!!\n";
  }

  if (merge.pltVerbose || true) {
    std::vector<double> test(merge.azimuthalsMs[0].size());
    plt.print1d(merge.azmReference, "./plots/reference");
    plt.printRC(merge.azimuthalAvg, "./plots/azmAvg");
    if (merge.didSMSnormalize) {
      plt.printRC(merge.azimuthalsMs, "./plots/azmAvgSMS");
    for (uint i=0; i<test.size(); i++) {
      test[i] = merge.azimuthalAvg[27][i]/merge.sMsAzmNorm[i];
      if (i<50) {
        test[i] = 0;
      }
    }
    plt.print1d(test, "test1d");
    }
  }

  save::saveDat<double>(merge.azimuthalAvg, 
      merge.mergeScansOutputDir + "data-"
      + merge.run + "-" + prefix + "azmAvgDiff["
      + to_string(merge.azimuthalAvg.size()) + ","
      + to_string(merge.azimuthalAvg[0].size()) + "].dat");
  if (merge.computeBootstrapSEM || !merge.useBootstrapSEM) {
    save::saveDat<double>(merge.runAzmSEM,
        merge.mergeScansOutputDir + "data-"
        + merge.run + "-" + prefix + "azmAvgSEM[" 
        + to_string(merge.runAzmSEM.size()) + ","
        + to_string(merge.runAzmSEM[0].size()) + "].dat");
  }

  save::saveDat<double>(merge.azmReference,
      merge.mergeScansOutputDir + "data-"
      + merge.run + "-" + prefix + "azmReference[" 
      + to_string(merge.azmReference.size()) + "].dat");
  for (auto & aItr : merge.azmIndReference) {
    save::saveDat<double>(aItr.second,
        merge.mergeScansOutputDir + "data-"
        + merge.run + "_" + prefix + "reference-"
        + to_string(aItr.first) + "_bins[" 
        + to_string(merge.NradAzmBins) + "].dat");
  }

  if (merge.smearedTime) {
    save::saveDat<double>(merge.smearedAzmAvg, 
        merge.mergeScansOutputDir + "data-"
        + merge.run + "-" + prefix + "tSmeared-"
        + to_string(merge.timeSmearSTD) + "-azmAvgDiff["
        + to_string(merge.smearedAzmAvg.size()) + ","
        + to_string(merge.smearedAzmAvg[0].size()) + "].dat");
    
    if (merge.didSMSnormalize) {
      save::saveDat<double>(merge.smearedAzmsMs, 
        merge.mergeScansOutputDir + "data-"
        + merge.run + "-" + prefix + "tSmeared-"
        + to_string(merge.timeSmearSTD) + "-sMsAzmAvgDiff["
        + to_string(merge.smearedAzmAvg.size()) + ","
        + to_string(merge.smearedAzmAvg[0].size()) + "].dat");
    }
  }

  if (merge.didSMSnormalize) {
    save::saveDat<double>(merge.azimuthalsMs, 
        merge.mergeScansOutputDir + "data-"
        + merge.run + "-" + prefix + "sMsAzmAvgDiff["
        + to_string(merge.azimuthalAvg.size()) + ","
        + to_string(merge.azimuthalAvg[0].size()) + "].dat");
    if (merge.computeBootstrapSEM || !merge.useBootstrapSEM) {
      save::saveDat<double>(merge.runsMsMeans,
          merge.mergeScansOutputDir + "data-"
          + merge.run + "-" + prefix + "sMsMean[" 
          + to_string(merge.runAzmMeans.size()) + ","
          + to_string(merge.runAzmMeans[0].size()) + "].dat");

      save::saveDat<double>(merge.runsMsSEM,
          merge.mergeScansOutputDir + "data-"
          + merge.run + "-" + prefix + "sMsSEM[" 
          + to_string(merge.runAzmSEM.size()) + ","
          + to_string(merge.runAzmSEM[0].size()) + "].dat");

      save::saveDat<double>(merge.runsMsRefMean,
          merge.mergeScansOutputDir + "data-"
          + merge.run + "-" + prefix
          + "referenceAzmsMsMean[" 
          + to_string(merge.NradAzmBins) + "].dat");

      save::saveDat<double>(merge.runsMsRefSEM,
          merge.mergeScansOutputDir + "data-"
          + merge.run + "-" + prefix
          + "referenceAzmsMsSEM[" 
          + to_string(merge.NradAzmBins) + "].dat");
    }
  }
 
  if (merge.didPairCorrSTD) {
    if (merge.computeBootstrapSEM || !merge.useBootstrapSEM) {
      save::saveDat<double>(merge.runPCorrSEM,
          merge.mergeScansOutputDir + "data-"
          + merge.run + "-" + prefix + "pairCorrSEM["
          + to_string(merge.runPCorrSEM.size()) + ","
          + to_string(merge.runPCorrSEM[0].size()) + "].dat");
    }
  }

  if (merge.saveMergeIntermediates) {
    for (auto & sRitr : merge.scanReferences) {
      for (auto & pItr : sRitr.second) {
        if (pItr.second.scale) {
          save::saveDat<double>(
              pItr.second.azmRef, 
              merge.saveMergeInterFolder 
                + "/azmAvgRef-" + merge.run + "_scan-" + to_string(sRitr.first)
                + "_stagePos-" + to_string(pItr.first) + ".dat");
        }
      }
    }

    for (auto & sAzml : merge.scanAzmAvg) {
      for (auto & pItr : merge.stagePosInds) {
        if (merge.scanScale[sAzml.first][pItr.second] != 0) {
          save::saveDat<double>(
              sAzml.second[pItr.second],
              merge.saveMergeInterFolder 
                + "/azmAvg-" + merge.run + "_scan-" + to_string(sAzml.first)
                + "_stagePos-" + to_string(pItr.first) + ".dat");
        }
      }
    }
  }


  if (scanSearch) {
    plt.printRC(merge.azimuthalAvg, 
        merge.mergeScansOutputDir + "/plots/data-"
        + merge.run + "-" + prefix + "azmAvgDiff", opts, vals);
  }


  /////////////////////////
  /////  Cleaning up  /////
  /////////////////////////

  if (merge.verbose) {
    std::cout << "Cleaning up" << endl;
  }

  return 0;
}
