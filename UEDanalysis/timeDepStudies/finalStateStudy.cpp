#include "../analysis.h" 


using namespace std;


int main(int argc, char* argv[]) {

  if (argc<2) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe 'fileList.txt' !!!"<<endl;
    cerr<<"         Can also run ./analysis.exe 'fileList.txt' 'treeName'"<<endl;
    exit(0);
  }


  ///// Load environment and get the number of events /////

  string filePath(argv[1]);
  auto sPos = filePath.find_last_of("/");
  std::string folder = filePath.substr(0, sPos);
  std::string fileName = filePath.substr(sPos+1, sPos-filePath.size());

  bool findBestTheoryMatch = false;

  // Get run name
  auto iPos = filePath.find("-20");
  iPos += 1;
  auto fPos = filePath.find("-", iPos);
  std::string runName = filePath.substr(iPos, fPos - iPos);
  cout<<"RunName: "<<runName<<endl;
  parameterClass params(runName);


  // Import Data
  std::vector<int> shape = save::getShape(folder, fileName);
  std::vector< std::vector<double> > tdDiffraction(shape[0]);
  for (int it=0; it<shape[0]; it++) {
     tdDiffraction[it].resize(shape[1], 0);
  }
  save::importDat<double>(tdDiffraction, filePath);

  // Import reference
  std::string simRefName =
      params.simOutputDir + "nitrobenzene_sMsPatternLineOut_Qmax-"
      + to_string(params.maxQazm) + "_Ieb-"
      + to_string(params.Iebeam) + "_scrnD-"
      + to_string(params.screenDist) + "_elE-"
      + to_string(params.elEnergy) + "_Bins["
      + to_string(shape[1]) + "].dat";

  std::vector<double> reference(shape[1]);
  save::importDat<double>(reference, simRefName);

  // Import final state
  std::string simFinalName =
      params.simOutputDir + params.finalState 
      + "_sMsPatternLineOut_Qmax-" 
      + to_string(params.maxQazm) + "_Ieb-"
      + to_string(params.Iebeam) + "_scrnD-"
      + to_string(params.screenDist) + "_elE-"
      + to_string(params.elEnergy) + "_Bins["
      + to_string(shape[1]) + "].dat";

  std::vector<double> simFinalDiffraction(shape[1]);
  save::importDat<double>(simFinalDiffraction, simFinalName);
  
  std::vector<double> simFinalState(shape[1]);
  for (int iq=0; iq<shape[1]; iq++) {
    simFinalState[iq] = simFinalDiffraction[iq] - reference[iq];
  }

  std::string finalStateSimFname = params.finalState 
      + "_sMsFinalState["
      + to_string(shape[1]) + "].dat";
  save::saveDat<double>(simFinalState, params.simOutputDir + finalStateSimFname); 
        
 

  std::vector<TH1*> hists(2);
  std::vector<PLOToptions> opts(2);
  std::vector<std::string> vals(2);
  opts[0] = xLabel;   vals[0] = "Q [#AA^{-1}]";
  opts[1] = xSpan;    vals[1] = to_string(0) + "," + to_string(params.maxQazm);
  std::vector<double> finalState(params.NradAzmBins, 0);
  ///// Loop through events in the file /////
  for (int it=shape[0]-1; it>shape[0]-1-params.NfinalPoints; it--) {
    //if (it == 22 || it == 24) {
    //  continue;
    //}
    std::transform(tdDiffraction[it].begin(), tdDiffraction[it].end(),
        finalState.begin(), finalState.begin(),
        std::plus<double>());
  //hists[0] = plt.print1d(simFinalState, "simFinalDiff");
  //hists[1] = plt.print1d(tdDiffraction[it], "dataFinalDiff");
  //plt.print1d(hists, "./compareFinalStateDiffraction_" + to_string(it), opts, vals);
  }
  delete hists[0];
  delete hists[1];

  for (int iq=0; iq<params.NradAzmBins; iq++) {
    finalState[iq] /= params.NfinalPoints;
  }

  std::string finalStateDataFname = 
      "data-" + runName 
      + "-sMsFinalState["
      + to_string(shape[1]) + "].dat";

  save::saveDat<double>(finalState, 
      params.mergeScansOutputDir + finalStateDataFname);

  Eigen::Matrix<double, Eigen::Dynamic, 2> Xq;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Yq;
  Xq.resize(shape[1] - params.NbinsSkip, 2);
  Yq.resize(shape[1] - params.NbinsSkip, 1);

  for (int iq=params.NbinsSkip; iq<shape[1]; iq++) {
    Xq(iq-params.NbinsSkip,0) = 1;
    Xq(iq-params.NbinsSkip,1) = simFinalState[iq];
    Yq(iq-params.NbinsSkip,0) = finalState[iq];
  }

  Eigen::MatrixXd fit = tools::normalEquation(Xq, Yq);

  for (int iq=0; iq<shape[1]; iq++) {
    //simFinalState[iq] = simFinalState[iq]*fit(1) + fit(0);
  }

  //save::saveDat<double>(simFinalState, "./results/" + finalStateSimFname);

  //hists[0] = plt.print1d(simFinalState, "simFinalDiff");
  //hists[1] = plt.print1d(finalState, "dataFinalDiff");
  //plt.print1d(hists, "./plots/compareFinalStateDiffraction_" + params.molName + params.finalState, opts, vals);

  ///////////////////////////////
  /////  Pair Correlations  /////
  ///////////////////////////////

  // Final State
  if (params.fillLowQtheory) {
    system(("./pairCorr.exe simulateReference -Idir "
        + params.simOutputDir + " -Fname "
        + finalStateSimFname + " -Osuf "
        + params.finalState + " -lowQtheory "
        + params.fillLowQfile).c_str());
    system(("./pairCorr.exe " + runName + " -Idir "
        + params.mergeScansOutputDir + " -Fname "
        + finalStateDataFname + " -Osuf -finalState"
        + " -lowQtheory " + params.fillLowQfile).c_str());

  }
  else {
    system(("./pairCorr.exe simulateReference -Idir "
        + params.simOutputDir + " -Fname "
        + finalStateSimFname + " -Osuf "
        + params.finalState).c_str());
    system(("./pairCorr.exe " + runName + " -Idir "
        + params.mergeScansOutputDir + " -Fname "
        + finalStateDataFname 
        + " -Osuf -finalState").c_str());
  }

  if (findBestTheoryMatch) {
    std::vector<double> scales = 
        {0.2, 0.3, 0.4, 0.5, 0.6, 
          0.7, 0.8, 0.9, 1, 1.1, 
          1.2, 1.3, 1.4, 1.5};
    std::vector<double> currentSim(params.NradAzmBins);
    std::string currentSimName = "./results/currentLowQSim["
        + to_string(params.NradAzmBins) + "].dat";
    for (auto& scl : scales) {
      for (int iq=0; iq<shape[1]; iq++) {
        currentSim[iq] = scl*simFinalState[iq];
      }
      save::saveDat<double>(currentSim, currentSimName);

      system(("./pairCorr.exe " + runName 
          + " -Osuf -FinalStateDiff -lowQtheory "
          + currentSimName).c_str());

      system(("mv ./results/data-" + runName 
          + "-FinalStateDiff-pairCorrEven[" 
          + to_string(shape[0]) + ","
          + to_string(params.maxRbins) 
          + "].dat ./results/data-" + runName
          + "-lowQscale-" 
          + to_string(scl) + "-pairCorrEven["
          + to_string(shape[0]) + ","
          + to_string(params.maxRbins) + "].dat").c_str());

      system(("mv ./results/data-" + runName 
          + "-FinalStateDiff-pairCorrOdd[" 
          + to_string(shape[0]) + ","
          + to_string(params.maxRbins) 
          + "].dat ./results/data-" + runName
          + "-lowQscale-" 
          + to_string(scl) + "-pairCorrOdd["
          + to_string(shape[0]) + ","
          + to_string(params.maxRbins) + "].dat").c_str());


      system(("rm " + currentSimName).c_str());
    }

    for (int iq=0; iq<shape[1]; iq++) {
      currentSim[iq] = params.lowQfillSimScale*simFinalState[iq];
    }
    save::saveDat<double>(currentSim, 
        "./results/sim-phenoxyRadicalLowQfill["
        + to_string(params.NradAzmBins) + "].dat");
  }

  //delete hists[0];
  //delete hists[1];
  return 1;
}
