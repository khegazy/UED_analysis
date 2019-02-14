#include "../analysis.h" 


using namespace std;


int main(int argc, char* argv[]) {

  if (argc<2) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe 'runName' !!!"<<endl;
    exit(0);
  }


  ///// Load environment and get the number of events /////

  string runName(argv[1]);
  bool findBestTheoryMatch = false;

  // Get parameters
  cout<<"RunName: "<<runName<<endl;
  parameterClass params(runName);


  // Import Data
  std::string diffFileName = "data-" + runName + "-azmAvgDiff";
  std::vector<int> shape  = save::getShape(
                                params.mergeScansOutputDir, 
                                diffFileName);
  diffFileName = 
      params.mergeScansOutputDir 
      + diffFileName + "["
      + to_string(shape[0]) + ","
      + to_string(shape[1]) + "].dat";
  std::string sMsFileName = 
      params.mergeScansOutputDir 
      + "data-" + runName + "-sMsAzmAvgDiff["
      + to_string(shape[0]) + ","
      + to_string(shape[1]) + "].dat";
  std::string pairCorrFileName = 
      "./results/data-" + runName
      + "_pairCorrOdd[" + to_string(shape[0]) 
      + ","  + to_string(params.maxRbins) + "].dat";
  std::vector< std::vector<double> > tmDpDiffraction(shape[0]);
  std::vector< std::vector<double> > tmDpSMS(shape[0]);
  std::vector< std::vector<double> > tmDpPairCorr(shape[0]);
  for (int it=0; it<shape[0]; it++) {
     tmDpDiffraction[it].resize(shape[1], 0);
     tmDpSMS[it].resize(shape[1], 0);
     tmDpPairCorr[it].resize(params.maxRbins, 0);
  }
  save::importDat<double>(tmDpDiffraction, diffFileName);
  save::importDat<double>(tmDpSMS, sMsFileName);
  save::importDat<double>(tmDpPairCorr, pairCorrFileName);

  // Import reference
  std::string simRefName =
      params.simOutputDir + "nitrobenzene_molDiffractionPatternLineOut_Qmax-"
      + to_string(params.maxQazm) + "_Ieb-"
      + to_string(params.Iebeam) + "_scrnD-"
      + to_string(params.screenDist) + "_elE-"
      + to_string(params.elEnergy) + "_Bins["
      + to_string(shape[1]) + "].dat";

  std::string simSMSrefName =
      params.simOutputDir + "nitrobenzene_sMsPatternLineOut_Qmax-"
      + to_string(params.maxQazm) + "_Ieb-"
      + to_string(params.Iebeam) + "_scrnD-"
      + to_string(params.screenDist) + "_elE-"
      + to_string(params.elEnergy) + "_Bins["
      + to_string(shape[1]) + "].dat";


  std::vector<double> reference(shape[1]);
  std::vector<double> sMsReference(shape[1]);
  save::importDat<double>(reference, simRefName);
  save::importDat<double>(sMsReference, simSMSrefName);
  plt.print1d(sMsReference, "ref");

  plt.print1d(reference, "refTest");
  plt.print1d(sMsReference, "refsMsTest");

  // Getting final state from data
  std::vector<TH1*> hists(2);
  std::vector<PLOToptions> opts(2);
  std::vector<std::string> vals(2);
  opts[0] = xLabel;   vals[0] = "Q [#AA^{-1}]";
  opts[1] = xSpan;    vals[1] = to_string(0) + "," + to_string(params.maxQazm);
  std::vector<double> finalState(params.NradAzmBins, 0);
  std::vector<double> sMsFinalState(params.NradAzmBins, 0);
  ///// Loop through events in the file /////
  for (int it=shape[0]-1; it>shape[0]-1-params.NfinalPoints; it--) {
    //if (it == 22 || it == 24) {
    //  continue;
    //}
    std::transform(tmDpDiffraction[it].begin(), tmDpDiffraction[it].end(),
        finalState.begin(), finalState.begin(),
        std::plus<double>());
    std::transform(tmDpSMS[it].begin(), tmDpSMS[it].end(),
        sMsFinalState.begin(), sMsFinalState.begin(),
        std::plus<double>());
  //hists[0] = plt.print1d(simFinalState, "simFinalDiff");
  //hists[0] = plt.print1d(simFinalState, "simFinalDiff");
  //hists[1] = plt.print1d(tmDpDiffraction[it], "dataFinalDiff");
  //plt.print1d(hists, "./compareFinalStateDiffraction_" + to_string(it), opts, vals);
  }
  delete hists[0];
  delete hists[1];

  for (int iq=0; iq<params.NradAzmBins; iq++) {
    finalState[iq] /= params.NfinalPoints;
    sMsFinalState[iq] /= params.NfinalPoints;
  }

  std::string finalStateDataFileName = 
      "data-" + runName 
      + "_diffFinalState["
      + to_string(shape[1]) + "].dat";
  std::string finalStateDataSMSfileName = 
      "data-" + runName 
      + "_sMsFinalState["
      + to_string(shape[1]) + "].dat";

  save::saveDat<double>(finalState, 
      params.mergeScansOutputDir + finalStateDataFileName);
  save::saveDat<double>(sMsFinalState, 
      params.mergeScansOutputDir + finalStateDataSMSfileName);


  // Simulated final states
  int fillQbegin = (int)(shape[1]*params.fsQfitBegin/params.maxQazm);
  int fillQend   = (int)(shape[1]*params.fsQfitEnd/params.maxQazm);
  int fillRbegin = (int)(params.maxRbins*params.fsRfitBegin/params.maxR);
  int fillRend   = (int)(params.maxRbins*params.fsRfitEnd/params.maxR);
  Eigen::Matrix<double, Eigen::Dynamic, 2> Xfit;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Xqs;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Xrs;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Yfit;
  Eigen::MatrixXd fit;
  Eigen::MatrixXd tmDpFit;
  Xqs.resize(fillQend - fillQbegin, params.finalStates.size() + 1);
  Xrs.resize(fillRend - fillRbegin, params.finalStates.size() + 1);
  for (int i=0; i<fillQend - fillQbegin; i++) {
    Xqs(i,0) = 1;
  }
  for (int i=0; i<fillRend - fillRbegin; i++) {
    Xrs(i,0) = 1;
  }
  std::vector<double> sMsFitCoeff(2,0), rawFitCoeff(2,0), pcFitCoeff(2,0);
  std::vector<double> pairCorrFinalState(params.maxRbins);
  std::vector<double> sMsFinalState_scaled(shape[1]);
  std::vector<double> finalState_scaled(shape[1]);
  std::vector<double> pairCorrFinalState_scaled(params.maxRbins);
  std::vector< std::vector<double> > simPairCorrFinalStates(params.finalStates.size());
  std::vector< std::vector<double> > simFinalDiffractions(params.finalStates.size());
  std::vector< std::vector<double> > simFinalSMSs(params.finalStates.size());
  for (int ifs=0; ifs<(int)params.finalStates.size(); ifs++) {
    simFinalDiffractions[ifs].resize(shape[1], 0);
    simFinalSMSs[ifs].resize(shape[1], 0);
    simPairCorrFinalStates[ifs].resize(params.maxRbins, 0);
  }
  for (int ifs=0; ifs<(int)params.finalStates.size(); ifs++) {
    std::string fsName = params.finalStates[ifs];
    std::string simFinalName =
        params.simOutputDir + fsName 
        + "_molDiffractionPatternLineOut_Qmax-" 
        + to_string(params.maxQazm) + "_Ieb-"
        + to_string(params.Iebeam) + "_scrnD-"
        + to_string(params.screenDist) + "_elE-"
        + to_string(params.elEnergy) + "_Bins["
        + to_string(shape[1]) + "].dat";

    std::string simSMSfinalName =
        params.simOutputDir + fsName 
        + "_sMsPatternLineOut_Qmax-" 
        + to_string(params.maxQazm) + "_Ieb-"
        + to_string(params.Iebeam) + "_scrnD-"
        + to_string(params.screenDist) + "_elE-"
        + to_string(params.elEnergy) + "_Bins["
        + to_string(shape[1]) + "].dat";


    save::importDat<double>(simFinalDiffractions[ifs], simFinalName);
    save::importDat<double>(simFinalSMSs[ifs], simSMSfinalName);
   
    // Simulation final state difference diffraction
    std::vector<double> simFinalState(shape[1]);
    std::vector<double> simSMSfinalState(shape[1]);
    for (int iq=0; iq<shape[1]; iq++) {
      simFinalState[iq] = simFinalDiffractions[ifs][iq] - reference[iq];
      simSMSfinalState[iq] = simFinalSMSs[ifs][iq] - sMsReference[iq];
    }
    plt.print1d(simFinalState, fsName+"_simFinTest");
    plt.print1d(simSMSfinalState, fsName+"_simsMsFinTest");


    std::string finalStateSimFileName = fsName 
        + "_diffFinalState["
        + to_string(shape[1]) + "].dat";
    save::saveDat<double>(
        simFinalState, 
        params.simOutputDir + finalStateSimFileName); 
    std::string finalStateSimSMSfileName = fsName 
        + "_sMsFinalState["
        + to_string(shape[1]) + "].dat";
    save::saveDat<double>(
        simSMSfinalState, 
        params.simOutputDir + finalStateSimSMSfileName); 
          
    // Fitting 
    Xfit.resize(fillQend - fillQbegin, 2);
    Yfit.resize(fillQend - fillQbegin, 1);
    for (int iq=fillQbegin; iq<fillQend; iq++) {
      Xfit(iq-fillQbegin,0)   = 1;
      Xfit(iq-fillQbegin,1)   = simSMSfinalState[iq];
      Xqs(iq-fillQbegin,ifs+1)= simSMSfinalState[iq];
      Yfit(iq-fillQbegin,0)   = sMsFinalState[iq];
    }

    fit = tools::normalEquation(Xfit, Yfit);
    for (int iq=0; iq<shape[1]; iq++) {
      sMsFinalState_scaled[iq] = simSMSfinalState[iq]*fit(1) + fit(0);
    }

    save::saveDat<double>(sMsFinalState_scaled, 
        "./results/sim-" + fsName 
          + "_sMsFinalState_scaled[" 
          + to_string(shape[1]) + "].dat");

  
    for (int iq=fillQbegin; iq<fillQend; iq++) {
      Xfit(iq-fillQbegin,0) = 1;
      Xfit(iq-fillQbegin,1) = simFinalState[iq];
      Yfit(iq-fillQbegin,0) = finalState[iq];
    }

    fit = tools::normalEquation(Xfit, Yfit);
    for (int iq=0; iq<shape[1]; iq++) {
      finalState_scaled[iq] = simFinalState[iq]*fit(1) + fit(0);
    }

    save::saveDat<double>(finalState_scaled, 
        "./results/sim-" + fsName 
          + "_diffFinalState_scaled[" 
          + to_string(shape[1]) + "].dat");


    ///////////////////////////////
    /////  Pair Correlations  /////
    ///////////////////////////////

    // Final State
    if (params.fillLowQtheory) {
      system(("./pairCorr.exe simulateReference -Idir "
          + params.simOutputDir + " -Fname "
          + finalStateSimSMSfileName + " -lowQtheory "
          + params.fillLowQfile + " -MolName "
          + fsName).c_str());
      system(("./pairCorr.exe " + runName + " -Idir "
          + params.mergeScansOutputDir + " -Fname "
          + finalStateDataSMSfileName + " -Osuf _finalState"
          + " -lowQtheory " + params.fillLowQfile).c_str());

    }
    else {
      system(("./pairCorr.exe simulateReference -Idir "
          + params.simOutputDir + " -Fname "
          + finalStateSimSMSfileName + " -MolName "
          + fsName).c_str());
      system(("./pairCorr.exe " + runName + " -Idir "
          + params.mergeScansOutputDir + " -Fname "
          + finalStateDataSMSfileName 
          + " -Osuf _finalState").c_str());
    }

    // Fitting 
    std::string simPairCorrName = 
        "./results/sim-" + fsName 
        + "_pairCorrOdd[" + to_string(params.maxRbins)
        + "].dat";
    
    std::string pairCorrName = 
        "./results/data-" + runName 
        + "_finalState_pairCorrOdd[" 
        + to_string(params.maxRbins) + "].dat";

    save::importDat<double>(simPairCorrFinalStates[ifs], simPairCorrName);
    save::importDat<double>(pairCorrFinalState, pairCorrName);

    Xfit.resize(fillRend - fillRbegin, 2);
    Yfit.resize(fillRend - fillRbegin, 1);
    for (int ir=fillRbegin; ir<fillRend; ir++) {
      Xfit(ir-fillRbegin,0)   = 1;
      Xfit(ir-fillRbegin,1)   = simPairCorrFinalStates[ifs][ir];
      Xrs(ir-fillRbegin,ifs+1)= simPairCorrFinalStates[ifs][ir];
      Yfit(ir-fillRbegin,0)   = pairCorrFinalState[ir];
    }

    fit = tools::normalEquation(Xfit, Yfit);
    for (int ir=0; ir<params.maxRbins; ir++) {
      /// FIXME: check if should be sim not paircorr
      pairCorrFinalState_scaled[ir] = pairCorrFinalState[ir]*fit(1) + fit(0);
    }

    save::saveDat<double>(pairCorrFinalState_scaled, 
        "./results/sim-" + fsName 
          + "_pairCorrFinalState_scaled[" 
          + to_string(params.maxRbins) + "].dat");


    // Searching for best low Q fill
    if (findBestTheoryMatch) {
      std::vector<double> scales = 
          {0.2, 0.3, 0.4, 0.5, 0.6, 
            0.7, 0.8, 0.9, 1, 1.1, 
            1.2, 1.3, 1.4, 1.5};
      std::vector<double> currentSim(params.NradAzmBins);
      for (auto& scl : scales) {
        for (int iq=0; iq<shape[1]; iq++) {
          currentSim[iq] = scl*simSMSfinalState[iq];
        }
      
        std::string currentSimName = 
          "./results/sim-" + fsName + "_LowQfill_scale-"
            + to_string(scl) + "_Bins["
            + to_string(params.NradAzmBins) + "].dat";
        save::saveDat<double>(currentSim, currentSimName);

        system(("./pairCorr.exe " + runName 
            + " -Osuf _finalStateDiff -lowQtheory "
            + currentSimName).c_str());

        system(("mv ./results/data-" + runName 
            + "_finalStateDiff-pairCorrEven[" 
            + to_string(shape[0]) + ","
            + to_string(params.maxRbins) 
            + "].dat ./results/data-" + runName
            + "_lowQscale-" 
            + to_string(scl) + "_pairCorrEven["
            + to_string(shape[0]) + ","
            + to_string(params.maxRbins) + "].dat").c_str());

        system(("mv ./results/data-" + runName 
            + "_finalStateDiff-pairCorrOdd[" 
            + to_string(shape[0]) + ","
            + to_string(params.maxRbins) 
            + "].dat ./results/data-" + runName
            + "_lowQscale-" 
            + to_string(scl) + "_pairCorrOdd["
            + to_string(shape[0]) + ","
            + to_string(params.maxRbins) + "].dat").c_str());


        //system(("rm " + currentSimName).c_str());
      }

      for (int iq=0; iq<shape[1]; iq++) {
        currentSim[iq] = params.lowQfillSimScale*simFinalState[iq];
      }
      save::saveDat<double>(currentSim, 
          "./results/sim-phenoxyRadicalLowQfill["
          + to_string(params.NradAzmBins) + "].dat");
    }
  }



  cout<<"starting fitting linear comb\n";
  ////////////////////////////////////////////////////////////
  /////  Fitting Linear Combination of All Final States  /////
  ////////////////////////////////////////////////////////////

  /////  Final State  /////

  ///  SMS  ///
  Yfit.resize(fillQend - fillQbegin, 1);
  for (int iq=fillQbegin; iq<fillQend; iq++) {
    Yfit(iq-fillQbegin,0)   = sMsFinalState[iq];
  }

  fit = tools::normalEquation(Xqs, Yfit);
  std::vector<double> combinedQ(shape[1]);
  for (int iq=0; iq<shape[1]; iq++) {
    combinedQ[iq] += fit(0);
    for (int ifs=0; ifs<(int)params.finalStates.size(); ifs++) {
      combinedQ[iq] += simFinalSMSs[ifs][iq]*fit(ifs+1);
    }
  }

  save::saveDat<double>(combinedQ, 
      "./results/sim-" + runName
      + "_sMsFinalState_scaledLinComb_Bins["
      + to_string(combinedQ.size()) + "].dat");


  ///  Pair Correlation  ///
  Yfit.resize(fillRend - fillRbegin, 1);
  for (int ir=fillRbegin; ir<fillRend; ir++) {
    Yfit(ir-fillRbegin,0)   = pairCorrFinalState[ir];
  }

  fit = tools::normalEquation(Xrs, Yfit);
  std::vector<double> combinedR(params.maxRbins);
  for (int ir=0; ir<params.maxRbins; ir++) {
    combinedR[ir] += fit(0);
    for (int ifs=0; ifs<(int)params.finalStates.size(); ifs++) {
      combinedR[ir] += simPairCorrFinalStates[ifs][ir]*fit(ifs+1);
    }
  }

  save::saveDat<double>(combinedR, 
      "./results/sim-" + runName
      + "_pairCorrFinalState_scaledLinComb_Bins["
      + to_string(combinedR.size()) + "].dat");


  cout<<"starting time dep"<<endl;
  /////  Time Dependent  /////
  Yfit.resize(fillQend - fillQbegin, 1);
  std::vector< std::vector<double> > combinedQs(tmDpSMS.size());
  std::vector< std::vector<double> > tmDpQfits(params.finalStates.size()+1);
  for (int i=0; i<(int)tmDpQfits.size(); i++) {
    tmDpQfits[i].resize(tmDpSMS.size());
  }
  for (int it=0; it<(int)tmDpSMS.size(); it++) {
    combinedQs[it].resize(params.NradAzmBins, 0);
    for (int iq=fillQbegin; iq<fillQend; iq++) {
      Yfit(iq-fillQbegin,0) = tmDpSMS[it][iq];
    }

    tmDpFit = tools::normalEquation(Xqs, Yfit);

    for (int ifs=0; ifs<(int)tmDpQfits.size(); ifs++) {
      tmDpQfits[ifs][it] = tmDpFit(ifs);
    }

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      combinedQs[it][iq] += tmDpFit(0);
      for (int ifs=0; ifs<(int)params.finalStates.size(); ifs++) {
        combinedQs[it][iq] += simFinalSMSs[ifs][iq]*tmDpFit(ifs+1);
      }
    }
  }
  save::saveDat<double>(tmDpQfits[0], 
      "./results/sim-" + runName
        + "-offset_sMsFitCoeffsLinComb_Bins[" 
        + to_string(tmDpSMS.size()) + "].dat");

  for (int ifs=0; ifs<(int)params.finalStates.size(); ifs++) {
    save::saveDat<double>(tmDpQfits[ifs+1], 
        "./results/sim-" + runName + "-" + params.finalStates[ifs]
          + "_sMsFitCoeffsLinComb_Bins[" 
          + to_string(tmDpSMS.size()) + "].dat");
  }

  save::saveDat<double>(combinedQs, 
      "./results/sim-" + runName
      + "_sMsFinalStateFitLinComb_Bins["
      + to_string(combinedQs.size()) + ","
      + to_string(combinedQs[0].size()) + "].dat");


cout<<"starting paircorr"<<endl;

  Yfit.resize(fillRend - fillRbegin, 1);
  std::vector< std::vector<double> > combinedRs(tmDpPairCorr.size());
  std::vector< std::vector<double> > tmDpRfits(params.finalStates.size()+1);
  for (int i=0; i<(int)tmDpRfits.size(); i++) {
    tmDpRfits[i].resize(tmDpPairCorr.size());
  }
  for (int it=0; it<(int)tmDpPairCorr.size(); it++) {
    cout<<"init"<<endl;
    combinedRs[it].resize(params.NradAzmBins, 0);
    for (int ir=fillRbegin; ir<fillRend; ir++) {
      Yfit(ir-fillRbegin,0) = tmDpPairCorr[it][ir];
    }

    tmDpFit = tools::normalEquation(Xrs, Yfit);

    for (int ifs=0; ifs<(int)tmDpRfits.size(); ifs++) {
      tmDpRfits[ifs][it] = tmDpFit(ifs);
    }

    for (int ir=0; ir<params.maxRbins; ir++) {
      combinedRs[it][ir] += tmDpFit(0);
      for (int ifs=0; ifs<(int)params.finalStates.size(); ifs++) {
        combinedRs[it][ir] += simPairCorrFinalStates[ifs][ir]*tmDpFit(ifs+1);
      }
    }
  }

  save::saveDat<double>(tmDpRfits[0], 
      "./results/sim-" + runName
        + "-offset_pairCorrFitCoeffsLinComb_Bins[" 
        + to_string(tmDpPairCorr.size()) + "].dat");

  for (int ifs=0; ifs<(int)params.finalStates.size(); ifs++) {
    save::saveDat<double>(tmDpRfits[ifs+1], 
        "./results/sim-" + runName + "-" +params.finalStates[ifs] 
          + "_pairCorrFitCoeffsLinComb_Bins[" 
          + to_string(tmDpSMS.size()) + "].dat");
  }

  save::saveDat<double>(combinedRs, 
      "./results/sim-" + runName 
      + "_pairCorrFinalStateFitLinComb_Bins["
      + to_string(combinedRs.size()) + ","
      + to_string(combinedRs[0].size()) + "].dat");







  //delete hists[0];
  //delete hists[1];
  return 1;
}
