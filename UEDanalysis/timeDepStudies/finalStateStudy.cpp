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

  ///  Get parameters  ///
  cout<<"RunName: "<<runName<<endl;
  parameterClass params(runName);
  int mInd;
  int offShift = 0;
  int NfitParams = params.finalStates.size();
  if (params.fsFitOffset) {
    NfitParams += 1;
    offShift = 1;
  }


  ///  Import Data  ///
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

  // Importing standard deviations
  std::string sMsSEMfileName = 
      params.mergeScansOutputDir 
      + "data-" + runName + "-sMsSEM["
      + to_string(shape[0]) + ","
      + to_string(shape[1]) + "].dat";
  std::string pairCorrSEMfileName = 
      params.mergeScansOutputDir 
      + "data-" + runName + "-pairCorrSEM["
      + to_string(shape[0]) + ","
      + to_string(params.maxRbins) + "].dat";


  std::vector< std::vector<double> > tmDpSMSSEM(shape[0]);
  std::vector< std::vector<double> > tmDpSMSvar(shape[0]);
  std::vector< std::vector<double> > tmDpPairCorrSEM(shape[0]);
  std::vector< std::vector<double> > tmDpPairCorrVar(shape[0]);
  for (int it=0; it<shape[0]; it++) {
    tmDpSMSSEM[it].resize(shape[1], 0);
    tmDpSMSvar[it].resize(shape[1], 0);
    tmDpPairCorrSEM[it].resize(params.maxRbins, 0);
    tmDpPairCorrVar[it].resize(params.maxRbins, 0);
  }
  save::importDat<double>(tmDpSMSSEM, sMsSEMfileName);
  save::importDat<double>(tmDpPairCorrSEM, pairCorrSEMfileName);

  for (int it=0; it<shape[0]; it++) {
    for (int iq=0; iq<shape[1]; iq++) {
      tmDpSMSvar[it][iq] = std::pow(tmDpSMSSEM[it][iq], 2);
    }

    for (int ir=0; ir<params.maxRbins; ir++) {
      tmDpPairCorrVar[it][ir] = std::pow(tmDpPairCorrSEM[it][ir], 2);
    }
  }

  // Getting final state from data
  std::vector<TH1*> hists(2);
  std::vector<PLOToptions> opts(2);
  std::vector<std::string> vals(2);
  opts[0] = xLabel;   vals[0] = "Q [#AA^{-1}]";
  opts[1] = xSpan;    vals[1] = to_string(0) + "," + to_string(params.maxQazm);
  std::vector<double> finalState(params.NradAzmBins, 0);
  std::vector<double> sMsFinalState(params.NradAzmBins, 0);
  std::vector<double> sMsFinalStateScaled(params.NradAzmBins, 0);
  std::vector<double> sMsFinalStateVar(params.NradAzmBins, 0);
  std::vector<double> sMsFinalStateVarScaled(params.NradAzmBins, 0);
  std::vector<double> sMsFinalStateSEM(params.NradAzmBins, 0);
  std::vector<double> sMsFinalStateSEMScaled(params.NradAzmBins, 0);
  std::vector<double> pairCorrFinalStateVar(params.maxRbins, 0);
  std::vector<double> pairCorrFinalStateSEM(params.maxRbins, 0);

  ///// Loop through events in the file /////
  for (int it=shape[0]-1; it>shape[0]-1-params.NfinalPoints; it--) {
    std::transform(tmDpDiffraction[it].begin(), tmDpDiffraction[it].end(),
        finalState.begin(), finalState.begin(),
        std::plus<double>());
    std::transform(tmDpSMS[it].begin(), tmDpSMS[it].end(),
        sMsFinalState.begin(), sMsFinalState.begin(),
        std::plus<double>());
    std::transform(tmDpSMSvar[it].begin(), tmDpSMSvar[it].end(),
        sMsFinalStateVar.begin(), sMsFinalStateVar.begin(),
        std::plus<double>());
    std::transform(tmDpPairCorrVar[it].begin(), tmDpPairCorrVar[it].end(),
        pairCorrFinalStateVar.begin(), pairCorrFinalStateVar.begin(),
        std::plus<double>());
    plt.print1d(tmDpSMSvar[it], "testVar_"+to_string(it));
  //hists[0] = plt.print1d(simFinalState, "simFinalDiff");
  //hists[0] = plt.print1d(simFinalState, "simFinalDiff");
  //hists[1] = plt.print1d(tmDpDiffraction[it], "dataFinalDiff");
  //plt.print1d(hists, "./compareFinalStateDiffraction_" + to_string(it), opts, vals);
  }
  delete hists[0];
  delete hists[1];

  for (int iq=0; iq<params.NradAzmBins; iq++) {
    finalState[iq]        /= params.NfinalPoints;
    sMsFinalState[iq]     /= params.NfinalPoints;
    sMsFinalStateVar[iq]  /= std::pow(params.NfinalPoints, 2); // Power of two for SEM
    sMsFinalStateSEM[iq]  = std::sqrt(sMsFinalStateVar[iq]);
    sMsFinalStateScaled[iq]     = sMsFinalState[iq]
                                    *exp(-1*std::pow(iq, 2)/(2*params.fsFilterVar));
    sMsFinalStateSEMScaled[iq]  = sMsFinalStateSEM[iq]
                                    *exp(-1*std::pow(iq, 2)/(2*params.fsFilterVar));
    sMsFinalStateVarScaled[iq]  = sMsFinalStateVar[iq]
                                    *std::pow(
                                        exp(-1*std::pow(iq, 2)/(2*params.fsFilterVar)), 2);
  }
    plt.print1d(sMsFinalStateVar, "testVar_Final");
  for (int ir=0; ir<params.maxRbins; ir++) {
    pairCorrFinalStateVar[ir]  /= std::pow(params.NfinalPoints, 2); // Power of two for SEM
    pairCorrFinalStateSEM[ir]  = std::sqrt(pairCorrFinalStateVar[ir]);
  }

  std::string finalStateDataFileName = 
      "data-" + runName 
      + "_diffFinalState["
      + to_string(shape[1]) + "].dat";
  std::string finalStateDataSMSfileName = 
      "data-" + runName 
      + "_sMsFinalState["
      + to_string(shape[1]) + "].dat";
  std::string finalStateDataSMSScaledfileName = 
      "data-" + runName 
      + "_sMsFinalStateScaled["
      + to_string(shape[1]) + "].dat";
  std::string finalStateDataSMSSEMFileName = 
      "data-" + runName 
      + "_sMsFinalStateSEM["
      + to_string(shape[1]) + "].dat";
  std::string finalStateDataSMSSEMScaledFileName = 
      "data-" + runName 
      + "_sMsFinalStateSEMScaled["
      + to_string(shape[1]) + "].dat";
  std::string finalStateDataPairCorrSEMfileName = 
      "data-" + runName 
      + "_pairCorrFinalStateSEM["
      + to_string(shape[1]) + "].dat";

  save::saveDat<double>(finalState, 
      params.mergeScansOutputDir + finalStateDataFileName);
  save::saveDat<double>(sMsFinalState, 
      params.mergeScansOutputDir + finalStateDataSMSfileName);
  save::saveDat<double>(sMsFinalStateScaled, 
      params.mergeScansOutputDir + finalStateDataSMSScaledfileName);
  save::saveDat<double>(sMsFinalStateSEM, 
      params.mergeScansOutputDir + finalStateDataSMSSEMFileName);
  save::saveDat<double>(sMsFinalStateSEMScaled, 
      params.mergeScansOutputDir + finalStateDataSMSSEMScaledFileName);
  save::saveDat<double>(pairCorrFinalStateSEM, 
      params.mergeScansOutputDir + finalStateDataPairCorrSEMfileName);


  ///  Import Simulations  ///
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


  // Simulated final states
  int fillQbegin = (int)(shape[1]*params.fsQfitBegin/params.maxQazm);
  int fillQend   = (int)(shape[1]*params.fsQfitEnd/params.maxQazm);
  int fillRbegin = (int)(params.maxRbins*params.fsRfitBegin/params.maxR);
  int fillRend   = (int)(params.maxRbins*params.fsRfitEnd/params.maxR);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Xqs;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Xrs;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Xfit;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Yfit;
  Eigen::Matrix<double, Eigen::Dynamic, 1> W;
  Eigen::MatrixXd fit;
  Eigen::MatrixXd tmDpFit;
  Xqs.resize(fillQend-fillQbegin, NfitParams);
  Xrs.resize(fillRend-fillRbegin, NfitParams);
  std::vector<double> sMsFitCoeff(2,0), rawFitCoeff(2,0), pcFitCoeff(2,0);
  std::vector<double> pairCorrFinalState(params.maxRbins);
  std::vector<double> sMsFinalState_scaled(shape[1]);
  std::vector<double> finalState_scaled(shape[1]);
  std::vector<double> pairCorrFinalState_scaled(params.maxRbins);
  std::vector< std::vector<double> > XqsOut(NfitParams);
  std::vector< std::vector<double> > XrsOut(NfitParams);
  std::vector< std::vector<double> > simPairCorrFinalStates(params.finalStates.size());
  std::vector< std::vector<double> > simFinalDiffractions(params.finalStates.size());
  std::vector< std::vector<double> > simFinalSMSs(params.finalStates.size());

  if (params.fsFitOffset) {
    XqsOut[0].resize(fillQend-fillQbegin, 1);
    XrsOut[0].resize(fillRend-fillRbegin, 1);
    for (int i=0; i<fillQend - fillQbegin; i++) {
      Xqs(i,0) = 1;
    }
    for (int i=0; i<fillRend - fillRbegin; i++) {
      Xrs(i,0) = 1;
    }
  }
  for (int ifs=0; ifs<(int)params.finalStates.size(); ifs++) {
    simFinalDiffractions[ifs].resize(shape[1], 0);
    simFinalSMSs[ifs].resize(shape[1], 0);
    simPairCorrFinalStates[ifs].resize(params.maxRbins, 0);
    XqsOut[ifs+offShift].resize(fillQend-fillQbegin, 0);
    XrsOut[ifs+offShift].resize(fillRend-fillRbegin, 0);
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
          

    /////  Fitting Data  /////

    Xfit.resize(fillQend-fillQbegin, 1+offShift);
    Yfit.resize(fillQend-fillQbegin, 1);
    W.resize(fillQend-fillQbegin, 1);
    std::vector<double> YfitOut(params.NradAzmBins, 0);
    std::vector<double> YfitSEMout(params.NradAzmBins, 0);
    for (int iq=0; iq<params.NradAzmBins; iq++) {
      if ((iq >= fillQbegin) && (iq < fillQend)) {
        mInd = iq - fillQbegin;
        Xfit(mInd, 0) = 1;
        Xfit(mInd, offShift) = simSMSfinalState[iq];
        Xqs(mInd, ifs+offShift)     = simSMSfinalState[iq];
        XqsOut[ifs+offShift][mInd]  = simSMSfinalState[iq];
        if (params.fsFilterSMS) {
          Yfit(mInd, 0) = sMsFinalStateScaled[iq];
          W(mInd, 0)    = 1./sMsFinalStateVarScaled[iq];
          Xfit(mInd, offShift)        *= exp(-1*std::pow(iq, 2)/(2*params.fsFilterVar));
          Xqs(mInd, ifs+offShift)     *= exp(-1*std::pow(iq, 2)/(2*params.fsFilterVar));
          XqsOut[ifs+offShift][mInd]  *= exp(-1*std::pow(iq, 2)/(2*params.fsFilterVar));
        }
        else {
          Yfit(mInd, 0) = sMsFinalState[iq];
          W(mInd, 0)    = 1./sMsFinalStateVar[iq];
        }
      }

      if (params.fsFilterSMS) {
        YfitOut[iq]     = sMsFinalStateScaled[iq];
        YfitSEMout[iq]  = sMsFinalStateSEMScaled[iq];
      }
      else {
        YfitOut[iq]     = sMsFinalState[iq];
        YfitSEMout[iq]  = sMsFinalStateSEM[iq];
      }
    }

    fit = tools::normalEquation(Xfit, Yfit, W);

    if (params.fsFitOffset) {
      std::fill(sMsFinalState_scaled.begin(), sMsFinalState_scaled.end(), fit(0));
    }
    else {
      std::fill(sMsFinalState_scaled.begin(), sMsFinalState_scaled.end(), 0);
    }
    for (int iq=0; iq<shape[1]; iq++) {
      if (params.fsFilterSMS) {
        sMsFinalState_scaled[iq] += simSMSfinalState[iq]*fit(offShift)
            *exp(-1*std::pow(iq, 2)/(2*params.fsFilterVar));
      }
      else {
        sMsFinalState_scaled[iq] += simSMSfinalState[iq]*fit(offShift);
      }
    }

    save::saveDat<double>(sMsFinalState_scaled, 
        "./results/sim-" + fsName 
          + "_sMsFinalState_scaled[" 
          + to_string(shape[1]) + "].dat");

    save::saveDat<double>(YfitOut,
        params.mergeScansOutputDir 
          + "/data-" + runName
          + "_sMsFinalStateFittedTo["
          + to_string(shape[1]) + "].dat");

    save::saveDat<double>(YfitSEMout,
        params.mergeScansOutputDir 
          + "/data-" + runName
          + "_sMsFinalStateSEMFittedTo["
          + to_string(shape[1]) + "].dat");


  
    /*
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
          */


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


    // Fitting 
    Xfit.resize(fillRend-fillRbegin, 1+offShift);
    Yfit.resize(fillRend-fillRbegin, 1);
    W.resize(fillRend-fillRbegin, 1);
    for (int ir=fillRbegin; ir<fillRend; ir++) {
      mInd = ir - fillRbegin;
      Xfit(mInd, 0) = 1;
      Xfit(mInd, offShift) = simPairCorrFinalStates[ifs][ir];
      Yfit(mInd, 0) = pairCorrFinalState[ir];
      W(mInd, 0)    = 1/pairCorrFinalStateVar[ir];
      Xrs(mInd, ifs+offShift)     = simPairCorrFinalStates[ifs][ir];
      XrsOut[ifs+offShift][mInd]  = simPairCorrFinalStates[ifs][ir];
    }

    fit = tools::normalEquation(Xfit, Yfit, W);
    if (params.fsFitOffset) {
      std::fill(
          pairCorrFinalState_scaled.begin(), 
          pairCorrFinalState_scaled.end(), 
          fit(0));
    }
    else {
      std::fill(
          pairCorrFinalState_scaled.begin(), 
          pairCorrFinalState_scaled.end(), 
          0);
    }
    for (int ir=0; ir<params.maxRbins; ir++) {
      pairCorrFinalState_scaled[ir] += pairCorrFinalState[ir]*fit(offShift);
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



  ////////////////////////////////////////////////////////////
  /////  Fitting Linear Combination of All Final States  /////
  ////////////////////////////////////////////////////////////

  std::string codeDir = "/reg/neh/home5/khegazy/baseTools/UEDanalysis/timeDepStudies/";
  std::string XrsFileName = "./results/fitRvaluesX_Bins["
      + to_string(XrsOut.size()) + ","
      + to_string(XrsOut[0].size()) + "].dat";
  std::string XqsFileName = "./results/fitSMSvaluesX_Bins["
      + to_string(XqsOut.size()) + ","
      + to_string(XqsOut[0].size()) + "].dat";
  std::string coeffsFileName = "./results/fitCoefficients_Bins["
      + to_string(XqsOut.size()) + "].dat";
  std::string coeffErrsFileName = "./results/fitCoefficientErrors_Bins["
      + to_string(XqsOut.size()) + "].dat";

  save::saveDat<double>(XrsOut, XrsFileName);
  save::saveDat<double>(XqsOut, XqsFileName);

  /////  Final State  /////
  std::vector<double> Wout;
  std::vector<double> Yout;
  std::vector<double> fitCoeffs(NfitParams);
  std::vector<double> fitCoeffErrs(NfitParams);
  std::string WqFileName = "./results/fitQvaluesW_Bins["
      + to_string(fillQend-fillQbegin) + "].dat";
  std::string WrFileName = "./results/fitRvaluesW_Bins["
      + to_string(fillRend-fillRbegin) + "].dat";
  std::string YqFileName = "./results/fitQvaluesY_Bins["
      + to_string(fillQend-fillQbegin) + "].dat";
  std::string YrFileName = "./results/fitRvaluesY_Bins["
      + to_string(fillRend-fillRbegin) + "].dat";

  ///  SMS  ///
  Yout.resize(fillQend-fillQbegin);
  Wout.resize(fillQend-fillQbegin);
  for (int iq=fillQbegin; iq<fillQend; iq++) {
    mInd = iq - fillQbegin;
    if (params.fsFilterSMS) {
      Yout[mInd]  = sMsFinalStateScaled[iq];
      Wout[mInd]  = sMsFinalStateSEMScaled[iq];
    }
    else {
      Yout[mInd]  = sMsFinalState[iq];
      Wout[mInd]  = sMsFinalStateSEM[iq];
    }
  }

  save::saveDat<double>(Wout, WqFileName);
  save::saveDat<double>(Yout, YqFileName);

  std::system(("python " + codeDir 
      + "fitLinCurveFit.py --X " + XqsFileName
      + " --Y " + YqFileName
      + " --W " + WqFileName).c_str());

  save::importDat<double>(fitCoeffs, coeffsFileName);
  save::importDat<double>(fitCoeffErrs, coeffErrsFileName);

  std::vector<double> combinedQ(shape[1]);
  for (int iq=0; iq<shape[1]; iq++) {
    for (int i=0; i<NfitParams; i++) {
      combinedQ[iq] += XqsOut[i][iq]*fitCoeffs[i];
    }
  }

  save::saveDat<double>(combinedQ, 
      "./results/sim-" + runName
      + "_sMsFinalState_scaledLinComb_Bins["
      + to_string(combinedQ.size()) + "].dat");


  ///  Pair Correlation  ///
  Yout.resize(fillRend-fillRbegin);
  Wout.resize(fillRend-fillRbegin);
  for (int ir=fillRbegin; ir<fillRend; ir++) {
    mInd = ir - fillRbegin;
    Yout[mInd]  = pairCorrFinalState[ir];
    Wout[mInd]  = pairCorrFinalStateVar[ir];
  }

  save::saveDat<double>(Wout, WrFileName);
  save::saveDat<double>(Yout, YrFileName);

  std::system(("python " + codeDir 
      + "fitLinCurveFit.py --X " + XrsFileName
      + " --Y " + YrFileName
      + " --W " + WrFileName).c_str());

  save::importDat<double>(fitCoeffs, coeffsFileName);
  save::importDat<double>(fitCoeffErrs, coeffErrsFileName);

  std::vector<double> combinedR(params.maxRbins);
  for (int ir=0; ir<params.maxRbins; ir++) {
    for (int i=0; i<NfitParams; i++) {
      combinedR[ir] += XrsOut[i][ir]*fitCoeffs[i];
    }
  }

  save::saveDat<double>(combinedR, 
      "./results/sim-" + runName
      + "_pairCorrFinalState_scaledLinComb_Bins["
      + to_string(combinedR.size()) + "].dat");


  /////  Time Dependent  /////
  Yout.resize(fillQend-fillQbegin);
  Wout.resize(fillQend-fillQbegin);
  std::vector< std::vector<double> > combinedQs(tmDpSMS.size());
  std::vector< std::vector<double> > tmDpQfits(NfitParams);
  std::vector< std::vector<double> > tmDpQfitErrs(NfitParams);

  ///  sMs  ///
  for (int i=0; i<(int)tmDpQfits.size(); i++) {
    tmDpQfits[i].resize(tmDpSMS.size());
    tmDpQfitErrs[i].resize(tmDpSMS.size());
  }
  for (int it=0; it<(int)tmDpSMS.size(); it++) {
    cout<<"time: "<<it<<endl;
    combinedQs[it].resize(params.NradAzmBins, 0);
    for (int iq=fillQbegin; iq<fillQend; iq++) {
      mInd = iq - fillQbegin;
      Yout[mInd]  = tmDpSMS[it][iq];
      Wout[mInd]  = tmDpSMSSEM[it][iq];
      if (params.fsFilterSMS) {
        Yout[mInd] *= exp(-1*std::pow(iq, 2)/(2*params.fsFilterVar));
        Wout[mInd] *= exp(-1*std::pow(iq, 2)/(2*params.fsFilterVar));
      }
    }

    save::saveDat<double>(Wout, WqFileName);
    save::saveDat<double>(Yout, YqFileName);

    cout<<"Starting python fit"<<endl;
    std::system(("python " + codeDir 
        + "fitLinCurveFit.py --X " + XqsFileName
        + " --Y " + YqFileName
        + " --W " + WqFileName//).c_str());
        + " --plot testFit_" + to_string(it)).c_str());

    save::importDat<double>(fitCoeffs, coeffsFileName);
    save::importDat<double>(fitCoeffErrs, coeffErrsFileName);

    for (int ifs=0; ifs<(int)tmDpQfits.size(); ifs++) {
      tmDpQfits[ifs][it]    = fitCoeffs[ifs];
      tmDpQfitErrs[ifs][it] = fitCoeffErrs[ifs];
    }

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      for (int i=0; i<NfitParams; i++) {
        combinedQs[it][iq] += XqsOut[i][iq]*fitCoeffs[i];
      }
    }
  }

  if (params.fsFitOffset) {
    save::saveDat<double>(tmDpQfits[0], 
        "./results/sim-" + runName
          + "-offset_sMsFitCoeffsLinComb_Bins[" 
          + to_string(tmDpSMS.size()) + "].dat");
    save::saveDat<double>(tmDpQfitErrs[0], 
        "./results/sim-" + runName
          + "-offset_sMsFitCoeffErrorsLinComb_Bins[" 
          + to_string(tmDpSMS.size()) + "].dat");
  }
  for (int ifs=0; ifs<(int)params.finalStates.size(); ifs++) {
    save::saveDat<double>(tmDpQfits[ifs+offShift], 
        "./results/sim-" + runName + "-" + params.finalStates[ifs]
          + "_sMsFitCoeffsLinComb_Bins[" 
          + to_string(tmDpSMS.size()) + "].dat");
    save::saveDat<double>(tmDpQfitErrs[ifs+offShift], 
        "./results/sim-" + runName + "-" + params.finalStates[ifs]
          + "_sMsFitCoeffErrorsLinComb_Bins[" 
          + to_string(tmDpSMS.size()) + "].dat");
  }

  save::saveDat<double>(combinedQs, 
      "./results/sim-" + runName
      + "_sMsFinalStateFitLinComb_Bins["
      + to_string(combinedQs.size()) + ","
      + to_string(combinedQs[0].size()) + "].dat");


cout<<"starting paircorr"<<endl;

  ///  Pair Correlation  ///
  Yout.resize(fillRend-fillRbegin);
  Wout.resize(fillRend-fillRbegin);
  std::vector< std::vector<double> > combinedRs(tmDpPairCorr.size());
  std::vector< std::vector<double> > tmDpRfits(NfitParams);
  std::vector< std::vector<double> > tmDpRfitErrs(NfitParams);
  for (int i=0; i<(int)tmDpRfits.size(); i++) {
    tmDpRfits[i].resize(tmDpPairCorr.size());
    tmDpRfitErrs[i].resize(tmDpPairCorr.size());
  }
  for (int it=0; it<(int)tmDpPairCorr.size(); it++) {
    combinedRs[it].resize(params.NradAzmBins, 0);
    for (int ir=fillRbegin; ir<fillRend; ir++) {
      mInd = ir - fillRbegin;
      Yout[mInd] = tmDpPairCorr[it][ir];
      Wout[mInd] = tmDpPairCorrSEM[it][ir];
    }

    save::saveDat<double>(Wout, WrFileName);
    save::saveDat<double>(Yout, YrFileName);

    std::system(("python " + codeDir 
        + "fitLinCurveFit.py --X " + XrsFileName
        + " --Y " + YrFileName
        + " --W " + WrFileName).c_str());

    save::importDat<double>(fitCoeffs, coeffsFileName);
    save::importDat<double>(fitCoeffErrs, coeffErrsFileName);

    for (int ifs=0; ifs<(int)tmDpRfits.size(); ifs++) {
      tmDpRfits[ifs][it]    = fitCoeffs[ifs];
      tmDpRfitErrs[ifs][it] = fitCoeffErrs[ifs];
    }


    for (int ir=0; ir<params.maxRbins; ir++) {
      for (int i=0; i<NfitParams; i++) {
        combinedRs[it][ir] += XrsOut[i][ir]*fitCoeffs[i];
      }
    }
  }

  if (params.fsFitOffset) {
    save::saveDat<double>(tmDpRfits[0], 
        "./results/sim-" + runName
          + "-offset_pairCorrFitCoeffsLinComb_Bins[" 
          + to_string(tmDpPairCorr.size()) + "].dat");
    save::saveDat<double>(tmDpRfitErrs[0], 
        "./results/sim-" + runName
          + "-offset_pairCorrFitCoeffErrorsLinComb_Bins[" 
          + to_string(tmDpPairCorr.size()) + "].dat");
  }

  for (int ifs=0; ifs<(int)params.finalStates.size(); ifs++) {
    save::saveDat<double>(tmDpRfits[ifs+offShift], 
        "./results/sim-" + runName + "-" +params.finalStates[ifs] 
          + "_pairCorrFitCoeffsLinComb_Bins[" 
          + to_string(tmDpSMS.size()) + "].dat");
    save::saveDat<double>(tmDpRfitErrs[ifs+offShift], 
        "./results/sim-" + runName + "-" +params.finalStates[ifs] 
          + "_pairCorrFitCoeffErrorsLinComb_Bins[" 
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
