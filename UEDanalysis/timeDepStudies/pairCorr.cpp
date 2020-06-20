#include "../analysis.h"
#include "/reg/neh/home/khegazy/baseTools/tools/parameters.h"

using namespace std;


int main(int argc, char* argv[]) {

  if (argc<2) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe 'runName' !!!"<<endl;
    cerr<<"         Can also run ./analysis.exe 'fileList.txt' 'treeName'"<<endl;
    exit(0);
  }


  cout.setf(ios::scientific);
  
  //////////////////////////////////
  /////  Setting up variables  /////
  //////////////////////////////////

  std::string runName(argv[1]);
  parameterClass params(runName);

  bool givenFileName      = false;
  bool subtractReference  = false;
  bool fitLowQ            = false;
  bool removeSkippedBins  = true;
  bool fitPC              = false;
  std::string fileName    = "NULL";
  std::string inputDir    = "NULL";
  std::string outPrefix   = "data-" + runName;
  std::string outSuffix   = "";
  std::string lowQtheory  = "NULL";
  std::string saveLowQtheory  = "NULL";
  std::string outputDir   = "./results/";
  int fitTstep = -1;

  /////  Importing variables from command line  /////
  if (params.verbose)
    std::cout << "Importing command line variables.\n";

  for (int iarg=2; iarg<argc; iarg+=2) {
    if (strcmp(argv[iarg], "-MolName") == 0) {
      string str(argv[iarg+1]); 
      params.molName = str;
    }
    else if (strcmp(argv[iarg], "-Odir") == 0) {
      string str(argv[iarg+1]); 
      params.mergeScansOutputDir = str;
    }
    else if (strcmp(argv[iarg], "-Idir") == 0) {
      string str(argv[iarg+1]);
      inputDir = str;
    }
    else if (strcmp(argv[iarg], "-Fname") == 0) {
      givenFileName = true;
      string str(argv[iarg+1]);
      fileName = str;
    }
    else if (strcmp(argv[iarg], "-fitTstep") == 0) {
      fitTstep = atoi(argv[iarg+1]);
    }
    else if (strcmp(argv[iarg], "-Osuf") == 0) {
      string str(argv[iarg+1]);
      outSuffix = str;
    }
    else if (strcmp(argv[iarg], "-OutDir") == 0) {
      cout<<"HCASDFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"<<endl;
      string str(argv[iarg+1]);
      outputDir = str;
    }
    else if (strcmp(argv[iarg], "-lowQtheory") == 0) {
      string str(argv[iarg+1]);
      params.fillLowQtheory = true;
      params.fillLowQfile   = str;
    }
    else if (strcmp(argv[iarg], "-saveLowQtheory") == 0) {
      string str(argv[iarg+1]);
      saveLowQtheory = str;
    }
    else if (strcmp(argv[iarg], "-FitPC") == 0) {
      string str(argv[iarg+1]);
      if (str.compare("true") == 0) fitPC = true;
      else {
        std::cerr << "ERROR: Do not understand FitPC option " << str << endl;
        exit(0);
      }
    }
 
    else if (strcmp(argv[iarg], "-FillQ") == 0) {
      string str(argv[iarg+1]);
      if (str.compare("false") == 0) params.fillLowQtheory = false;
      else if (str.compare("zeros") == 0) {
        params.fillLowQzeros = true;
        params.fillLowQtheory = false;
      }
      else if (str.compare("fit") == 0) {
        fitLowQ = true;
        params.fillLowQtheory = false;
      }
      else if (str.compare("theory") == 0) params.fillLowQtheory = true;
      else {
        std::cerr << "ERROR: Do not understand FillQ option " << str << endl;
        exit(0);
      }
    }
    else if (strcmp(argv[iarg], "-SubR") == 0) {
      string str(argv[iarg+1]);
      if (str.compare("true") == 0) subtractReference = true;
      else if (str.compare("false") == 0) ;
      else {
        std::cerr << "ERROR: Do not understand SubR option " << str << endl;
        exit(0);
      }
    }

    else {
      cerr<<"ERROR!!! Option "<<argv[iarg]<<" does not exist!"<<endl;
      exit(0);
    }
  }

  if ((inputDir.find("sim") != string::npos) 
      || (inputDir.find("Sim") != string::npos)) {
    outPrefix = "sim-" + params.molName;
    removeSkippedBins = false;
  }

  /////  Plotting variables  /////
  std::vector<PLOToptions> opts(6);
  std::vector<std::string> vals(6);
  std::vector<PLOToptions> oppts(3);
  std::vector<std::string> vaals(3);
  oppts[0] = yLabel;   vaals[0] = "Time [ps]";
  oppts[1] = xLabel;   vaals[1] = "Scattering Q [arb units]";
  oppts[2] = draw;     vaals[2] = "COLZ";
  opts[0] = yLabel;   vals[0] = "Time [ps]";
  opts[1] = xLabel;   vals[1] = "Scattering Q [arb units]";
  opts[2] = draw;     vals[2] = "COLZ";
  opts[3] = minimum;  vals[3] = "-2e-2";
  opts[4] = maximum;  vals[4] = "2e-2";
  opts[5] = xSpan;    


  ////////////////////////////
  /////  Importing data  /////
  ////////////////////////////

  if (params.verbose)
    std::cout << "Importing data.\n";

  std::vector<int> shape;
 
  // Importing Data
  if (!givenFileName) {
    shape = save::getShape( 
        params.mergeScansOutputDir,
        "data-" + runName + "-sMsAzmAvgDiff[");

    assert(shape[1] == params.NradAzmBins);

    std::string fileType = "sMsAzmAvgDiff";
    if (params.useFilledHole) {
      fileType += "_filledHole";
    }
    fileName = params.mergeScansOutputDir 
        + "data-" + runName + "-" + fileType +"[" 
        + to_string(shape[0])
        + "," + to_string(shape[1]) + "].dat";
    cout<<"FILENAME "<<fileName<<endl;
  }
  else {
    assert(strcmp(inputDir.c_str(), "NULL") != 0);

    shape = save::getShape(inputDir, fileName);
    cout<<"got shape "<<shape[0]<<endl;
    if (shape.size() == 1) {
      shape.insert(shape.begin(), 1);
    }

    params.NradAzmBins = shape[1];
    
    fileName = inputDir + "/" + fileName;
  }

  std::vector< std::vector<double> > tmdDiff;
  tmdDiff.resize(shape[0]);
  for (auto& v : tmdDiff) {
    v.resize(shape[1], 0);
  }

  if (params.verbose)
    std::cout << "Retriving data from " << fileName << endl; 

  std::cout << "Retriving data from " << fileName << endl; 
  save::importDat<double>(tmdDiff, fileName);
   
  double pCorrQcut = 8;
  int NqBins = (int)(pCorrQcut*params.NradAzmBins/params.maxQazm);
  for (uint tm=0; tm<tmdDiff.size(); tm++) {
    tmdDiff[tm].resize(NqBins);
  }

  if (params.pltVerbose) {
    std::vector<PLOToptions> pOpts(2);
    std::vector<std::string> pVals(2);
    pOpts[0] = minimum;  pVals[0] = "-2e-2";
    pOpts[1] = maximum;  pVals[1] = "2e-2";
    plt.printRC(tmdDiff, "./plots/importedTDdiff", pOpts, pVals);
  }


  ///////////////////////////////////
  /////  Importing Simulations  /////
  ///////////////////////////////////

  std::string sMsSimName;
  std::vector<double> sMsGroundState(shape[1]), sMsSimDer(shape[1]);
  if (subtractReference || params.simHotFinalState) {
    // Importing Simulation
    sMsSimName = 
        params.simOutputDir + "nitrobenzene_sMsPatternLineOut_Qmax-"
        + to_string(params.maxQazm) + "_Ieb-"
        + to_string(params.Iebeam) + "_scrnD-"
        + to_string(params.screenDist) + "_elE-"
        + to_string(params.elEnergy) + "_Bins["
        + to_string(shape[1]) + "].dat";
    
    save::importDat<double>(sMsGroundState, sMsSimName);
  }
  if (subtractReference) {
    for (int itm=0; itm<shape[0]; itm++) {
      for (int iq=0; iq<(int)tmdDiff[itm].size(); iq++) {
        tmdDiff[itm][iq] -= sMsGroundState[iq];
      }
    }
  }
 
  std::vector<double> lowQsim(shape[1]);
  if (params.fillLowQtheory) { 
    save::importDat<double>(lowQsim, params.fillLowQfile);
    cout<<"shouldn't be here"<<endl;
    if (params.verbose)
      std::cout << "Importing simulation from " << params.fillLowQfile << endl;
  }


  // Saving lowQ theory
  if (saveLowQtheory.compare("NULL") != 0) {
    std::string fName = outputDir + "/theoryLowQ_" 
                          + params.molName
                          + outSuffix;
    if (tmdDiff.size() == 1) {
      save::saveDat(tmdDiff[0], 
          fName + "[" + to_string(tmdDiff[0].size()) + "].dat");
    }
    else {
      save::saveDat(tmdDiff, fName 
          + "[" + to_string(tmdDiff[0].size())
          + to_string(tmdDiff[0].size()) + "].dat");
    }
  }


  // Importing low Q fit coefficients
  if (params.verbose)
    std::cout << "Importing low Q coefficients\n";

  std::string lowQcoeffName = outputDir + "/lowQcoefficients_"
      + runName + "["
      + to_string(shape[0]) + ",4].dat";
  std::string lowQwName = outputDir + "/lowQfrequencies_"
      + runName + "["
      + to_string(shape[0]) + ",4].dat";
  bool lowQcoeffExists = tools::fileExists(lowQcoeffName);
  std::vector< std::vector<double> > lowQcoeff, lowQw;
  if (lowQcoeffExists && fitLowQ) {
    lowQcoeff.resize(shape[0]);
    lowQw.resize(shape[0]);
    for (int i=0; i<shape[0]; i++) {
      lowQcoeff[i].resize(4);
      lowQw[i].resize(4);
    }

    save::importDat<double>(lowQcoeff, lowQcoeffName);
    save::importDat<double>(lowQw, lowQwName);
  }
  else {
    lowQcoeff.resize(shape[0]);
    if (fitTstep != -1) {
      if (fitTstep >= shape[0]) {
        exit(1);
      }  
    }

    // Combine results of fitting individual images to single file
    std::string coeffFileName = outputDir + "/lowQfitCoeff_" + runName + "_time-";
    std::string wFileName = outputDir + "/lowQfitW_" + runName + "_time-";
    if (tools::fileExists(coeffFileName + "0.dat") && (fitTstep == -1)) {
      lowQw.resize(shape[0]);
      std::vector<double> coeffs(4), ws(4);
      for (int itm=0; itm<shape[0]; itm++) {
        save::importDat<double>(coeffs, coeffFileName + to_string(itm) + ".dat");
        save::importDat<double>(ws, wFileName + to_string(itm) + ".dat");
        lowQcoeff[itm].resize(4);
        lowQw[itm].resize(4);
        for (uint i=0; i<coeffs.size(); i++) {
          lowQcoeff[itm][i] = coeffs[i];
          lowQw[itm][i] = ws[i];
        }
      }

      save::saveDat<double>(lowQcoeff, lowQcoeffName);
      save::saveDat<double>(lowQw, lowQwName);
      lowQcoeffExists = true;
    }
  }

  std::vector< std::vector<double> > fitTmdDiff(tmdDiff.size());
  if (params.fillLowQfitTheory) {
    std::string simSMSrefName =
        params.simOutputDir
        + params.radicalNames[params.initialState]
        + "_sMsPatternLineOut_Qmax-"
        + to_string(params.maxQazm) + "_Ieb-"
        + to_string(params.Iebeam) + "_scrnD-"
        + to_string(params.screenDist) + "_elE-"
        + to_string(params.elEnergy) + "_Bins["
        + to_string(params.NradAzmBins) + "].dat";

    std::vector<double> sMsReference(params.NradAzmBins);
    save::importDat<double>(sMsReference, simSMSrefName);
    
    std::vector< std::vector<double> > simFinalSMSs(params.finalStates.size());
    std::vector< std::vector<double> > tmdFitCoeffs(params.finalStates.size());
    for (int ifs=0; ifs<(int)params.finalStates.size(); ifs++) {
      std::string fsName = params.finalStates[ifs];
 
      std::string simSMSfinalName =
          params.simOutputDir + fsName
          + "_sMsPatternLineOut_Qmax-"
          + to_string(params.maxQazm) + "_Ieb-"
          + to_string(params.Iebeam) + "_scrnD-"
          + to_string(params.screenDist) + "_elE-"
          + to_string(params.elEnergy) + "_Bins["
          + to_string(params.NradAzmBins) + "].dat";
 
      simFinalSMSs[ifs].resize(params.NradAzmBins);
      save::importDat<double>(simFinalSMSs[ifs], simSMSfinalName);

      for (int iq=0; iq<params.NradAzmBins; iq++) {
        simFinalSMSs[ifs][iq] -= sMsReference[iq];
      }

      std::string fitCoeffFileName = 
        outputDir + "/sim-" + runName + "-" + params.finalStates[ifs]
        + "_sMsFitCoeffsLinComb_Bins["
        + to_string(tmdDiff.size()) + "].dat";

      tmdFitCoeffs[ifs].resize(tmdDiff.size());
      save::importDat<double>(tmdFitCoeffs[ifs], fitCoeffFileName);
    }

    for (uint i=0; i<tmdDiff.size(); i++) {
      fitTmdDiff[i].resize(params.NradAzmBins, 0);
    }
    for (int itm=0; itm<(int)tmdDiff.size(); itm++) {
      for (int ifs=0; ifs<(int)params.finalStates.size(); ifs++) {
        for (int iq=0; iq<params.NradAzmBins; iq++) {
          fitTmdDiff[itm][iq] += simFinalSMSs[ifs][iq]*tmdFitCoeffs[ifs][itm];
        }
      }
      plt.print1d(fitTmdDiff[itm], "testing-" + to_string(itm));
    }
  }

    





  ///////////////////////////////////////
  /////  Pair correlation function  /////
  ///////////////////////////////////////
  
  if (params.verbose) {
    std::cout << "Begin calculating pair correlation functions\n";
  }

  bool mirrorDiffPattern = false;
  int inpFFTsize;
  if (mirrorDiffPattern) inpFFTsize = params.NradAzmBins*2 + 2*params.NautCpadding - 1;
  else inpFFTsize = params.NradAzmBins + params.NautCpadding;
  int outFFTsize = inpFFTsize/2 + 1;
  int outSize = params.rMaxAzmRat*outFFTsize;

  double* qSpace = (double*) fftw_malloc(sizeof(double)*inpFFTsize);
  fftw_complex* rSpace =
        (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*outFFTsize);
  fftw_plan fftB;
  fftB = fftw_plan_dft_r2c_1d(inpFFTsize, qSpace, rSpace, FFTW_MEASURE);

  int FFTcent = inpFFTsize/2;
  std::vector<double> fftInp(inpFFTsize);
  std::vector<double> smoothed(NqBins, 0);
  std::vector<double> filter(NqBins, 0);
  std::vector< std::vector<double> > results;
  std::vector< std::vector<double> > pairCorrEven(shape[0]);
  std::vector< std::vector<double> > pairCorrOdd(shape[0]);
  for (int i=0; i<shape[0]; i++) {
    pairCorrEven[i].resize(outSize);
    pairCorrOdd[i].resize(outSize);
  }


  double filterScale = 1;
  int maxFitInd = 100; 
  Eigen::MatrixXd fit, bestFit;
  Eigen::Matrix<double, Eigen::Dynamic, 4> X;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Y;
  X.resize(maxFitInd - params.NbinsSkip, 4);
  Y.resize(maxFitInd - params.NbinsSkip, 1);

  std::vector<double> lowPassFilter(NqBins);
  if (params.pCorrButterFilter) {
    std::string filterName =
        "/reg/neh/home/khegazy/analysis/filters/" + params.pCorrFilterType
        + "Filter_Bins-" + to_string(NqBins)
        + "_order-" + to_string(params.pCorrFilterOrder)
        + "_WnHigh-"+ to_string(params.pCorrWnHigh) + ".dat";
    if (!tools::fileExists(filterName)) {
      cout << "\n\nINFO: Making new filter\n";
      system(("python /reg/neh/home/khegazy/analysis/filters/makeFilters.py --Nbins "
            + to_string(NqBins)
            + " --Ftype " + params.pCorrFilterType
            + " --Order " + to_string(params.pCorrFilterOrder)
            + " --WnHigh " + to_string(params.pCorrWnHigh)).c_str());
    }
    save::importDat<double>(lowPassFilter, filterName);
  }


  for (int itm=0; itm<shape[0]; itm++) {

    if (fitTstep != -1) {
      if (itm != fitTstep) {
        continue;
      }
    }

    int iq = 0;
    if (removeSkippedBins && !params.useFilledHole) {
      for ( ; iq<params.NbinsSkip; iq++) 
        tmdDiff[itm][iq] = 0;
    }
    iq = params.NbinsSkip;

    if (!lowQcoeffExists && fitLowQ) {
      double bestLoss = 1e100;
      double rScaling = 0.1;
      double fitVal   = 0;
      int Nsteps      = 50;
      std::vector<double> wCur(4);
      std::vector<double> w(4);
      
      for (iq=params.NbinsSkip; iq<maxFitInd; iq++) {
        Y(iq-params.NbinsSkip,0) = tmdDiff[itm][iq];
      }

      for (int i=0; i<Nsteps; i++) {
        wCur[0] = params.QperPix*(1 + rScaling*i);
        for (iq=params.NbinsSkip; iq<maxFitInd; iq++) {
          X(iq-params.NbinsSkip,0) = sin(wCur[0]*iq);
        }
        for (int j=i+1; j<Nsteps; j++) {
          if (i == j) continue;
          wCur[1] = params.QperPix*(1 + rScaling*j);
          for (iq=params.NbinsSkip; iq<maxFitInd; iq++) {
            X(iq-params.NbinsSkip,1) = sin(wCur[1]*iq);
          }
          for (int k=j+1; k<Nsteps; k++) {
            if (k == j) continue;
            wCur[2] = params.QperPix*(1 + rScaling*k);
            for (iq=params.NbinsSkip; iq<maxFitInd; iq++) {
              X(iq-params.NbinsSkip,2) = sin(wCur[2]*iq);
            }
            for (int l=k+1; l<Nsteps; l++) {
              if (l == k) continue;
              wCur[3] = params.QperPix*(1 + rScaling*l);
              for (iq=params.NbinsSkip; iq<maxFitInd; iq++) {
                X(iq-params.NbinsSkip,3) = sin(wCur[3]*iq);
              }

              fit = tools::normalEquation(X,Y);

              // Skip fits with positive contributions
              try {
                for (int iqq=0; iqq<params.NbinsSkip*0.75; iqq++) {
                  fitVal = 0;
                  for (int n=0; n<fit.rows(); n++) {
                    fitVal += fit(n)*sin(wCur[n]*iqq);
                  }
                  if (fitVal > 0) {
                    throw 0;
                  }
                }
              }
              catch (int e) {
                continue;
              }

              // Check if this is the best fit
              if (bestLoss > (X*fit - Y).squaredNorm()) {
                cout<<"fitting: "<<i<<"  "<<j<<"  "<<k<<"  "<<l<<"  "<<(X*fit - Y).squaredNorm()<<"   "<<fit<<endl;
                w = wCur;
                bestFit = fit;
                bestLoss = (X*fit - Y).squaredNorm();
              }
            }
          }
        }
      }
      lowQw.push_back(w);
      lowQcoeff[itm].resize(bestFit.rows());
      for (int i=0; i<bestFit.rows(); i++) {
        lowQcoeff[itm][i] = bestFit(i);
      }

      if (fitTstep != -1) {
        save::saveDat<double>(lowQw[0], outputDir + "/lowQfitW_" 
            + runName + "_time-" + to_string(itm) + ".dat");
        save::saveDat<double>(lowQcoeff[itm], outputDir + "/lowQfitCoeff_" 
            + runName + "_time-" + to_string(itm) + ".dat");
        exit(1);
      }
    }


    /////  Filling FFT Input  /////
    for (iq=0; iq<NqBins; iq++) {
      if (!params.useFilledHole) {
        if (params.fillLowQtheory) {
          if (iq < params.NbinsSkip) {
            tmdDiff[itm][iq] = lowQsim[iq];
                //*(tmdDiff[0][params.NbinsSkip]/lowQsim[params.NbinsSkip]);
          }
        }
        else if (params.fillLowQzeros) {
          if (iq < params.NbinsSkip) {
            tmdDiff[itm][iq] = 0;
          }
        }
        else if (fitLowQ) {
          if (iq < params.NbinsSkip) {
            for (uint i=0; i<lowQw[itm].size(); i++) {
              if (i == 0) {
                tmdDiff[itm][iq] = lowQcoeff[itm][i]*sin(lowQw[itm][i]*iq); 
                continue;
              }
              tmdDiff[itm][iq] += lowQcoeff[itm][i]*sin(lowQw[itm][i]*iq); 
            }
          }
        }
        else if (params.fillLowQsine) {
          if (iq < params.NbinsSkip) {
            tmdDiff[itm][iq] = std::sin(iq*PI/(2*50))
                                *tmdDiff[itm][params.NbinsSkip]
                                /std::sin(params.NbinsSkip*PI/(2*50));
          }
        }
        else if (params.fillLowQfitTheory) {
          if (iq < params.NbinsSkip) {
            tmdDiff[itm][iq] = fitTmdDiff[itm][iq];
          }
        }
      }


      //  Apply Filter
      if (params.pCorrGaussFilter) {
        filter[iq] = exp(-1*std::pow(iq, 2)
                            /(2*params.filterVar));
      }
      else if (params.pCorrButterFilter) {
        filter[iq] = lowPassFilter[iq];
      }

      if (!params.simHotFinalState 
          || (runName.compare("simulateReference") != 0)) {
        smoothed[iq] = tmdDiff[itm][iq]*filter[iq];
      }
      else if (!params.useFilledHole && (!params.fillLowQzeros 
          || (params.fillLowQzeros && (iq > params.NbinsSkip)))) {
        smoothed[iq] = (tmdDiff[itm][iq] + sMsGroundState[iq])
                        *exp(-1*std::pow(iq, 2)
                            /(2*params.hotFStdepVar));
        smoothed[iq] -= sMsGroundState[iq] 
                        *exp(-1*std::pow(iq, 2)
                            /(2*params.hotFSrefVar));
      }
      else {
        smoothed[iq] = 0;
      }
    }

    if (params.pltVerbose) {
      plt.print1d(tmdDiff[itm], "./plots/inpDiff_" + runName + "_" + to_string(itm));
      plt.print1d(smoothed, "./plots/smoothed_" + runName + "_" + to_string(itm));
      plt.print1d(filter, "./plots/filter");
    }


    /////  Apply FFT  /////
    if (mirrorDiffPattern) {
      fftInp[FFTcent] = smoothed[0];
      for (iq=1; iq<NqBins; iq++) {
        fftInp[FFTcent+iq] = smoothed[iq];
        fftInp[FFTcent-iq] = smoothed[iq];
      }
      results = tools::fft1dRtoC(fftInp, fftB, qSpace, rSpace, 
                    false, true, NULL, "fftInp_" + to_string(itm));
    }
    else {
      for (iq=0; iq<NqBins; iq++) {
        fftInp[iq] = smoothed[iq];
      }
      results = tools::fft1dRtoC(fftInp, fftB, qSpace, rSpace, 
                    false, false, NULL, "fftInp_" + to_string(itm));
    }

    if (params.pltVerbose)
      plt.print1d(fftInp, "./plots/fftFuncInp_" + to_string(itm));

    for (int ir=0; ir<outSize; ir++) {
      pairCorrEven[itm][ir] = std::pow(results[0][ir], 1);
      pairCorrOdd[itm][ir] = std::pow(results[1][ir], 1);
    }
  }

  if (params.pltVerbose) {
    vals[5] = "0," + to_string(params.rMaxAzmRat*outFFTsize*(2*PI/(params.QperPix*inpFFTsize)));
    if (pairCorrEven.size() == 1) {
      plt.print1d(pairCorrEven[0], "./plots/pairCorrEven");
      plt.print1d(pairCorrOdd[0], "./plots/pairCorrOdd");
    }
    else {
      plt.printRC(pairCorrEven, "./plots/pairCorrEven", opts, vals);
      plt.printRC(pairCorrOdd, "./plots/pairCorrOdd", opts, vals);
    }
  }


  ////////////////////
  /////  Saving  /////
  ////////////////////

  if (params.verbose)
    std::cout << "Saving.\n";

  if (!lowQcoeffExists && fitLowQ) {
    for (auto& v : lowQw[0]) cout<<"testing: " << v/params.QperPix<<endl;
    save::saveDat<double>(lowQcoeff, lowQcoeffName);
    save::saveDat<double>(lowQw, lowQwName);
  }

  if (shape[0] != 1) {
    save::saveDat<double>(pairCorrEven, 
        outputDir + "/" + outPrefix 
        + outSuffix + "_pairCorrEven["
        + to_string(pairCorrEven.size()) + ","
        + to_string(pairCorrEven[0].size()) + "].dat");
    save::saveDat<double>(pairCorrOdd, 
        outputDir + "/" + outPrefix   
        + outSuffix + "_pairCorrOdd["
        + to_string(pairCorrOdd.size()) + ","
        + to_string(pairCorrOdd[0].size()) + "].dat");
  }
  else {
    cout<<"HERE SAVING " << outputDir + "/" + outPrefix+ outSuffix + "_pairCorrEven["+ to_string(pairCorrEven.size()) + ","+ to_string(pairCorrEven[0].size()) + "].dat"<<endl;
    save::saveDat<double>(pairCorrEven[0], 
        outputDir + "/" + outPrefix 
        + outSuffix + "_pairCorrEven["
        + to_string(pairCorrEven[0].size()) + "].dat");
    save::saveDat<double>(pairCorrOdd[0], 
        outputDir + "/" + outPrefix 
        + outSuffix + "_pairCorrOdd["
        + to_string(pairCorrOdd[0].size()) + "].dat");
    if (fitPC) {
      save::saveDat<double>(pairCorrOdd, 
          outputDir + "/fitPairCorr"    
          + outSuffix + "_maxR-"
          + to_string(params.maxR) + "["
          + to_string(pairCorrOdd[0].size()) + "].dat");
    }

  }


  /////////////////////////
  /////  Cleaning up  /////
  /////////////////////////

  if (params.verbose) 
    std::cout << "Cleaning up" << endl;

  fftw_free(rSpace);
  fftw_free(qSpace);
  fftw_destroy_plan(fftB);
  cout<<"done"<<endl;
  return 1;
}
