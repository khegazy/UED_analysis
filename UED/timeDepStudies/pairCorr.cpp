#include "../analysis.h"
#include "../../parameters.h"

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

  bool givenFileName = false;
  std::string fileName = "NULL";
  std::string inputDir = "NULL";
  std::string outPrefix = "data-";
  std::string outSuffix = "";

  //double filterVar = std::pow(params.NradAzmBins/3, 2);
  double filterVar = std::pow(params.NradAzmBins/2.5, 2);
  cout<<"filterVar: "<<filterVar<<endl;
  /////  Importing variables from command line  /////
  for (int iarg=2; iarg<argc; iarg+=2) {
    if (strcmp(argv[iarg], "-Odir") == 0) {
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
    else if (strcmp(argv[iarg], "-Osuf") == 0) {
      string str(argv[iarg+1]);
      outSuffix = str;
    }
    else {
      cerr<<"ERROR!!! Option "<<argv[iarg]<<" does not exist!"<<endl;
      exit(0);
    }
  }

  if ((inputDir.find("sim") != string::npos) 
      || (inputDir.find("Sim") != string::npos)) {
    outPrefix = "sim-";
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

  std::vector<int> shape;
 
  // Importing Data
  if (!givenFileName) {
    shape = save::getShape( 
        params.mergeScansOutputDir,
        "data-" + runName + "-sMsAzmAvgDiff");

    assert(shape[1] == params.NradAzmBins);

    fileName = params.mergeScansOutputDir 
        + "data-" + runName + "-sMsAzmAvgDiff[" 
        + to_string(shape[0])
        + "," + to_string(shape[1]) + "].dat";
  }
  else {
    assert(strcmp(inputDir.c_str(), "NULL") != 0);

    shape = save::getShape(inputDir, fileName);
    if (shape.size() == 1) {
      shape.insert(shape.begin(), 1);
    }

    params.NradAzmBins = shape[1];
    cout<<"shape: "<<shape[0]<<"  "<<shape[1]<<endl;
    
    fileName = inputDir + "/" + fileName;
  }

  std::vector< std::vector<double> > tmdDiff;
  tmdDiff.resize(shape[0]);
  for (auto& v : tmdDiff) {
    v.resize(shape[1], 0);
  }

  if (params.verbose)
    std::cout << "Retriving data from " << fileName << endl; 

  save::importDat<double>(tmdDiff, fileName);
    
  if (params.pltVerbose) {
    std::vector<PLOToptions> pOpts(2);
    std::vector<std::string> pVals(2);
    pOpts[0] = minimum;  pVals[0] = "-2e-2";
    pOpts[1] = maximum;  pVals[1] = "2e-2";
    plt.printRC(tmdDiff, "./plots/importedTDdiff", pOpts, pVals);
  }

  // Importing Simulation
  std::vector<double> sMsSim(shape[1]), sMsSimDer(shape[1]);
  std::string sMsSimName = 
      params.simReferenceDir + "nitrobenzene_sMsPatternLineOut_Qmax-"
      + to_string(params.maxQazm) + "_Ieb-"
      + to_string(params.Iebeam) + "_scrnD-"
      + to_string(params.screenDist) + "_elE-"
      + to_string(params.elEnergy) + "_Bins["
      + to_string(shape[1]) + "].dat";
  if (!tools::fileExists(sMsSimName)) {
    std::cerr << "ERROR: sMs simulation file " + sMsSimName + " does not exist!!!\n";
    exit(0);
  }
  
  save::importDat<double>(sMsSim, sMsSimName);

  for (int iq=1; iq<shape[1]-1; iq++) {
    sMsSimDer[iq] = (sMsSimDer[iq+1] - sMsSimDer[iq-1])/2;
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
  double w = PI/params.NradAzmBins;
  std::vector<double> fftInp(inpFFTsize);
  std::vector<double> smoothed(params.NradAzmBins, 0);
  std::vector<double> filter(params.NradAzmBins, 0);
  std::vector< std::vector<double> > results;
  std::vector< std::vector<double> > pairCorr(shape[0]);
  std::vector< std::vector<double> > pairCorrOdd(shape[0]);
  for (uint i=0; i<shape[0]; i++) {
    pairCorr[i].resize(outSize);
    pairCorrOdd[i].resize(outSize);
  }

  double filterScale = 1;
  for (int itm=0; itm<shape[0]; itm++) {

    /////  Filling FFT Input  /////
    for (int iq=0; iq<params.NradAzmBins; iq++) {
      if (iq < params.NbinsSkip) {
        tmdDiff[itm][iq] = sMsSim[iq]
            *(tmdDiff[itm][params.NbinsSkip]/sMsSim[params.NbinsSkip]);
      }

      //  Apply Filter
      filterScale = exp(-1*std::pow(iq, 2)
                        /(2*filterVar));
      smoothed[iq] = tmdDiff[itm][iq]*filterScale;
      filter[iq] = filterScale;
    }
    if (params.pltVerbose) {
      plt.print1d(tmdDiff[itm], "./plots/inpDiff_" + runName + "_" + to_string(itm));
      plt.print1d(smoothed, "./plots/smoothed_" + runName + "_" + to_string(itm));
      plt.print1d(filter, "./plots/filter");
    }

    /////  Apply FFT  /////
    if (mirrorDiffPattern) {
      fftInp[FFTcent] = smoothed[0];
      for (int iq=1; iq<params.NradAzmBins; iq++) {
        fftInp[FFTcent+iq] = smoothed[iq];
        fftInp[FFTcent-iq] = smoothed[iq];
      }
      results = tools::fft1dRtoC(fftInp, fftB, qSpace, rSpace, 
                    false, true, NULL, "fftInp_" + to_string(itm));
    }
    else {
      for (int iq=0; iq<params.NradAzmBins; iq++) {
        fftInp[iq] = smoothed[iq];
      }
      results = tools::fft1dRtoC(fftInp, fftB, qSpace, rSpace, 
                    false, false, &plt, "fftInp_" + to_string(itm));
    }
    if (params.pltVerbose)
      plt.print1d(fftInp, "./plots/fftFuncInp_" + to_string(itm));

    for (int ir=0; ir<outSize; ir++) {
      pairCorr[itm][ir] = std::pow(results[0][ir], 1);
      pairCorrOdd[itm][ir] = std::pow(results[1][ir], 1);
    }
  }

  if (params.pltVerbose) {
    vals[5] = "0," + to_string(params.rMaxAzmRat*outFFTsize*(2*PI/(params.QperPix*inpFFTsize)));
    plt.printRC(pairCorr, "./plots/pairCorr", opts, vals);
    plt.printRC(pairCorrOdd, "./plots/pairCorrOdd", opts, vals);
  }


  ////////////////////
  /////  Saving  /////
  ////////////////////

  if (shape[0] != 1) {
    save::saveDat<double>(pairCorr, 
        "./results/" + outPrefix 
        + runName + outSuffix + "-pairCorr["
        + to_string(pairCorr.size()) + ","
        + to_string(pairCorr[0].size()) + "].dat");
    save::saveDat<double>(pairCorrOdd, 
        "./results/" + outPrefix   
        + runName + outSuffix + "-pairCorrOdd["
        + to_string(pairCorr.size()) + ","
        + to_string(pairCorr[0].size()) + "].dat");
  }
  else {
    save::saveDat<double>(pairCorr[0], 
        "./results/" + outPrefix 
        + runName + outSuffix + "-pairCorr["
        + to_string(pairCorr[0].size()) + "].dat");
    save::saveDat<double>(pairCorrOdd[0], 
        "./results/" + outPrefix 
        + runName + outSuffix + "-pairCorrOdd["
        + to_string(pairCorr[0].size()) + "].dat");
  }


  /////////////////////////
  /////  Cleaning up  /////
  /////////////////////////

  if (params.verbose) 
    std::cout << "Cleaning up" << endl;

  fftw_free(rSpace);
  fftw_free(qSpace);
  fftw_destroy_plan(fftB);
  return 1;
}
