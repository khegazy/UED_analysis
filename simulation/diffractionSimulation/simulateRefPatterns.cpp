#include "/reg/neh/home/khegazy/baseScripts/atomClass.h"
#include "/reg/neh/home/khegazy/baseScripts/moleculeClass.h"
#include "/reg/neh/home/khegazy/baseScripts/molEnsembleMC.h"
#include "/reg/neh/home/khegazy/baseScripts/diffractionClass.h"
#include "/reg/neh/home/khegazy/baseScripts/plotClass.h"
#include "/reg/neh/home/khegazy/baseScripts/imageProcessing.h"
#include "/reg/neh/home/khegazy/baseScripts/saveClass.h"
#include "../../parameters.h"


int main(int argc, char* argv[]) {

  ////////////////////////
  /////  Parameters  /////
  ////////////////////////

  PLOTclass plt;
  parameterClass params("simulateReference");
  uint Nbins = 2*params.NradSimBins - 1;
  int Nmols = 1;
  int index = 0;
  double maxQ = params.NradSimBins*params.QperPix;
  double seed = (double)clock();
  std::string fileName = params.radicalNames[params.molecule];
  if (params.verbose)
    std::cout << "INFO: Parsing command line input\n";
  if (argc > 1) {
    for (int iarg=1; iarg<argc; iarg+=2) {
      if (strcmp(argv[iarg],"-Nmols")==0) 
          {Nmols = atoi(argv[iarg+1]);}
      else if (strcmp(argv[iarg],"-Nbins")==0) 
          {params.NradSimBins = atoi(argv[iarg+1]);
           Nbins = 2*params.NradSimBins - 1;}
      else if (strcmp(argv[iarg],"-Ofile")==0) 
          {string str(argv[iarg+1]); fileName=str;}
      else if (strcmp(argv[iarg],"-Odir")==0) 
          {string str(argv[iarg+1]); params.simOutputDir=str;}
      else if (strcmp(argv[iarg],"-QperPix")==0) 
          {params.QperPix = atof(argv[iarg+1]); 
           maxQ = params.NradSimBins*params.QperPix;}
      else if (strcmp(argv[iarg],"-maxQ")==0) 
          {maxQ = atof(argv[iarg+1]);}
      else if (strcmp(argv[iarg],"-Index")==0) 
          {index = atoi(argv[iarg+1]);}
      else {
        std::cerr << "ERROR!!! Option " << argv[iarg] << " does not exist!\n";
        exit(0);
      }
    }
  }


  //////////////////////////////////////////////////
  /////  Begin Diffraction Pattern Simulation  /////
  //////////////////////////////////////////////////

  /////  Building molecular ensembles  /////
  if (params.verbose)
    std::cout << "INFO: Building molecular ensembles\n";
  std::vector<MOLENSEMBLEMCclass*> molMCs;
  for (auto& fileName : params.xyzFiles) {
    MOLENSEMBLEMCclass* molMC = new MOLENSEMBLEMCclass(seed, 
        params.xyzDir + "/" + fileName);
    molMC->Nmols = Nmols;
    molMC->useOrientationMC = false;
    molMC->makeMolEnsemble();
    molMCs.push_back(molMC);
    cout<<"made nbzmol"<<endl;
  }

  /////  Make diffraction patterns and lineouts  /////
  if (params.verbose)
    std::cout << "INFO: Looping over MCs to make lineouts\n";
  //ofstream outtxt("cppSim.txt");
  std::vector<double> pairCorr;
  std::vector< std::vector<double> > curLineOuts;
  std::map< std::string, std::vector<double> > lineOuts;
  std::map< std::string, std::vector< std::vector<double> > > diffPatterns;
  for (auto mc : molMCs) {
    DIFFRACTIONclass diffP(mc, 
        maxQ, 
        params.Iebeam, 
        params.screenDist, 
        params.elEnergy, 
        Nbins,
        "/reg/neh/home/khegazy/simulations/scatteringAmplitudes/3.7MeV/");
    
    curLineOuts.clear();
    curLineOuts = diffP.azimuthalAvg_uniform();
    //diffP.diffPatternCalc_uniform();

    // fill diffraction patterns
    if (diffPatterns.size() == 0) {
      diffPatterns["diffractionPattern"]    = diffP.diffPattern;
      diffPatterns["molDiffractionPattern"] = diffP.diffMolPattern;
      diffPatterns["atmDiffractionPattern"] = diffP.diffAtmPattern;
      diffPatterns["sPattern"]              = diffP.sPattern;
    }
    else {
      for (uint ir=0; ir<diffP.diffPattern.size(); ir++) {
        for (uint ic=0; ic<diffP.diffPattern[ir].size(); ic++) {
          diffPatterns["diffractionPattern"][ir][ic]    += diffP.diffPattern[ir][ic];
          diffPatterns["molDiffractionPattern"][ir][ic] += diffP.diffMolPattern[ir][ic];
          diffPatterns["atmDiffractionPattern"][ir][ic] += diffP.diffAtmPattern[ir][ic];
        }
      }
    }

    cout<<"lineOuts"<<endl;

    // fill lineouts
    if (lineOuts.size() == 0) {
      lineOuts["diffractionPattern"]    = curLineOuts[0];
      lineOuts["molDiffractionPattern"] = curLineOuts[1];
      lineOuts["atmDiffractionPattern"] = curLineOuts[2];
      lineOuts["sPattern"]              = curLineOuts[3];
      /*
      for (uint i=0; i<lineOuts["diffractionPattern"].size(); i++) {
        outtxt<<i<<"  ";
        outtxt << lineOuts["sPattern"][i] << " ";
        outtxt << lineOuts["atmDiffractionPattern"][i] << "  ";
        outtxt << lineOuts["molDiffractionPattern"][i] << endl;
      }
    outtxt.close();
    */
    }
    else {
      for (uint i=0; i<lineOuts["diffractionPattern"].size(); i++) {
        lineOuts["diffractionPattern"][i]    += curLineOuts[0][i];
        lineOuts["molDiffractionPattern"][i] += curLineOuts[1][i];
        lineOuts["atmDiffractionPattern"][i] += curLineOuts[2][i];
      }
    }

    if (params.simPairCorr) {
      std::vector<double> pCorr = mc->simulatePairCorr(500, params.maxR, false, NULL);

      if (!pairCorr.size()) {
        pairCorr.resize(pCorr.size(), 0);
      }
      for (uint i=0; i<pCorr.size(); i++) {
        pairCorr[i] += pCorr[i];
      }
    }
    plt.print1d(pairCorr, "./plots/testInit");

  }

  /////  Compute sMs pattern  /////
  lineOuts["sMsPattern"].resize(lineOuts["sPattern"].size(), 0);
  for (uint i=0; i<lineOuts["sPattern"].size(); i++) {
    lineOuts["sMsPattern"][i] = lineOuts["sPattern"][i]
        *lineOuts["molDiffractionPattern"][i]
        /lineOuts["atmDiffractionPattern"][i];
  }

  cout<<"SIZES: "<<curLineOuts[0].size()<<endl;

  /////////////////////////////////////////
  /////  Plotting and Saving Results  /////
  /////////////////////////////////////////

  cout<<"testing qmax: "<<maxQ<<"  "<<params.NradSimBins<<"  "<<params.QperPix<<endl;
  std::string prefix = params.simOutputDir + "/" + params.molName + "_";
  std::string suffixLO = 
      "Qmax-" + to_string(maxQ)
      + "_Ieb-" + to_string(params.Iebeam) 
      + "_scrnD-" + to_string(params.screenDist) 
      + "_elE-" + to_string(params.elEnergy)
      + "_Bins[" + to_string(params.NradSimBins) + "]"; 
  std::string suffixDP = 
      "Qmax-" + to_string(maxQ)
      + "_Ieb-" + to_string(params.Iebeam) 
      + "_scrnD-" + to_string(params.screenDist) 
      + "_elE-" + to_string(params.elEnergy)
      + "_Bins[" + to_string(Nbins) + "]";

  /////  Plotting  /////
  if (params.simPltVerbose) {
    std::string lineOutSpan = "0," + to_string(maxQ);
    std::vector<PLOToptions> opts(3);
    std::vector<std::string> vals(3);
    opts[0] = xSpan;      vals[0] = to_string(-maxQ) + "," + to_string(maxQ);
    opts[1] = ySpan;      vals[1] = to_string(-maxQ) + "," + to_string(maxQ);
    opts[2] = fileType;   vals[2] = "png";

    delete (TH2F*)plt.printRC(diffPatterns["diffractionPattern"], 
        prefix + "diffractionPattern_" + suffixDP, opts, vals);
    delete (TH2F*)plt.printRC(diffPatterns["molDiffractionPattern"], 
        prefix + "molDiffractionPattern_" + suffixDP, opts, vals);
    delete (TH2F*)plt.printRC(diffPatterns["atmDiffractionPattern"], 
        prefix + "atmDiffractionPattern_" + suffixDP, opts, vals);
    delete (TH2F*)plt.printRC(diffPatterns["sPattern"], 
        prefix + "sPattern_" + suffixDP, opts, vals);
    delete (TH1F*)plt.print1d(lineOuts["diffractionPattern"], 
        prefix + "diffractionPatternLineOut_" + suffixLO, xSpan, lineOutSpan);
    delete (TH1F*)plt.print1d(lineOuts["molDiffractionPattern"], 
        prefix + "molDiffractionPatternLineOut_" + suffixLO, xSpan, lineOutSpan);
    delete (TH1F*)plt.print1d(lineOuts["atmDiffractionPattern"], 
        prefix + "atmDiffractionPatternLineOut_" + suffixLO, xSpan, lineOutSpan);
    delete (TH1F*)plt.print1d(lineOuts["sPattern"], 
        prefix + "sPatternLineOut_" + suffixLO, xSpan, lineOutSpan);
  }


  /////  Saving  /////

  save::saveDat<double>(lineOuts["diffractionPattern"], 
      prefix + "diffractionPatternLineOut_" + suffixLO + ".dat");
  save::saveDat<double>(lineOuts["molDiffractionPattern"],
      prefix + "molDiffractionPatternLineOut_" + suffixLO + ".dat");
  save::saveDat<double>(lineOuts["atmDiffractionPattern"],
      prefix + "atmDiffractionPatternLineOut_" + suffixLO + ".dat");
  save::saveDat<double>(lineOuts["sPattern"],
      prefix + "sPatternLineOut_" + suffixLO + ".dat");
  save::saveDat<double>(lineOuts["sMsPattern"],
      prefix + "sMsPatternLineOut_" + suffixLO + ".dat");
  save::saveDat<double>(pairCorr,
      prefix + "pairCorr_Bins[" + to_string(pairCorr.size()) + "].dat");


  //////////////////////
  /////  Clean up  /////
  //////////////////////
 
  cout<<"cleaning"<<endl;
  for (auto& mc : molMCs) {
    delete mc;
  }

  return 1;
}
