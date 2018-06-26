#include "simulations.h"


int main(int argc, char* argv[]) {

  radicalEnums molecule = nitrobenzene; //phenylRadical; //nitrobenzene; //phenoxyRadical;
  //radicalEnums molecule = phenylRadical; //nitrobenzene; //phenoxyRadical;
  //radicalEnums molecule = phenoxyRadical;

  // Diffraction Pattern
  uint NradBins = 30;
  uint Nbins = 2*NradBins - 1;
  //double maxQ = 11.3*((double)Nbins)/1000.0;
  double maxQ = 11.3*835/1024*(9.5/9.0);
  double elEnergy = 3.7e6;
  double Iebeam = 5;
  double screenDist = 4;


  //PLOTclass plt;

  double seed = (double)clock();
  int Nmols = 1;
  int index = 0;
  std::string fileName = radicalNames[molecule];
  std::string outputDir = "output/references";
  //string outputDir = "/reg/neh/home/khegazy/simulations/n2o/diffractionPatterns/output";
  std::string xyzDir = "/reg/neh/home5/khegazy/analysis/nitroBenzene/simulation/XYZfiles/";
  if (argc > 1) {
    for (int iarg=1; iarg<argc; iarg+=2) {
      if (strcmp(argv[iarg],"-Nmols")==0) {Nmols = atoi(argv[iarg+1]);}
      else if (strcmp(argv[iarg],"-Nbins")==0) {NradBins = atoi(argv[iarg+1]);}
      else if (strcmp(argv[iarg],"-Ofile")==0) {string str(argv[iarg+1]); fileName=str;}
      else if (strcmp(argv[iarg],"-Odir")==0) {string str(argv[iarg+1]); outputDir=str;}
      else if (strcmp(argv[iarg],"-maxQ")==0) {maxQ = atof(argv[iarg+1]);}
      else if (strcmp(argv[iarg],"-Index")==0) {index = atoi(argv[iarg+1]);}
      else {
        cerr<<"ERROR!!! Option "<<argv[iarg]<<" does not exist!"<<endl;
        exit(0);
      }
    }
  }

  ////////////////////////////////////////
  /////  Making Diffraction Pattern  /////
  ////////////////////////////////////////

  std::vector<MOLENSEMBLEMCclass*> molMCs;

  switch(molecule) {
    case nitrobenzene: {
      /////  nitrobenzene  /////
      MOLENSEMBLEMCclass* NBZmc = new MOLENSEMBLEMCclass(seed, 
          xyzDir + "17050202_Nitrobenzene_opt_B3LYP_6-31G.xyz");
      NBZmc->Nmols = Nmols;
      molMCs.push_back(NBZmc);
      break;
    }

    case phenoxyRadical: {  
      /////  phenoxy radical  /////
      MOLENSEMBLEMCclass* PNOXYmc = new MOLENSEMBLEMCclass(seed, 
          xyzDir + "18062101_phenyloxy_opt_B3LYP_6-31G.xyz");
      PNOXYmc->Nmols = Nmols;
      molMCs.push_back(PNOXYmc);

      MOLENSEMBLEMCclass* NOmc = new MOLENSEMBLEMCclass(seed,
          xyzDir + "18062102_NO_opt_B3LYP_6-31G.xyz");
      NOmc->Nmols = Nmols;
      molMCs.push_back(NOmc);
      break;
    }

    case phenylRadical: {     
      /////  phenyl radical  /////
      PNLRadMCclass* PNLmc = new PNLRadMCclass(seed);
      PNLmc->Nmols = Nmols;
      PNLmc->NmolAtoms = 11;
      PNLmc->atomTypes.push_back(C);
      PNLmc->atomTypes.push_back(H);
      molMCs.push_back(PNLmc);

      NO2MCclass* NO2mc = new NO2MCclass(seed);
      NO2mc->Nmols = Nmols;
      NO2mc->NmolAtoms = 3;
      NO2mc->atomTypes.push_back(O);
      NO2mc->atomTypes.push_back(N);
      molMCs.push_back(NO2mc);
      break;
    }

    default:                
      cerr << "ERROR: do not recognize molecule enum!!!\n";
      exit(0);
  }


  std::vector< std::vector<double> > curLineOuts;
  std::map< std::string, std::vector<double> > lineOuts;
  std::map< std::string, std::vector< std::vector<double> > > diffPatterns;
  for (auto mc : molMCs) {
    mc->useOrientationMC = false;
    mc->makeMolEnsemble();

    DIFFRACTIONclass diffP(mc, maxQ, Iebeam, screenDist, elEnergy, Nbins,
        "/reg/neh/home/khegazy/simulations/scatteringAmplitudes/3.7MeV/");
    
    curLineOuts.clear();
    curLineOuts = diffP.azimuthalAvg_uniform();
    diffP.diffPatternCalc_uniform();

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

    // fill lineouts
    if (lineOuts.size() == 0) {
      lineOuts["diffractionPattern"]    = curLineOuts[0];
      lineOuts["molDiffractionPattern"] = curLineOuts[1];
      lineOuts["atmDiffractionPattern"] = curLineOuts[2];
      lineOuts["sPattern"]              = curLineOuts[3];
    }
    else {
      for (uint i=0; i<lineOuts["diffractionPattern"].size(); i++) {
        lineOuts["diffractionPattern"][i]    += curLineOuts[0][i];
        lineOuts["molDiffractionPattern"][i] += curLineOuts[1][i];
        lineOuts["atmDiffractionPattern"][i] += curLineOuts[2][i];
      }
    }
  }


  /////////////////////////////////////////
  /////  Plotting and Saving Results  /////
  /////////////////////////////////////////

  std::string prefix = outputDir + "/" + radicalNames[molecule] + "_";
  std::string suffixLO = "Bins-" + to_string(NradBins) + "_Qmax-" + to_string(maxQ)
                        + "_Ieb-" + to_string(Iebeam) + "_scrnD-"
                        + to_string(screenDist) + "_elE-" + to_string(elEnergy);
  std::string suffixDP = "Bins-" + to_string(Nbins) + "_Qmax-" + to_string(maxQ)
                        + "_Ieb-" + to_string(Iebeam) + "_scrnD-"
                        + to_string(screenDist) + "_elE-" + to_string(elEnergy);

  std::string lineOutSpan = "0," + to_string(maxQ);
  std::vector<PLOToptions> opts(3);
  std::vector<std::string> vals(3);
  opts[0] = xSpan;      vals[0] = to_string(-maxQ) + "," + to_string(maxQ);
  opts[1] = ySpan;      vals[1] = to_string(-maxQ) + "," + to_string(maxQ);
  opts[2] = fileType;   vals[2] = "png";


  /*
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
      */


  cout<<"opening files"<<endl;
  FILE* otpDp = fopen((prefix + "diffractionPattern_" + suffixDP + ".dat").c_str(), "wb");
  FILE* otpAp = fopen((prefix + "atmDiffractionPattern_" + suffixDP + ".dat").c_str(), "wb");
  FILE* otpMp = fopen((prefix + "molDiffractionPattern_" + suffixDP + ".dat").c_str(), "wb");
  FILE* otpSp = fopen((prefix + "sPattern_" + suffixDP + ".dat").c_str(), "wb");
  FILE* otpDpLo = fopen((prefix + "diffractionPatternLineOut_" + suffixLO + ".dat").c_str(), "wb");
  FILE* otpMpLo = fopen((prefix + "molDiffractionPatternLineOut_" + suffixLO + ".dat").c_str(), "wb");
  FILE* otpApLo = fopen((prefix + "atmDiffractionPatternLineOut_" + suffixLO + ".dat").c_str(), "wb");
  FILE* otpSpLo = fopen((prefix + "sPatternLineOut_" + suffixLO + ".dat").c_str(), "wb");

  cout<<"writing to files"<<endl;
  for (uint ir=0; ir<diffPatterns["diffractionPattern"].size(); ir++) {
    fwrite(&diffPatterns["diffractionPattern"][ir][0], sizeof(double), 
        diffPatterns["diffractionPattern"][ir].size(), otpDp);
    fwrite(&diffPatterns["molDiffractionPattern"][ir][0], sizeof(double), 
        diffPatterns["molDiffractionPattern"][ir].size(), otpMp);
    fwrite(&diffPatterns["atmDiffractionPattern"][ir][0], sizeof(double), 
        diffPatterns["atmDiffractionPattern"][ir].size(), otpAp);
    fwrite(&diffPatterns["sPattern"][ir][0], sizeof(double), 
        diffPatterns["sPattern"][ir].size(), otpSp);
  }

  cout<<"closing files"<<endl;
  fwrite(&lineOuts["diffractionPattern"][0], sizeof(double), 
      lineOuts["diffractionPattern"].size(), otpDpLo);
  fwrite(&lineOuts["molDiffractionPattern"][0], sizeof(double), 
      lineOuts["molDiffractionPattern"].size(), otpMpLo);
  fwrite(&lineOuts["atmDiffractionPattern"][0], sizeof(double), 
      lineOuts["atmDiffractionPattern"].size(), otpApLo);
  fwrite(&lineOuts["sPattern"][0], sizeof(double), 
      lineOuts["sPattern"].size(), otpSpLo);


  //////////////////////
  /////  Clean up  /////
  //////////////////////
 
  cout<<"cleaning"<<endl;
  for (auto& mc : molMCs) {
    delete mc;
  }
  cout<<"files"<<endl;

  fclose(otpDp);
  fclose(otpAp);
  fclose(otpMp);
  fclose(otpSp);
  fclose(otpDpLo);
  fclose(otpApLo);
  fclose(otpMpLo);
  fclose(otpSpLo);
  cout<<"done"<<endl;


  return 1;
}
