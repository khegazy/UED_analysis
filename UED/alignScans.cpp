#include "plotLegendres.h"
#include "../simulation/diffractionPattern/simulations.h"
#include <TLegend.h>


using namespace std;


int main(int argc, char* argv[]) {

  /////  Flags  /////
  bool compareFinalStates = false;



  /////  Compare <t0 and simulation  /////
  //std::pair<double, double> t0SimCompTimes = {

  if (argc<2) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe 'fileList.txt' !!!"<<endl;
    cerr<<"         Can also run ./analysis.exe 'fileList.txt' 'treeName'"<<endl;
    exit(0);
  }

  string fileList(argv[1]);
  string treeName("physics");
  if (argc==3) string treeName(argv[2]);
  analysisClass analysis(fileList, treeName);  // Use this when specifying treeName


  ///// Load environment and get the number of events /////

  uint64_t Nentries;
  Nentries = analysis.setupEnvironment();  // Alter this function as needed for specific setup

  cout.setf(ios::scientific);
  

  /////  Setting up variables  /////

  bool verbose = false;

  const double maxQ = 9.726264; //11.3;
  const int NradBins = 30;
  const int Nlegendres = 5;
  std::vector<double> atmDiff(NradBins);
  std::vector<double> molDiff(NradBins);
 
  string fileName = "alignment.txt";
  string outputDir = "output/data/";
  for (int iarg=2; iarg<argc; iarg+=2) {
    if (strcmp(argv[iarg],"-Ofile")==0) {
      string str(argv[iarg+1]); 
      fileName = str + ".txt";
    }
    else if (strcmp(argv[iarg],"-Odir")==0) {
      string str(argv[iarg+1]); 
      outputDir = str;
    }
    else {
      cerr<<"ERROR!!! Option "<<argv[iarg]<<" does not exist!"<<endl;
      exit(0);
    }
  }

  double timeZero = 0.5;
  bool subtractT0 = true;

  int minImgCount = 5;
  double stdScale = 1.5;

  int curRun = -1;
  string runInd = "";
  string curDate = "";
  string curScan = "";
  std::vector<PLOToptions> opts(5);
  std::vector<std::string> vals(5);
  std::vector<PLOToptions> oppts(3);
  std::vector<std::string> vaals(3);
  oppts[0] = yLabel;   vaals[0] = "Time [ps]";
  oppts[1] = xLabel;   vaals[1] = "Scattering Q [arb units]";
  oppts[2] = draw;     vaals[2] = "COLZ";
  opts[0] = yLabel;   vals[0] = "Time [ps]";
  opts[1] = xLabel;   vals[1] = "Scattering Q [arb units]";
  opts[4] = draw;     vals[4] = "COLZ";
  opts[2] = minimum;
  opts[3] = maximum;

  std::map<string, double > cbar;
  cbar["20161102_LongScan1_Leg0"] = 0.2;
  cbar["20161102_LongScan1_Leg2"] = 0.2;


  std::map< string, std::vector<double> > refAutCor, refLineOut;

  double Iebeam     = 5;
  double elEnergy   = 3.7e6;
  double screenDist = 4;
  std::string simReferenceDir = "/reg/ued/ana/scratch/nitroBenzene/simulations/references/";
  //std::string simReferenceDir = "../simulation/diffractionPattern/output/references/";

  double stdev = 0.025;

  std::vector< std::vector< std::vector<double> > > legendres(Nlegendres), smearedImg(Nlegendres);
  std::map< int32_t, int > stagePosInds;
  std::map< int, std::vector<double> > scanCounts;
  std::map< int, std::vector< std::vector<double> > > scanLgndrs;


  
  
  std::map< int32_t, std::vector<double> > diffPs, legCoeff_map, pairCorr_map;
  std::map< int32_t, std::vector< std::vector<double> > > avgImgs_map;
  std::map< int32_t, double > counts;
  std::map< int32_t, int > Nimages;
  std::map< int32_t, std::vector<int32_t> > allPos;

  int NautCpadding = 1000;
  int FFTsize = NradBins*2 + NautCpadding + 1;
  //int FFTsize = (NdiffInds + NautCpadding)*2 + 1;
  int FTsize = (int)(FFTsize/2) + 1;
  double holeRat = 0.15;
  double rMaxRat = 0.75;
  double padDecayRat = 0.5;
  int indHR =  holeRat*NradBins;
  int outSize = FTsize*rMaxRat;
  double rMax = (rMaxRat*NradBins)*(2*PI/(2*11.3));
  cout<<"rmax  "<<rMax<<endl;

  double* qSpace = (double*) fftw_malloc(sizeof(double)*FFTsize);
  fftw_complex* rSpace = 
        (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(int)((FFTsize/2) + 1));
  fftw_plan fftB;
  fftB = fftw_plan_dft_r2c_1d(FFTsize, qSpace, rSpace, FFTW_MEASURE);


  ////////////////////////////////////////////
  /////  Importing simulated references  /////
  ////////////////////////////////////////////

  std::string fileNameSuffix = "_Bins-" + to_string(NradBins) + "_Qmax-" + to_string(maxQ)
                                + "_Ieb-" + to_string(Iebeam) + "_scrnD-"
                                + to_string(screenDist) + "_elE-" + to_string(elEnergy) + ".dat";

  ///  Nitrobenzene reference  ///
  FILE* atmFile = fopen((simReferenceDir + "/nitrobenzene_atmDiffractionPatternLineOut"
                           + fileNameSuffix).c_str(), "rb");
  fread(&atmDiff[0], sizeof(double), atmDiff.size(), atmFile);
  fclose(atmFile);
  plt.print1d(atmDiff, "testATMNBZ");
  FILE* molFile = fopen((simReferenceDir + "/nitrobenzene_molDiffractionPatternLineOut"
                           + fileNameSuffix).c_str(), "rb");
  fread(&molDiff[0], sizeof(double), molDiff.size(), molFile);
  fclose(molFile);

  // sMs normalization
  std::vector<double> sMsNorm(NradBins);
  for (int ir=0; ir<NradBins; ir++) {
    sMsNorm[ir] = (maxQ*(ir+0.5)/NradBins)/(1e20*atmDiff[ir]);
  }

  ///  Various possible final states  ///
  if (compareFinalStates) {
    for (auto name : radicalNames) {

      cout<<"importing refs"<<endl;
      // Import data
      refLineOut[name].resize(NradBins, 0);

      FILE* refFile = fopen((simReferenceDir + "/" + name + "_molDiffractionPatternLineOut" 
                          + fileNameSuffix).c_str(), "rb");
      fread(&refLineOut[name][0], sizeof(double), NradBins, refFile);
      fclose(refFile);

      //delete plt.print1d(refLineOut[name], "./plots/sim_" + name + "_diffraction");
      // Normalize by atomic scattering
      for (int ir=0; ir<NradBins; ir++) {
        refLineOut[name][ir] *= sMsNorm[ir];
      }
      save::saveDat<double>(refLineOut[name], "./plots/data/sim_" + name + "Diffraction["
            + to_string(NradBins) + "].dat");
      //delete plt.print1d(refLineOut[name], "./sim_" + name + "_diffractionsMs");

      cout<<"start pair corr"<<endl;
      // Pair correlation 
      std::vector<double> inpDiff(NradBins*2+1, 0);
      std::vector< std::vector<double> > autCor1d;
      int centI = (int)(inpDiff.size()/2);
      cout<<"filling inpdiff"<<endl;
      for (int ir=0; ir<NradBins; ir++) {
        if (ir < indHR) {
          inpDiff[centI+1+ir] = refLineOut[name][indHR]
                *pow(sin((PI/2)*((ir+1)/((double)(indHR+1)))), 2);
          //inpDiff[centI+1+ir] = smearedImg[0][it][ir]
          //      *pow(sin((PI/2)*((ir+1)/((double)(indHR+1)))), 2);
          inpDiff[centI-1-ir] = inpDiff[centI+1+ir];
        }
        else {
          inpDiff[centI+1+ir] = refLineOut[name][ir];
          inpDiff[centI-1-ir] = refLineOut[name][ir];
        }
      }
      cout<<"filled inpdiff"<<endl;
      //cout<<"filled"<<endl;
      //plt.print1d(inpDiff, "inpDiff");

      cout<<"autcor"<<endl;
      autCor1d = tools::fft1dRtoC(inpDiff, rMaxRat, 
              NautCpadding, padDecayRat, fftB, qSpace, rSpace, false);
      refAutCor[name] = autCor1d[0];
      delete plt.print1d(refAutCor[name], name + "_pairCorrSim");
      save::saveDat<double>(refAutCor[name], "./plots/data/sim_" + name + "PairCorr["
              + to_string(refAutCor[name].size()) + "].dat");
    }

    ///  Create difference patterns  ///
    std::vector<double> autDiff(refAutCor[radicalNames[0]].size(), 0);
    std::vector<double> diffDiff(NradBins, 0);
    for (uint i=1; i<Nradicals; i++) {
      for (uint ir=0; ir<autDiff.size(); ir++) {
        autDiff[ir] = refAutCor[radicalNames[i]][ir]
                      - refAutCor[radicalNames[0]][ir];
      }
      save::saveDat(autDiff, "./plots/data/sim_" + radicalNames[i] 
          + "PairCorrDiff[" + to_string(autDiff.size()) + "].dat");

      for (uint ir=0; ir<NradBins; ir++) {
        diffDiff[ir] = refLineOut[radicalNames[i]][ir]
                       - refLineOut[radicalNames[0]][ir];
      }
      save::saveDat(diffDiff, "./plots/data/sim_" + radicalNames[i] 
          + "DiffractionDiff[" + to_string(diffDiff.size()) + "].dat");
    }
  }



  ///////////////////////////
  /////  Aligning Runs  /////
  ///////////////////////////

  if (verbose)  {
    std::cout << "Begin to align runs" << endl;
  }

  std::vector<string> runInds;
  std::map< string, std::map< double, double* > >  diffP_arrays;
  std::map< string, int > runShifts;

  ///// Loop through events in the file and saving to maps  /////
  curRun = -1;
  int pInd;
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);

    // Ignore reference images taken before scan
    if (imgIsRef) {
      continue;
    }
    if (stagePos < (t0StagePos - 0.021)*10000) {
      continue;
    }

    ///  Initialize variables for new run  ///
    if (runNum != curRun) {
      curRun = runNum;

      // Initialize legendres
      uint timeSize = stagePosInds.size();
      const std::vector<double> emptyVec((*legCoeffs).size(), 0.0);
      scanLgndrs[runNum].resize(0);
      for (uint i=0; i<timeSize; i++) {
        scanLgndrs[runNum].push_back(emptyVec);
      }

      // Initialize counts
      scanCounts[runNum].resize(timeSize, 0.0);
    }

    ///  Get index of stage position  ///
    auto pos = stagePosInds.find(stagePos);
    pInd = pos->second;

    /// Add new stage position to index lookup table and all runs  ///
    if (pos == stagePosInds.end()) {
      int ind = 0;
      stagePosInds[stagePos] = -1;
      for (auto& itr : stagePosInds) {

        if (itr.second == -1) {
          pInd = ind;
          auto sCitr = scanCounts.begin();
          auto sLitr = scanLgndrs.begin();
          const std::vector<double> emptyVec((*legCoeffs).size(), 0.0);
          while (sCitr != scanCounts.end()) {
            sCitr->second.insert(sCitr->second.begin()+ind, 0);
            sLitr->second.insert(sLitr->second.begin()+ind, emptyVec);
            sCitr++; sLitr++;
          }
        }
        itr.second = ind;
        ind++;
      }
    }

    if (verbose) {
      std::cout << "ind/sizes: " << pInd << "  " 
        << scanCounts[runNum].size() << "  " 
        << scanLgndrs[runNum].size() << "  "
        << scanLgndrs[runNum][pInd].size() << std::endl;
    }

    ///  Add coefficients and counts to maps  ///
    scanCounts[runNum][pInd] = imgNorm;
    for (uint i=0; i<(*legCoeffs).size(); i++) {
      scanLgndrs[runNum][pInd][i] = (*legCoeffs)[i];
    }
  }




  ///////////////////////////////////////////////////////////
  /////  Prune data for outliers and sparse time steps  /////
  ///////////////////////////////////////////////////////////


  /////  Remove time bins with few images and get mean  /////
  std::vector<int32_t> removePos;
  std::vector< std::vector<double> > runMeans(stagePosInds.size());
  std::vector< std::vector<double> > runStdev(stagePosInds.size());
  for (uint i=0; i<runMeans.size(); i++) {
    runMeans[i].resize((*legCoeffs).size(), 0);
    runStdev[i].resize((*legCoeffs).size(), 0);
  }

  // looping over time bins
  for (auto pItr : stagePosInds) {
    int count = 0;
    double norm = 0;
    for (auto sLitr : scanLgndrs) {
      // mean
      for (uint i=0; i<(*legCoeffs).size(); i++) {
        runMeans[pItr.second][i] += sLitr.second[pItr.second][i];
      }

      // counting non zero bins
      if (scanCounts[sLitr.first][pItr.second] > 0) {
        count++;
      }
      norm += 1; //scanCounts[sLitr.first][pItr.second];
    }

    // calculate mean
    for (uint i=0; i<(*legCoeffs).size(); i++) {
      runMeans[pItr.second][i] /= norm;
    }

    // save time bins to cut
    cout<<"COUNTING: "<<pItr.first<<"   "<<count<<endl;
    if (count < minImgCount) {
      removePos.push_back(pItr.first);
    }
  }

  //plt.printXY(runMeans, "testRunMeans", maximum, "5e2");

  ///  Remove time bins  ///
  if (verbose) {
    std::cout << "Removing " << removePos.size() << " bins\n";
    std::cout << "Before removal: " << stagePosInds.size() << "  " 
      << scanLgndrs.begin()->second.size() << "  " 
      << scanCounts.begin()->second.size() << "  "
      << runMeans.size() << "  " << runStdev.size() << std::endl;
  }
  /*
  for (int i=((int)removePos.size())-1; i>=0; i--) {
    int pInd = stagePosInds[removePos[i]];
    for (auto& sLitr : scanLgndrs) {
      sLitr.second.erase(sLitr.second.begin() + pInd);
    }
    for (auto& sCitr : scanCounts) {
      sCitr.second.erase(sCitr.second.begin() + pInd);
    }
    stagePosInds.erase(stagePosInds.find(removePos[i]));

    runMeans.erase(runMeans.begin() + pInd);
    runStdev.erase(runStdev.begin() + pInd);
  }
  */
  int pCount = 0;
  for (auto& pItr : stagePosInds) {
    pItr.second = pCount;
    pCount++;
  }
  if (verbose) {
    std::cout << "After removal: " << stagePosInds.size() << "  " 
      << scanLgndrs.begin()->second.size() << "  " 
      << scanCounts.begin()->second.size() << "  "
      << runMeans.size() << "  " << runStdev.size() << std::endl;
  }
  
  // Calculating time delays from stage position
  if (verbose) {
    std::cout << "calculating time delays: " << runMeans.size() << std::endl;
  }

  int NtimeSteps = stagePosInds.size();
  int tInd = 1;
  double *timeDelays = new double[NtimeSteps + 1];
  for (auto pItr : stagePosInds) {
    timeDelays[tInd] = 2*pItr.first/(3e3);
    cout<<"time: "<<timeDelays[tInd]<<endl;
    tInd++;
  }
  timeDelays[0] = 2*timeDelays[1] - timeDelays[2];
  for (int i=NtimeSteps; i>=0; i--) {
    timeDelays[i] -= timeZero + timeDelays[0];
  }

  save::saveDat<double>(timeDelays, NtimeSteps + 1, 
      "./plots/data/timeDelays["
      + to_string(stagePosInds.size() + 1) + "].dat");
 


  ///  Calculate bin standard deviation  ///
  for (auto pItr : stagePosInds) {
    double norm = 0;
    for (auto sLitr : scanLgndrs) {
      // mean
      for (uint i=0; i<(*legCoeffs).size(); i++) {
        runStdev[pItr.second][i] 
            += std::pow((sLitr.second[pItr.second][i]
                    - runMeans[pItr.second][i]), 2);
        /*
            += scanCounts[sLitr.first][pItr.second]
                  *std::pow((sLitr.second[pItr.second][i]
                    /scanCounts[sLitr.first][pItr.second]
                    - runMeans[pItr.second][i]), 2);
                    */
      }
      norm += 1; //scanCounts[sLitr.first][pItr.second];
    }

    cout<<"\n\nTIME: "<<pItr.first<<endl<<endl;
    // calculate standard deviation
    for (uint i=0; i<(*legCoeffs).size(); i++) {
      runStdev[pItr.second][i] = std::sqrt(runStdev[pItr.second][i]/norm);
      cout<<runMeans[pItr.second][i]<<"/"<<runStdev[pItr.second][i]<<"\t"<<endl;
    }

  }

  
  std::vector<double> unPumped(NradBins, 0);
  int NrefBins = 0;
  if (subtractT0) {
    while (timeDelays[NrefBins+1] < 0) {
      NrefBins++;
    }
  }
  else {
    NrefBins = NtimeSteps;
  }
  cout<<"NrefBins: "<<NrefBins<<endl;
  int rInd;
  double norm;
  for (uint ilg=0; ilg<Nlegendres; ilg++) {
    legendres[ilg].resize(stagePosInds.size());
    cout<<"starting to fill"<<endl;
    auto pItr = stagePosInds.begin();
    for (uint it=0; it<NtimeSteps; it++) {
      legendres[ilg][it].resize(NradBins, 0.0);
      for (uint ir=0; ir<NradBins; ir++) {
        rInd = ilg*NradBins + ir;
        norm = 0;
        for (auto sLiter : scanLgndrs) {
          cout<<"filling: "<<ilg<<"  "<<it<<"  "<<ir<<"  "<<sLiter.first<<"  "<<sLiter.second[it][rInd]<<"  "<<sLiter.second[it][rInd] - runMeans[it][rInd]<<"  "<< stdScale*runStdev[it][rInd]<<endl;
          if (fabs(sLiter.second[it][rInd] - runMeans[it][rInd]) 
              < stdScale*runStdev[it][rInd]) {
            legendres[ilg][it][ir] += sLiter.second[it][rInd];
            norm += scanCounts[sLiter.first][it];
          }
        }
        legendres[ilg][it][ir] /= norm;
      }
      pItr++;
    }
    cout<<"filled for lg: "<<ilg<<endl;

    // Mean/t0 subtraction and normalizing by atomic scattering
    // Raw data
    std::fill(unPumped.begin(), unPumped.end(), 0.0);
    for (int ir=0; ir<NradBins; ir++) {
      for (int tm=0; tm<NrefBins; tm++) {
        unPumped[ir] += legendres[ilg][tm][ir];
      }
      unPumped[ir] /= NrefBins;
    }
    cout<<"unp filled"<<endl;
    for (int ir=0; ir<NradBins; ir++) {
      for (uint tm=0; tm<legendres[ilg].size(); tm++) {
        legendres[ilg][tm][ir] -= unPumped[ir];
        legendres[ilg][tm][ir] *= sMsNorm[ir];
      }
    }
    cout<<"scaled"<<endl;
    if (ilg==0) {
      plt.print1d(unPumped, "testUnp");
    }

    cout<<"saving"<<endl;
    ///  Saving data  ///
    save::saveDat<double>(legendres[ilg], "./plots/data/data_sMsL" 
        + to_string(ilg) + "Diff["
        + to_string(legendres[ilg].size()) + ","
        + to_string(legendres[ilg][0].size()) + "].dat");
    save::saveDat<double>(legendres[ilg][legendres[ilg].size()-1], "./plots/data/data_sMsFinalL" 
        + to_string(ilg) + "Diff["
        + to_string(legendres[ilg][0].size()) + "].dat");


    // Save raw unpumped data to fit background
    if (ilg == 0) {
      save::saveDat<double>(unPumped, "./qScale/results/unPumpedDiffractionL0["
          + to_string(unPumped.size()) + "].dat");
      plt.print1d(unPumped, "unPumpedRaw");
    }


    /////  Smearing time dependence  /////
    if (verbose) {
      std::cout << "Begin smearing legendres in time\n";
    }

    double sum = 0;
    smearedImg[ilg].resize(NtimeSteps);
    for (uint ir=0; ir<NtimeSteps; ir++) {
      smearedImg[ilg][ir].resize(NradBins, 0);
    }
    for (int ir=0; ir<NradBins; ir++) {
      for (int tm=0; tm<NtimeSteps; tm++) {
        sum = 0;
        for (int tmm=0; tmm<NtimeSteps; tmm++) {
          smearedImg[ilg][tm][ir] += legendres[ilg][tmm][ir]
                                    *std::exp(-1*std::pow(timeDelays[tm] - timeDelays[tmm], 2)/(2*std::pow(stdev, 2)));
          sum += exp(-1*std::pow(timeDelays[tm] - timeDelays[tmm], 2)/(2*std::pow(stdev, 2)));
        }
        smearedImg[ilg][tm][ir] /= sum;
      }
    }

    ///  Subtract before t0  ///
    std::fill(unPumped.begin(), unPumped.end(), 0.0);
    for (int ir=0; ir<NradBins; ir++) {
      for (int tm=0; tm<NrefBins; tm++) {
        unPumped[ir] += smearedImg[ilg][tm][ir];
      }
      unPumped[ir] /= (double)NrefBins;
    }
    for (int ir=0; ir<NradBins; ir++) {
      for (uint tm=0; tm<NtimeSteps; tm++) {
        smearedImg[ilg][tm][ir] -= unPumped[ir];
      }
    }

    if (ilg == 0) {
      save::saveDat<double>(unPumped, "./plots/data/data_unPumpedDiffractionL0Smeared["
          + to_string(NradBins) + "].dat");
    }

    save::saveDat<double>(smearedImg[ilg], "./plots/data/data_sMsL"
        + to_string(ilg) + "DiffSmear"
        + to_string(stdev) + "["
        + to_string(smearedImg[ilg].size()) + ","
        + to_string(smearedImg[ilg][0].size()) + "].dat");
    
    save::saveDat<double>(smearedImg[ilg][NtimeSteps-1], 
        "./plots/data/data_sMsFinalL" 
        + to_string(ilg) + "DiffSmear"
        + to_string(stdev) + "["
        + to_string(smearedImg[ilg][0].size()) + "].dat");

    cout<<"last entry"<<endl;
    cout<<"plotting"<<endl;

    /*
      vals[2] = "-8";
      vals[3] = "8";
      plt.printXY(img[i], timeDelays, "plotRun" + to_string(curRun) + "_" + curDate + "_" + curScan 
            + "_Leg" + to_string(i), oppts, vaals);
    */
      
    ///  Plotting time dependend legendre signal  ///
    opts.resize(3);
    vals.resize(3);
    opts[0] = yLabel;   vals[0] = "Time [ps]";
    opts[1] = xLabel;   vals[1] = "Scattering Q [arb units]";
    opts[2] = draw;     vals[2] = "COLZ";


    std::string cbarLegName = curDate + "_" + curScan + "_Leg" + to_string(ilg); 
    if (cbar.find(cbarLegName) != cbar.end()) {
      opts.push_back(minimum);
      vals.push_back(to_string(-1*cbar[cbarLegName]));
      opts.push_back(maximum);
      vals.push_back(to_string(cbar[cbarLegName]));
    }

    plt.printRC(legendres[ilg], timeDelays, 0, maxQ, curDate + "_" + curScan 
          + "_Leg" + to_string(ilg), opts, vals);
          //+ "_Leg" + to_string(ilg) + "_run" + to_string(curRun), oppts, vaals);
    plt.printRC(smearedImg[ilg], timeDelays, 0, maxQ, "smeared_" + curDate + "_" + curScan
          + "_Leg" + to_string(ilg), opts, vals);
          //+ "_Leg" + to_string(ilg) + "_run" + to_string(curRun), oppts, vaals);

    
    cout<<"plotted"<<endl;
  }

  cout<<"start autocorrelation"<<endl;

  ///////////////////////////////////////
  /////  Pair correlation function  /////
  ///////////////////////////////////////
  if (verbose) {
    std::cout << "Begin calculating pair correlation functions\n";
  }

  std::vector< std::vector<double> > pairCorr(smearedImg[0].size()), pairCorr1d;
  std::vector<double> powerSpct(outSize);
  for (int it=0; it<NtimeSteps; it++) {
    pairCorr[it].resize(outSize);

    std::vector<double> inpDiff(smearedImg[0][0].size()*2+1, 0);
    int centI = (int)(inpDiff.size()/2);
    int indHR =  holeRat*smearedImg[0][it].size();
    for (int ir=0; ir<(int)smearedImg[0][it].size(); ir++) {
      if (ir < indHR) {
        inpDiff[centI+1+ir] = smearedImg[0][it][indHR]
              *pow(sin((PI/2)*((ir+1)/((double)(indHR+1)))), 2);
        //inpDiff[centI+1+ir] = smearedImg[0][it][ir]
        //      *pow(sin((PI/2)*((ir+1)/((double)(indHR+1)))), 2);
        inpDiff[centI-1-ir] = inpDiff[centI+1+ir];
      }
      else {
        inpDiff[centI+1+ir] = smearedImg[0][it][ir];
        inpDiff[centI-1-ir] = smearedImg[0][it][ir];
      }
    }

    pairCorr1d = tools::fft1dRtoC(inpDiff, rMaxRat, NautCpadding, 
          padDecayRat, fftB, qSpace, rSpace, false);

    // Retrieve result and save results
    for (int ir=0; ir<outSize; ir++) {
      pairCorr[it][ir] = pairCorr1d[0][ir]; 
    }
  }

  // Saving pair correlation
  save::saveDat<double>(pairCorr, "./plots/data/data_pairCorrDiffSmear"
      + to_string(stdev) + "["
      + to_string(pairCorr.size()) + ","
      + to_string(pairCorr[0].size()) + "].dat");
  save::saveDat<double>(pairCorr[pairCorr.size()-1], 
      "./plots/data/data_pairCorrFinalDiffSmear" 
      + to_string(stdev) + "["
      + to_string(pairCorr[0].size()) + "].dat");
 


  /////////////////////////
  /////  Cleaning up  /////
  /////////////////////////

  if (verbose) {
    std::cout << "Cleaning up" << endl;
  }

  fftw_destroy_plan(fftB);
  fftw_free(rSpace);
  fftw_free(qSpace);
  return 1;
}
