#include "plotLegendres.h"
#include "../simulation/diffractionPattern/simulations.h"
#include <TLegend.h>


using namespace std;


int main(int argc, char* argv[]) {

  /////  Flags  /////
  bool compareFinalStates = true;



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

  const double maxQ = 9.726264; //11.3;
  const int NradBins = 30;
  const int Nlegendres = 5;
  const int fitRange = NradBins - 5;
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

  double dshift = 4e-3;
  double seed = clock();


  int ir, ic, itr = 0;
  int curRun = -1;
  int Nscans = 0;
  string runInd = "";
  string curDate = "";
  string curScan = "";
  double curPosition = 0;
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
  bool plotSimDiffs = true;

  double Iebeam     = 5;
  double elEnergy   = 3.7e6;
  double screenDist = 4;
  std::string simReferenceDir = "../simulation/diffractionPattern/output/references/";

  int arrSize = 188448;
  bool filledRun = false;
  double* refImg = NULL;
  double prevStagePos, stageDiff;
  double stageCut = 1.8e-2; //0.0041;
  double compare;
  double stdev = 0.75;

  bool newEntry;
  double rad;
  std::map< int32_t, int > stagePosInds;
  std::map< int, std::vector<double> > scanCounts;
  std::map< int, std::vector< std::vector<double> > > scanLgndrs;


  
  
  std::map< int32_t, std::vector<double> > diffPs, legCoeff_map, autCorr_map;
  std::map< int32_t, std::vector< std::vector<double> > > avgImgs_map;
  std::map< int32_t, double > counts;
  std::map< int32_t, int > Nimages;
  std::map< int32_t, std::vector<int32_t> > allPos;

  int NimgRows = (*imgSubBkg).size();
  int NimgCols = (*imgSubBkg)[0].size();
  int NdiffInds = (int)(NimgRows/2);

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
  for (ir=0; ir<NradBins; ir++) {
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
      for (ir=0; ir<NradBins; ir++) {
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
      for (ir=0; ir<NradBins; ir++) {
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

  std::vector<string> runInds;
  std::map< string, std::map< double, double* > >  diffP_arrays;
  std::map< string, int > runShifts;

  cout<<"begin adding to maps"<<endl;
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

    cout<<"ind/sizes: "<<pInd<<"  "<<scanCounts[runNum].size()<<"  "<<scanLgndrs[runNum].size()<<"  "<<scanLgndrs[runNum][pInd].size()<<endl;
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
      norm += scanCounts[sLitr.first][pItr.second];
    }

    // calculate mean
    for (uint i=0; i<runMeans.size(); i++) {
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
  cout<<"Removing "<<removePos.size()<<" bins"<<endl;
  cout<<"Before removal: "<<stagePosInds.size()<<"  "<<scanLgndrs.begin()->second.size()<<"  "<<scanCounts.begin()->second.size()<<"  "<<runMeans.size()<<"  "<<runStdev.size()<<endl;
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
  int ind = 0;
  for (auto& pItr : stagePosInds) {
    pItr.second = ind;
    ind++;
  }
  cout<<"After removal: "<<stagePosInds.size()<<"  "<<scanLgndrs.begin()->second.size()<<"  "<<scanCounts.begin()->second.size()<<"  "<<runMeans.size()<<"  "<<runStdev.size()<<endl;

  ///  Calculate bin standard deviation  ///
  for (auto pItr : stagePosInds) {
    double norm = 0;
    for (auto sLitr : scanLgndrs) {
      // mean
      for (uint i=0; i<(*legCoeffs).size(); i++) {
        runStdev[pItr.second][i] 
            += std::pow((sLitr.second[pItr.second][i] - runMeans[pItr.second][i]), 2);
      }
      norm += scanCounts[sLitr.first][pItr.second];
    }

    cout<<"\n\nTIME: "<<pItr.first<<endl<<endl;
    // calculate standard deviation
    for (uint i=0; i<runMeans.size(); i++) {
      runStdev[pItr.second][i] = std::sqrt(runStdev[pItr.second][i]/norm);
      cout<<runMeans[pItr.second][i]<<"/"<<runStdev[pItr.second][i]<<"\t"<<endl;
    }

  }



  for (uint i=0; i<8; i++) {
    runMeans.erase(--runMeans.end());
    cout<<"reduce size: "<<runMeans.size()<<endl;
    stagePosInds.erase(--stagePosInds.end());
  }



  cout<<"calculating time  "<<runMeans.size()<<endl;
  // Calculating time delays from stage position
  double *timeDelays = new double[runMeans.size() + 1];
  ind = 1;
  for (auto pItr : stagePosInds) {
    timeDelays[ind] = 2*pItr.first/(3e3);
    cout<<"time: "<<timeDelays[ind]<<endl;
    ind++;
  }
  timeDelays[0] = 2*timeDelays[1] - timeDelays[2];
  for (int i=runMeans.size(); i>=0; i--) {
    timeDelays[i] -= timeZero + timeDelays[0];
  }

  save::saveDat<double>(timeDelays, runMeans.size() + 1, 
      "./plots/data/timeDelays["
      + to_string(runMeans.size() + 1) + "].dat");
 
 /////  Filling array of time dependend legendre signal  /////
 std::vector<double> unPumped(NradBins, 0);
 int NrefBins = 0;
 if (subtractT0) {
   while (timeDelays[NrefBins+1] < 0) {
     NrefBins++;
   }
 }
 else {
   NrefBins = runMeans.size();
 }
 cout<<"NrefBins"<<endl;

 cout<<"starting ilg loop"<<endl;
 std::vector< std::vector< std::vector<double> > > img(Nlegendres); //legCoeff_map.size());
 for (uint ilg=0; ilg<Nlegendres; ilg++) {
   img[ilg].resize(runMeans.size());
   for (uint it=0; it<runMeans.size(); it++) {
     img[ilg][it].resize(NradBins, 0);
     for (uint ir=0; ir<NradBins; ir++) {
       img[ilg][it][ir] = runMeans[it][ilg*NradBins + ir];
     }
   }
   cout<<"filled img"<<endl;

   // Mean/t0 subtraction and normalizing by atomic scattering
   // Raw data
   std::fill(unPumped.begin(), unPumped.end(), 0.0);
   for (ir=0; ir<NradBins; ir++) {
     for (uint tm=0; tm<NrefBins; tm++) {
       unPumped[ir] += img[ilg][tm][ir];
     }
     unPumped[ir] /= NrefBins;
   }
   for (ir=0; ir<NradBins; ir++) {
     for (uint tm=0; tm<img[ilg].size(); tm++) {
       img[ilg][tm][ir] -= unPumped[ir];
       img[ilg][tm][ir] *= sMsNorm[ir];
     }
   }
   if (ilg==0) {
     plt.print1d(unPumped, "testUnp");
   }

   ///  Saving data  ///
   save::saveDat<double>(img[ilg], "./plots/data/data_sMsL" 
       + to_string(ilg) + "Diff["
       + to_string(img[ilg].size()) + ","
       + to_string(img[ilg][0].size()) + "].dat");
   save::saveDat<double>(img[ilg][img[ilg].size()-1], "./plots/data/data_sMsFinalL" 
       + to_string(ilg) + "Diff["
       + to_string(img[ilg][0].size()) + "].dat");


   // Save raw unpumped data to fit background
   if (ilg == 0) {
     save::saveDat<double>(unPumped, "./qScale/results/unPumpedDiffractionL0["
         + to_string(unPumped.size()) + "].dat");
     plt.print1d(unPumped, "unPumpedRaw");
   }

 }



  exit(0);






  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);


    /////  Make new root file and initialize variables  /////
    if ((curDate != *date) || (curScan != *scan) || (ievt == Nentries-1)/* || (curRun != runNum)*/) {

      if  ((curDate != *date) || (curScan != *scan) || (ievt == Nentries-1)/* || (curRun != runNum)*/) {

        ///  Prune results  ///
        auto litr = legCoeff_map.begin();
        auto aitr = avgImgs_map.begin();
        auto citr = counts.begin();
        auto iitr = Nimages.begin();
        while (citr != counts.end()) {
          cout<<"counts: "<<citr->first<<"   "<<iitr->second<<"  "<<citr->second<<endl;
          if (citr->second < 10) {
            //litr--;
            legCoeff_map.erase(litr);
            //aitr--;
            avgImgs_map.erase(aitr);
            //citr--;
            counts.erase(citr);
          }
          litr++;
          aitr++;
          citr++;
          iitr++;
        }


        cout<<"start plot"<<endl;
        if ((curRun != -1) && (legCoeff_map.size() > 1)) {
          cout<<"Plotting"<<endl;
          legCoeff_map.erase(legCoeff_map.begin());
          avgImgs_map.erase(avgImgs_map.begin());
          legCoeff_map.erase(--legCoeff_map.end());
          avgImgs_map.erase(--avgImgs_map.end());
          legCoeff_map.erase(--legCoeff_map.end());
          avgImgs_map.erase(--avgImgs_map.end());
          legCoeff_map.erase(--legCoeff_map.end());
          avgImgs_map.erase(--avgImgs_map.end());
          legCoeff_map.erase(--legCoeff_map.end());
          avgImgs_map.erase(--avgImgs_map.end());
          legCoeff_map.erase(--legCoeff_map.end());
          avgImgs_map.erase(--avgImgs_map.end());
          legCoeff_map.erase(--legCoeff_map.end());
          avgImgs_map.erase(--avgImgs_map.end());
          legCoeff_map.erase(--legCoeff_map.end());
          avgImgs_map.erase(--avgImgs_map.end());

          int ind = 0;
   //       double *timeDelays = new double[legCoeff_map.size() + 1];
     //     std::vector< std::vector< std::vector<double> > > img(Nlegendres); //legCoeff_map.size());
          std::vector< std::vector< std::vector<double> > > smearedImg(Nlegendres); //legCoeff_map.size());

          // Calculating time delays from stage position
          ind = 1;
          for (auto coeffs : legCoeff_map) {
            timeDelays[ind] = 2*coeffs.first/(3e3) - timeZero;
            cout<<"Times: "<<ind<<"  "<<timeDelays[ind]<<"  "<<coeffs.first<<endl;
            ind++;
          }
          timeDelays[0] = 2*timeDelays[1] - timeDelays[2];

          save::saveDat<double>(timeDelays, legCoeff_map.size() + 1, 
              "./plots/data/timeDelays["
              + to_string(legCoeff_map.size() + 1) + "].dat");
              
            opts.resize(2);
            vals.resize(2);
            opts[0] = yLabel;   vals[0] = "Time [ps]";
            opts[1] = xLabel;   vals[1] = "Scattering Q [arb units]";

          for (uint k=0; k<legCoeff_map.size() + 1; k++) {
            cout<<"Time "<<k<<"   "<<timeDelays[k]<<endl;
          }


          /////  Filling array of time dependend legendre signal  /////
          std::vector<double> unPumped(NradBins, 0);
          int NrefBins = 0;
          if (subtractT0) {
            while (timeDelays[NrefBins+1] < 0) {
              NrefBins++;
            }
          }
          else {
            NrefBins = legCoeff_map.size();
          }

          for (uint ilg=0; ilg<Nlegendres; ilg++) {
            ind = 0;
            img[ilg].resize(legCoeff_map.size());
            for (auto coeffs : legCoeff_map) {
              img[ilg][ind].resize(NradBins, 0);
              for (uint ir=0; ir<NradBins; ir++) {
                img[ilg][ind][ir] = coeffs.second[ilg*NradBins + ir]/counts[coeffs.first];
              }
              ind++;
            }

            // Mean/t0 subtraction and normalizing by atomic scattering
            // Raw data
            std::fill(unPumped.begin(), unPumped.end(), 0.0);
            for (ir=0; ir<NradBins; ir++) {
              for (uint tm=0; tm<NrefBins; tm++) {
                unPumped[ir] += img[ilg][tm][ir];
              }
              unPumped[ir] /= NrefBins;
            }
            for (ir=0; ir<NradBins; ir++) {
              for (uint tm=0; tm<img[ilg].size(); tm++) {
                img[ilg][tm][ir] -= unPumped[ir];
                img[ilg][tm][ir] *= sMsNorm[ir];
              }
            }

            ///  Saving data  ///
            save::saveDat<double>(img[ilg], "./plots/data/data_sMsL" 
                + to_string(ilg) + "Diff["
                + to_string(img[ilg].size()) + ","
                + to_string(img[ilg][0].size()) + "].dat");
            save::saveDat<double>(img[ilg][img[ilg].size()-1], "./plots/data/data_sMsFinalL" 
                + to_string(ilg) + "Diff["
                + to_string(img[ilg][0].size()) + "].dat");


            // Save raw unpumped data to fit background
            if (ilg == 0) {
              save::saveDat<double>(unPumped, "./qScale/results/unPumpedDiffractionL0["
                  + to_string(unPumped.size()) + "].dat");
              plt.print1d(unPumped, "unPumpedRaw");
            }


            /*
            if (i == 0) {
            Eigen::Matrix<double, fitRange, 3> X;
            Eigen::Matrix<double, fitRange, 1> Y;
            int shift = NradBins-fitRange;
            for (ir=0; ir<fitRange; ir++) {
              X(ir,0) = 1;
              X(ir,1) = atmDiff[shift+ir];
              X(ir,2) = molDiff[shift+ir];
              Y(ir,0) = unPumped[shift+ir];
            }

            //Eigen::MatrixXd nrmInp = X.transpose()*X;
            //Eigen::MatrixXd norm = tools::SVDinvert(nrmInp);
            //Eigen::MatrixXd weights = norm*(X.transpose()*Y);
            Eigen::MatrixXd weights = tools::normalEquation(X, Y);

            for (uint ir=0; ir<NradBins; ir++) {
              unPumped[ir] -= weights(0) + weights(1)*atmDiff[ir];
            }
 
              save::saveDat<double>(unPumped, "./plots/data/unPumpedAvg[" 
                  + to_string(NradBins) + "].dat");
              plt.print1d(unPumped, "unPumped");
              plt.printRC(img[i], timeDelays, 0, maxQ, "testNorm", opts, vals);
            }
            */


           cout<<"222"<<endl;

            ///  Smearing time dependence  ///
            double sum = 0;
            smearedImg[ilg].resize(legCoeff_map.size());
            for (uint j=0; j<legCoeff_map.size(); j++) {
              smearedImg[ilg][j].resize(NradBins, 0);
            }
            for (uint ir=0; ir<img[ilg][0].size(); ir++) {
              for (uint tm=0; tm<img[ilg].size(); tm++) {
                sum = 0;
                //cout<<"delay: "<<timeDelays[ix]<<endl;
                for (uint tmm=0; tmm<img[ilg].size(); tmm++) {
                  smearedImg[ilg][tm][ir] += img[ilg][tmm][ir]
                                            *std::exp(-1*std::pow(timeDelays[tm] - timeDelays[tmm], 2)/(2*std::pow(stdev, 2)));
                  sum += exp(-1*std::pow(timeDelays[tm] - timeDelays[tmm], 2)/(2*std::pow(stdev, 2)));
                }
                smearedImg[ilg][tm][ir] /= sum;
              }
            }

            ///  Subtract before t0 and normalize to sMs  ///
            std::fill(unPumped.begin(), unPumped.end(), 0.0);
            for (ir=0; ir<NradBins; ir++) {
              for (uint tm=0; tm<NrefBins; tm++) {
                unPumped[ir] += smearedImg[ilg][tm][ir];
              }
              unPumped[ir] /= (double)NrefBins;
            }
            for (ir=0; ir<NradBins; ir++) {
              for (uint tm=0; tm<img[ilg].size(); tm++) {
                smearedImg[ilg][tm][ir] -= unPumped[ir];
              }
            }

            if (ilg == 0) {
              for (uint tm=0; tm<NrefBins; tm++) {
                plt.print1d(smearedImg[ilg][tm], "testsmeart0_"+to_string(tm));
              }
              save::saveDat<double>(unPumped, "./plots/data/data_unPumpedDiffractionL0Smeared["
                  + to_string(NradBins) + "].dat");
              plt.print1d(unPumped, "unPumpedSmeared");
              cout<<"Unpump Smear: "<<unPumped.size()<<endl;
              for (uint k=0; k<unPumped.size(); k++) {
                cout<<"smear  "<<k<<"  "<<unPumped[k]<<endl;
              }
            }
            if (ilg==0 || ilg==2) {
              plt.printRC(smearedImg[ilg], timeDelays, 0, maxQ, "testSmearRaw" + to_string(ilg), opts, vals);
            }

            save::saveDat<double>(smearedImg[ilg], "./plots/data/data_sMsL"
                + to_string(ilg) + "DiffSmear"
                + to_string(stdev) + "["
                + to_string(smearedImg[ilg].size()) + ","
                + to_string(smearedImg[ilg][0].size()) + "].dat");
            
            save::saveDat<double>(smearedImg[ilg][smearedImg[ilg].size()-1], 
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

            plt.printRC(img[ilg], timeDelays, 0, maxQ, curDate + "_" + curScan 
                  + "_Leg" + to_string(ilg), opts, vals);
                  //+ "_Leg" + to_string(ilg) + "_run" + to_string(curRun), oppts, vaals);
            plt.printRC(smearedImg[ilg], timeDelays, 0, maxQ, "smeared_" + curDate + "_" + curScan
                  + "_Leg" + to_string(ilg), opts, vals);
                  //+ "_Leg" + to_string(ilg) + "_run" + to_string(curRun), oppts, vaals);

            
            cout<<"plotted"<<endl;
          }

          cout<<"start autocorrelation"<<endl;

          /////  Pair correlation function  /////
          std::vector< std::vector<double> > autCorr(smearedImg[0].size()), autCorr1d;
          std::vector<double> powerSpct(outSize);
          for (uint it=0; it<smearedImg[0].size(); it++) {
            autCorr[it].resize(outSize);

            std::vector<double> inpDiff(smearedImg[0][0].size()*2+1, 0);
            int centI = (int)(inpDiff.size()/2);
            int indHR =  holeRat*smearedImg[0][it].size();
            for (ir=0; ir<smearedImg[0][it].size(); ir++) {
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
            //cout<<"filled"<<endl;
            //plt.print1d(inpDiff, "inpDiff");

            autCorr1d = tools::fft1dRtoC(inpDiff, rMaxRat, NautCpadding, 
                  padDecayRat, fftB, qSpace, rSpace, false);

            // Retrieve result and save results
            for (ir=0; ir<outSize; ir++) {
              autCorr[it][ir] = autCorr1d[0][ir]; 
              //if (ir < 3) {
              //  autCorr[it][ir] = 0;
              //}
            }
          }
 
          // Saving pair correlation
          save::saveDat<double>(autCorr, "./plots/data/data_pairCorrDiffSmear"
              + to_string(stdev) + "["
              + to_string(autCorr.size()) + ","
              + to_string(autCorr[0].size()) + "].dat");
          save::saveDat<double>(autCorr[autCorr.size()-1], 
              "./plots/data/data_pairCorrFinalDiffSmear" 
              + to_string(stdev) + "["
              + to_string(autCorr[0].size()) + "].dat");
          /*
          for (auto& pInd : autCorrLOinds) {
            cout<<"IND: "<<pInd<<"   "<<timeDelays[pInd]<<endl;
            save::saveDat<double>(autCorr[pInd], "./plots/data/data_pairCorrLO_smear"
                  + to_string(stdev) + "_time"
                  + to_string(timeDelays[pInd]) + "_["
                  + to_string(autCorr[pInd].size()) + "].dat");
          }
          */


          /*
          cout<<"start opts/vals"<<endl;
          opts.resize(2);
          vals.resize(2);
          opts[0] = yLabel;     vals[0] = "Time [ps]";
          opts[1] = xLabel;     vals[1] = "R [Ang]";
          //opts[2] = draw;       vals[2] = "CONT4Z";
 
          cout<<"if opts/vals"<<endl;
          std::string cbarAutName = curDate + "_" + curScan + "_Aut"; 
          if (cbar.find(cbarAutName) != cbar.end()) {
            opts.push_back(minimum);
            vals.push_back(to_string(-1*cbar[cbarAutName]));
            opts.push_back(maximum);
            vals.push_back(to_string(cbar[cbarAutName]));
          }

          cout<<"plotting "<<opts.size()<<"   "<<vals.size()<<endl;
          plt.printRC(autCorr, timeDelays, 0, rMax, curDate + "_" + curScan + "_AutCor", opts, vals);

          opts.resize(2);
          vals.resize(2);
          opts[0] = xSpan;    vals[0] = "0," + to_string(rMax);
          opts[1] = xLabel;   vals[1] = "R [Angs]";

          cout<<"Time Size: "<<autCorr.size()<<endl;
          std::vector<TH1F*> hists;
          for (auto& pInd : autCorrLOinds) {
            cout<<"size: "<<autCorr[pInd].size()<<endl;
            hists.push_back(plt.plot1d(autCorr[pInd], "autCorLO" + to_string(pInd), opts, vals));
          }
          TCanvas* Canv = new TCanvas("Canv","Canv", 800, 600);
          //TLegend* lgnd = new TLegend(0.6, 0.15, 0.8, 0.55);
          TLegend* lgnd = new TLegend(0.6, 0.6, 0.8, 0.95);
          hists[0]->SetMinimum(-1e-3);
          hists[0]->SetMaximum(1e-3);
          hists[0]->Draw("l");
          lgnd->AddEntry(hists[0], to_string(timeDelays[autCorrLOinds[0]]).c_str(), "l");
          for (uint i=1; i<hists.size(); i++) {
            hists[i]->SetLineColor(i+1);
            hists[i]->Draw("lSAME");
            lgnd->AddEntry(hists[i], to_string(timeDelays[autCorrLOinds[i]]).c_str(), "l");
          }
          lgnd->Draw("SAME");
          Canv->Print("autCorrDelays.png");
          cout<<"plotted corrdels"<<endl;

          cout<<"deletin 1"<<endl;
          for (auto& hst : hists) {
            delete hst;
          }
          cout<<"del 2"<<endl;
          delete Canv;
          cout<<"del 3"<<endl;

          */
          unPumped.resize(autCorr[0].size(), 0);
          std::fill(unPumped.begin(), unPumped.end(), 0.0);
          for (ir=0; ir<NradBins; ir++) {
            for (uint tm=0; tm<NrefBins; tm++) {
              unPumped[ir] += autCorr[tm][ir];
            }
            unPumped[ir] /= (double)NrefBins;
          }
          save::saveDat<double>(unPumped, "./plots/data/data_unPumpedPairCorr["
              + to_string(autCorr[0].size()) + "].dat");
          plt.print1d(unPumped, "unPumpedSmearedAutCorr");


          delete[] timeDelays;
          cout<<"del 4"<<endl;
          cout<<"del 5"<<endl;
        }
        
        cout<<"clearing"<<endl;
        legCoeff_map.clear();
        counts.clear();
        avgImgs_map.clear();
        autCorr_map.clear();
      }

          cout<<"reset"<<endl;
      filledRun = false;
      curPosition = 0;
      curRun  = runNum;
      curDate = (*date);
      curScan = (*scan);

      runInd = curDate + "_" + curScan + "_" + to_string(curRun);
      runInds.push_back(runInd);
      cout<<runInd<<endl;
      runShifts[runInd] = 0;
      prevStagePos = 1e10;
      Nscans = 0;


    }

    // Don't fill run with large time steps
    if (filledRun) {
      //continue;
    }

    cout<<"StagePos"<<stagePos<<endl;
    // Need initial stage position (skips first image)
    if (prevStagePos > stagePos) {
      Nscans++;
      prevStagePos = stagePos;
      curPosition = 0;
    }

    // Round the stage difference to the correct decimal place
    /*
    compare = 1;
    itr = -1;
    while ((fabs(compare) > 0.01) && (itr < 4)) {
      itr++;
      compare = (tools::round((stagePos - prevStagePos)*pow(10, itr))/pow(10, itr) 
          - (stagePos - prevStagePos))
          /(stagePos - prevStagePos);
    }
    */
    //stageDiff = ((double)tools::round((stagePos - prevStagePos)*pow(10, itr)))/pow(10, itr);
    stageDiff = stagePos - prevStagePos;
    prevStagePos = stagePos;

    //stageDiff = ((double)((int)(1000*stageDiff)))/1000;

    curPosition += stageDiff;
    //cout<<"stage diff: "<<prevStagePos<<"  "<<stagePos<<"  "<<stageDiff<<"      "<<curPosition<<endl;
    //cout<<endl;

    // Make sure we are still in regime of small step sizes
    if (stageDiff > stageCut) {
      filledRun = true;
      //continue;
    }

    //plt.print1d(*legCoeffs, "testLC"+to_string(runNum)+"_"+to_string(ievt));

    //curPosition = ((double)((int)(1000*curPosition)))/1000.;
    //if (legCoeff_map.find(curPosition) == legCoeff_map.end()) {
    newEntry = true;
    for (auto itr = legCoeff_map.begin(); itr!=legCoeff_map.end(); itr++) {
      if (curPosition == itr->first) {
        curPosition = itr->first;
        newEntry = false;
      }
    }
    if (newEntry) {
      if (Nscans > 3) {
        cout << "WARNING: Skipping position " << curPosition << endl;
        continue;
      }

      cout<<"ADDED "<<curPosition<<endl;
      legCoeff_map[curPosition].resize((*legCoeffs).size(), 0);
      avgImgs_map[curPosition].resize(NimgRows);
      for (int ir=0; ir<NimgRows; ir++) {
        avgImgs_map[curPosition][ir].resize(NimgCols);
      }
      counts[curPosition] = 0;
      Nimages[curPosition] = 0;
    }
    for (uint i=0; i<(*legCoeffs).size(); i++) {
      //legCoeff_map[curPosition][i] += imgNorm*(*legCoeffs)[i]; CHANGED
      legCoeff_map[curPosition][i] += (*legCoeffs)[i];
    }
    for (uint ir=0; ir<(*imgSubBkg).size(); ir++) {
      for (uint ic=0; ic<(*imgSubBkg)[ir].size(); ic++) {
        //avgImgs_map[curPosition][ir][ic] += imgNorm*(*imgSubBkg)[ir][ic]; CHANGED
        avgImgs_map[curPosition][ir][ic] += (*imgSubBkg)[ir][ic];
      }
    }

    counts[curPosition] += imgNorm;
    Nimages[curPosition] += 1;
    allPos[curPosition].push_back(stagePos);
  }

  // Check stage position binning
  for (auto itr : allPos) {
    for (uint k=1; k<itr.second.size(); k++) {
      if (itr.second[k-1] != itr.second[k]) {
        cerr << "WARNING: Stage positions "
          + to_string(itr.second[k-1]) + " and "
          + to_string(itr.second[k]) + " do no match in bin "
          + to_string(itr.first) << endl;
      }
    }
  }
  // Clean up

  for (auto &runItr : diffP_arrays) {
    for (auto &imgItr : runItr.second) {
      delete[] imgItr.second;
    }
  }

  fftw_destroy_plan(fftB);
  fftw_free(rSpace);
  fftw_free(qSpace);
  return 1;
}
