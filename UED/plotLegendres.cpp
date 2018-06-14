#include "alignRuns.h"
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

  double maxQ = 9.726264; //11.3;
  int NradBins = 30;
  std::vector<double> atmDiff(NradBins);
 
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

  double dshift = 4e-3;
  double seed = clock();


  ////////////////////////////////////////////
  /////  Finding Regions that Fluctuate  /////
  ////////////////////////////////////////////

  int ir, ic, itr = 0;
  int curRun = -1;
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
  oppts[2] = draw;     vaals[2] = "CONT4Z";
  opts[0] = yLabel;   vals[0] = "Time [ps]";
  opts[1] = xLabel;   vals[1] = "Scattering Q [arb units]";
  opts[4] = draw;     vals[4] = "CONT4Z";
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
  double stdev = 1.0;

  bool newEntry;
  double rad;
  std::map< double, std::vector<double> > diffPs, legCoeff_map, autCorr_map;
  std::map< double, std::vector< std::vector<double> > > avgImgs_map;
  std::map< double, double > counts;

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

  double* qSpace = (double*) fftw_malloc(sizeof(double)*FFTsize);
  fftw_complex* rSpace = 
        (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(int)((FFTsize/2) + 1));
  fftw_plan fftB;
  fftB = fftw_plan_dft_r2c_1d(FFTsize, qSpace, rSpace, FFTW_MEASURE);


  //int NrefSub = 30;
  int NrefSub = 20;
  ////////////////////////////////////////////
  /////  Importing simulated references  /////
  ////////////////////////////////////////////

  std::string fileNameSuffix = "_Bins-" + to_string(NradBins) + "_Qmax-" + to_string(maxQ)
                                + "_Ieb-" + to_string(Iebeam) + "_scrnD-"
                                + to_string(screenDist) + "_elE-" + to_string(elEnergy) + ".dat";

  ///  Nitrobenzene reference  ///
  cout << "fileName: " << simReferenceDir + "/nitrobenzene_atmDiffractionPatternLineOut" + fileNameSuffix<<endl;
  FILE* atmFile = fopen((simReferenceDir + "/nitrobenzene_atmDiffractionPatternLineOut"
                           + fileNameSuffix).c_str(), "rb");
  fread(&atmDiff[0], sizeof(double), atmDiff.size(), atmFile);
  fclose(atmFile);
  plt.print1d(atmDiff, "testATMNBZ");

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

      delete plt.print1d(refLineOut[name], name + "_diffractionSimOrig");
      // Normalize by atomic scattering
      for (ir=0; ir<NradBins; ir++) {
        refLineOut[name][ir] /= 1e20*atmDiff[ir]*(maxQ*(ir+0.5)/NradBins);
      }
      save::saveDat<double>(refLineOut[name], "plots/data/" + name + "SimGrnd["
            + to_string(NradBins) + "].dat");
      delete plt.print1d(refLineOut[name], name + "_diffractionSim");

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
      save::saveDat<double>(refAutCor[name], name + "AutCorr["
              + to_string(refAutCor[name].size()) + "].dat");
      cout<<"autcorred"<<endl;

    }

    if (plotSimDiffs) {
      opts.resize(2);
      vals.resize(2);
      //opts[0] = draw;     vals[0] = "l";
      opts[0] = xSpan;    vals[0] = "0," + to_string(rMax);
      opts[1] = xLabel;   vals[1] = "R [Angs]";

      std::vector<double> autDiff(refAutCor[radicalNames[0]].size(), 0);
      for (uint i=1; i<Nradicals; i++) {
        for (uint ir=0; ir<autDiff.size(); ir++) {
          autDiff[ir] = refAutCor[radicalNames[i]][ir]
                        - refAutCor[radicalNames[0]][ir];
        }
        delete plt.print1d(autDiff, "autCorDiff_" + radicalNames[i], opts, vals);
      }
    }
  }



  ///////////////////////////
  /////  Aligning Runs  /////
  ///////////////////////////

  std::vector<string> runInds;
  std::map< string, std::map< double, double* > >  diffP_arrays;
  std::map< string, int > runShifts;

  ///// Loop through events in the file /////
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);

    /*
    if (ievt != 11) {
      continue;
    }
    else {
      cout<<"image "<<imgNum<<endl;
      cout<<"size "<<(*legCoeffs).size()<<endl;
      for (int i=0; i<(*legCoeffs).size(); i++) {
          cout<<(*legCoeffs)[i]<<endl;
      }
    }
    */
    //if ((*date != "20161105")) {
    //if ((*date != "20161104") || (*scan != "LongScan3")) {
    //  continue;
    //}
    //cout << imgNum << "  " << runNum << "  "<< ievt << " "<<Nentries<<endl;
    // Expect the reference image to be the first image
    if (imgIsRef) {
      continue;
    }

    // Ignore reference images taken before scan
    if (stagePos < (t0StagePos - 0.021)) {
      continue;
    }

    /////  Make new root file and initialize variables  /////
    if ((curDate != *date) || (curScan != *scan) || (ievt == Nentries-1)/* || (curRun != runNum)*/) {

      if  ((curDate != *date) || (curScan != *scan) || (ievt == Nentries-1)/* || (curRun != runNum)*/) {
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
          for (auto itr = legCoeff_map.begin(); itr!=legCoeff_map.end(); itr++) {
            cout<<"test pos: "<<itr->first<<endl;
          }
          int ind = 0;
          double *legDelays = new double[legCoeff_map.size() + 1];
          std::vector< std::vector< std::vector<double> > > img(3); //legCoeff_map.size());
          std::vector< std::vector< std::vector<double> > > smearedImg(3); //legCoeff_map.size());

          // Calculating time delays from stage position
          ind = 1;
          for (auto coeffs : legCoeff_map) {
            legDelays[ind] = 2*(coeffs.first/(3e-1)-0.5);
            ind++;
          }
          legDelays[0] = 2*legDelays[1] - legDelays[2];
            opts.resize(2);
            vals.resize(2);
            opts[0] = yLabel;   vals[0] = "Time [ps]";
            opts[1] = xLabel;   vals[1] = "Scattering Q [arb units]";


          /////  Filling array of time dependend legendre signal  /////
          std::vector<double> unPumped(NradBins, 0);
          std::vector<int> autCorrLOinds = {2,45,55,60,65,70};
          for (uint i=0; i<3; i++) {
            cout<<"i: "<<i<<endl;
            ind = 0;
            img[i].resize(legCoeff_map.size());
            for (auto coeffs : legCoeff_map) {
              img[i][ind].resize(NradBins, 0);
              if (i==2) cout<<"counts: "<<coeffs.first<<"   "<<counts[coeffs.first]<<endl;
              for (uint ir=0; ir<NradBins; ir++) {
                img[i][ind][ir] = coeffs.second[i*NradBins + ir]/counts[coeffs.first];
              }
              ind++;
            }
            if (i==0) {
              plt.printRC(img[i], legDelays, 0, maxQ, "testRaw", opts, vals);
            }
            cout<<"111"<<endl;

            // Mean subtraction and normalizing by atomic scattering
            for (uint iy=0; iy<img[i][0].size(); iy++) {
              double sum = 0;
              for (uint ix=0; ix<img[i].size(); ix++) {
                sum += img[i][ix][iy];
              }
              sum /= (double)img[i].size();
              for (uint ix=0; ix<img[i].size(); ix++) {
                img[i][ix][iy] -= sum;
                img[i][ix][iy] /= 1e20*atmDiff[iy]*(maxQ*(iy+0.5)/NradBins);
                if (iy < 3) {
                  img[i][ix][iy] = 0;
                }
              }
            }

           cout<<"222"<<endl;

            ///  Smearing time dependence  ///
            double sum = 0;
            smearedImg[i].resize(legCoeff_map.size());
            for (uint j=0; j<legCoeff_map.size(); j++) {
              smearedImg[i][j].resize(NradBins, 0);
            }
            for (uint iy=0; iy<img[i][0].size(); iy++) {
              for (uint ix=0; ix<img[i].size(); ix++) {
                sum = 0;
                //cout<<"delay: "<<legDelays[ix]<<endl;
                for (uint ixx=0; ixx<img[i].size(); ixx++) {
                  smearedImg[i][ix][iy] += img[i][ixx][iy]
                                            *std::exp(-1*std::pow(legDelays[ix] - legDelays[ixx], 2)/(2*std::pow(stdev, 2)));
                  sum += exp(-1*std::pow(legDelays[ix] - legDelays[ixx], 2)/(2*std::pow(stdev, 2)));
                }
                smearedImg[i][ix][iy] /= sum;
              }

              // Mean subtraction
              sum = 0;
              for (uint ix=0; ix<img[i].size(); ix++) {
                sum += smearedImg[i][ix][iy];
              }
              sum /= (double)smearedImg[i].size();
              for (uint ix=0; ix<smearedImg[i].size(); ix++) {
                smearedImg[i][ix][iy] -= sum;
              }
            }

            ///  Subtract before t0  ///
            // Raw data
            unPumped.resize(NradBins,0);
            for (ir=0; ir<NradBins; ir++) {
              for (uint tm=0; tm<NrefSub; tm++) {
                unPumped[ir] += img[i][tm][ir];
              }
              unPumped[ir] /= NrefSub;
            }
            for (ir=0; ir<NradBins; ir++) {
              for (uint tm=0; tm<img[i].size(); tm++) {
                img[i][tm][ir] -= unPumped[ir];
              }
            }

            if (i == 0) {
              save::saveDat<double>(unPumped, "./plots/data/unPumpedAvg[" 
                  + to_string(NradBins) + "].dat");
              plt.print1d(unPumped, "unPumped");
              plt.printRC(img[i], legDelays, 0, maxQ, "testNorm", opts, vals);
            }

 
            // Smeared data
            unPumped.resize(NradBins,0);
            for (ir=0; ir<NradBins; ir++) {
              for (uint tm=0; tm<NrefSub; tm++) {
                unPumped[ir] += smearedImg[i][tm][ir];
              }
              unPumped[ir] /= NrefSub;
            }
            for (ir=0; ir<NradBins; ir++) {
              for (uint tm=0; tm<img[i].size(); tm++) {
                smearedImg[i][tm][ir] -= unPumped[ir];
              }
            }

            if (i == 0) {
              save::saveDat<double>(unPumped, "./plots/data/unPumpedSmearedAvg["
                  + to_string(NradBins) + "].dat");
              plt.print1d(unPumped, "unPumpedSmeared");
            }
            if (i==0 || i==2) {
              plt.printRC(smearedImg[i], legDelays, 0, maxQ, "testSmearRaw" + to_string(i), opts, vals);
            }
            cout<<"last entry"<<endl;
            cout<<"plotting"<<endl;

            /*
              vals[2] = "-8";
              vals[3] = "8";
              plt.printXY(img[i], legDelays, "plotRun" + to_string(curRun) + "_" + curDate + "_" + curScan 
                    + "_Leg" + to_string(i), oppts, vaals);
            */
              
            ///  Plotting time dependend legendre signal  ///
            opts.resize(3);
            vals.resize(3);
            opts[0] = yLabel;   vals[0] = "Time [ps]";
            opts[1] = xLabel;   vals[1] = "Scattering Q [arb units]";
            opts[2] = draw;     vals[2] = "CONT4Z";


            std::string cbarLegName = curDate + "_" + curScan + "_Leg" + to_string(i); 
            if (cbar.find(cbarLegName) != cbar.end()) {
              opts.push_back(minimum);
              vals.push_back(to_string(-1*cbar[cbarLegName]));
              opts.push_back(maximum);
              vals.push_back(to_string(cbar[cbarLegName]));
            }

            plt.printRC(img[i], legDelays, 0, maxQ, curDate + "_" + curScan 
                  + "_Leg" + to_string(i), opts, vals);
                  //+ "_Leg" + to_string(i) + "_run" + to_string(curRun), oppts, vaals);
            plt.printRC(smearedImg[i], legDelays, 0, maxQ, "smeared_" + curDate + "_" + curScan
                  + "_Leg" + to_string(i), opts, vals);
                  //+ "_Leg" + to_string(i) + "_run" + to_string(curRun), oppts, vaals);

            
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
              if (ir < 3) {
                autCorr[it][ir] = 0;
              }
              else {
                autCorr[it][ir] = autCorr1d[0][ir]; 
              }
            }
          }
 
          // Saving pair correlation line outs
          for (auto& pInd : autCorrLOinds) {
            save::saveDat<double>(autCorr[pInd], "data_pairCorrLO_smear"
                  + to_string(stdev) + "_time"
                  + to_string(legDelays[autCorrLOinds[0]]) + "_["
                  + to_string(autCorr[pInd].size()) + "].dat");
          }


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
          plt.printRC(autCorr, legDelays, 0, rMax, curDate + "_" + curScan + "_AutCor", opts, vals);

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
          lgnd->AddEntry(hists[0], to_string(legDelays[autCorrLOinds[0]]).c_str(), "l");
          for (uint i=1; i<hists.size(); i++) {
            hists[i]->SetLineColor(i+1);
            hists[i]->Draw("lSAME");
            lgnd->AddEntry(hists[i], to_string(legDelays[autCorrLOinds[i]]).c_str(), "l");
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



          unPumped.resize(autCorr[0].size(),0);
          for (ir=0; ir<NradBins; ir++) {
            for (uint tm=0; tm<NrefSub; tm++) {
              unPumped[ir] += autCorr[tm][ir];
            }
            unPumped[ir] /= NrefSub;
          }
          save::saveDat<double>(unPumped, "./plots/data/unPumpedAutCorr["
              + to_string(autCorr[0].size()) + "].dat");
          plt.print1d(unPumped, "unPumpedSmeared");


          delete[] legDelays;
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
      prevStagePos = -1;


    }

    // Don't fill run with large time steps
    if (filledRun) {
      //continue;
    }

    // Need initial stage position (skips first image)
    if (prevStagePos == -1) {
      prevStagePos = stagePos;
      continue;
    }

    // Round the stage difference to the correct decimal place
    compare = 1;
    itr = -1;
    while ((fabs(compare) > 0.01) && (itr < 4)) {
      itr++;
      compare = (tools::round((stagePos - prevStagePos)*pow(10, itr))/pow(10, itr) 
          - (stagePos - prevStagePos))
          /(stagePos - prevStagePos);
    }
    stageDiff = ((double)tools::round((stagePos - prevStagePos)*pow(10, itr)))/pow(10, itr);
    //cout<<"stage diff: "<<stagePos<<"  "<<stagePos - prevStagePos<<"  "<<stageDiff;
    prevStagePos = stagePos;

    stageDiff = ((double)((int)(1000*stageDiff)))/1000;

    curPosition += stageDiff;
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
      //cout<<"adding : "<<curPosition<<"   "<<endl;
          for (auto itr = legCoeff_map.begin(); itr!=legCoeff_map.end(); itr++) {
            //cout<<"search pos: "<<itr->first<<endl;
            if (fabs((curPosition - itr->first)/curPosition) < 0.01) {
              curPosition = itr->first;
              newEntry = false;
            }
          }
    if (newEntry) {
      //cout<<"ADDED"<<endl;
      legCoeff_map[curPosition].resize((*legCoeffs).size(), 0);
      avgImgs_map[curPosition].resize(NimgRows);
      for (int ir=0; ir<NimgRows; ir++) {
        avgImgs_map[curPosition][ir].resize(NimgCols);
      }
      counts[curPosition] = 0;
    }
    for (uint i=0; i<(*legCoeffs).size(); i++) {
      //legCoeff_map[curPosition][i] += imgNorm*(*legCoeffs)[i]; CHANGED
      legCoeff_map[curPosition][i] += imgNorm*(*legCoeffs)[i];
    }
    for (uint ir=0; ir<(*imgSubBkg).size(); ir++) {
      for (uint ic=0; ic<(*imgSubBkg)[ir].size(); ic++) {
        //avgImgs_map[curPosition][ir][ic] += imgNorm*(*imgSubBkg)[ir][ic]; CHANGED
        avgImgs_map[curPosition][ir][ic] += (*imgSubBkg)[ir][ic];
      }
    }

    counts[curPosition] += imgNorm;
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
