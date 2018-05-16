#include "alignRuns.h" 


using namespace std;


int main(int argc, char* argv[]) {

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

  double qMax = 11.3;
  int NradBins = 30;
  std::vector<double> atmDiff(NradBins);
  FILE* atmFile = fopen("../simulation/diffractionPattern/output/references/NBZrefDiff_atmDiffractionPatternLineOut_Bins-30_Qmax-11.300000_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000.dat", "rb");
  fread(&atmDiff[0], sizeof(double), atmDiff.size(), atmFile);
  fclose(atmFile);
 


  /////  Flags  /////
  bool findVarRegions = true;

  
  /////  Setting up variables  /////
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

  std::map<string, map<string, vector<double> > > cbar;
  cbar["20161102"]["LongScan1"].push_back(13);
  cbar["20161102"]["LongScan2"].push_back(25);
  cbar["20161104"]["LongScan2"].push_back(7);
  cbar["20161106"]["LongScan1"].push_back(6);
  cbar["20161102"]["LongScan1"].push_back(9);
  cbar["20161102"]["LongScan2"].push_back(25);
  cbar["20161104"]["LongScan2"].push_back(10);
  cbar["20161106"]["LongScan1"].push_back(10);


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
  double* qSpace = (double*) fftw_malloc(sizeof(double)*FFTsize);
  fftw_complex* rSpace = 
        (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(int)((FFTsize/2) + 1));
  fftw_plan fftB;
  fftB = fftw_plan_dft_r2c_1d(FFTsize, qSpace, rSpace, FFTW_MEASURE);

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

          ind = 1;
          for (auto coeffs : legCoeff_map) {
            cout<<"pos: "<<coeffs.first<<endl;
            legDelays[ind] = 2*(coeffs.first/(3e-1)-0.5);
            ind++;
          }
          legDelays[0] = 2*legDelays[1] - legDelays[2];

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
            cout<<"111"<<endl;

            for (uint iy=0; iy<img[i][0].size(); iy++) {
              double sum = 0;
              for (uint ix=0; ix<img[i].size(); ix++) {
                sum += img[i][ix][iy];
              }
              sum /= (double)img[i].size();
              for (uint ix=0; ix<img[i].size(); ix++) {
                img[i][ix][iy] -= sum;
                img[i][ix][iy] /= 1e20*atmDiff[iy]*(qMax*(iy+0.5)/NradBins);
                if (iy < 3) {
                  img[i][ix][iy] = 0;
                }
              }
            }
            cout<<"222"<<endl;

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

              // Subtract mean
              sum = 0;
              for (uint ix=0; ix<img[i].size(); ix++) {
                sum += smearedImg[i][ix][iy];
              }
              sum /= (double)smearedImg[i].size();
              for (uint ix=0; ix<smearedImg[i].size(); ix++) {
                smearedImg[i][ix][iy] -= sum;
              }
            }
            cout<<"last entry"<<endl;
            cout<<"plotting"<<endl;

            /*
              vals[2] = "-8";
              vals[3] = "8";
              plt.printXY(img[i], legDelays, "plotRun" + to_string(curRun) + "_" + curDate + "_" + curScan 
                    + "_Leg" + to_string(i), oppts, vaals);
            */
              
            if (false && (i==0 || i==2) && (cbar.find(curDate)!=cbar.end()) && (cbar[curDate].find(curScan)!=cbar[curDate].end())) {
              vals[2] = to_string(-1*cbar[curDate][curScan][i/2]);
              vals[3] = to_string(cbar[curDate][curScan][i/2]);
              plt.printRC(img[i], legDelays, 0, qMax, curDate + "_" + curScan 
                    + "_Leg" + to_string(i), opts, vals);
                    //+ "_Leg" + to_string(i) + "_run" + to_string(curRun), opts, vals);
            }
            else {
              vals[2] = "-4";
              vals[3] = "4";
              std::vector<PLOToptions> oopts;
              std::vector<std::string> vvals;
              //plt.printXY(img[i], legDelays, curDate + "_" + curScan 
              //      + "_Leg" + to_string(i), oopts,vvals);
              plt.printRC(img[i], legDelays, 0, qMax, curDate + "_" + curScan 
                    + "_Leg" + to_string(i), opts, vals);
                    //+ "_Leg" + to_string(i) + "_run" + to_string(curRun), oppts, vaals);
              plt.printRC(smearedImg[i], legDelays, 0, qMax, "smeared_" + curDate + "_" + curScan
                    + "_Leg" + to_string(i), opts, vals);
                    //+ "_Leg" + to_string(i) + "_run" + to_string(curRun), oppts, vaals);

            }
            
            cout<<"plotted"<<endl;
          }

          cout<<"start autocorrelation"<<endl;

          double rMaxRat = 0.75;
          double padDecayRat = 0.5;
          double *autDelays = new double[avgImgs_map.size() + 1];
          ind = 1;
          for (auto coeffs : avgImgs_map) {
            autDelays[ind] = 2*(coeffs.first/(3e-1)-0.5);
            ind++;
          }
          autDelays[0] = 2*autDelays[1] - autDelays[2];
          for (uint k=0; k<avgImgs_map.size(); k++) {
                //cout<<"timing: "<<autDelays[k]<<endl;
                if (autDelays[k] > autDelays[k+1]) {
                        cout<<"WRONG!!!!!"<<endl;
                }
          }




          int outSize = FTsize*rMaxRat;
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
                inpDiff[centI-1-ir] = smearedImg[0][it][indHR]
                      *pow(sin((PI/2)*((ir+1)/((double)(indHR+1)))), 2);
              }
              else {
                inpDiff[centI+1+ir] = smearedImg[0][it][ir];
                inpDiff[centI-1-ir] = smearedImg[0][it][ir];
              }
                //inpDiff[centI+1+ir] = smearedImg[0][it][ir];
                //inpDiff[centI-1-ir] = smearedImg[0][it][ir];
            }
            //cout<<"filled"<<endl;
            //plt.print1d(inpDiff, "inpDiff");

            autCorr1d = tools::fft1dRtoC(inpDiff, rMaxRat, NautCpadding, 
                  padDecayRat, fftB, qSpace, rSpace, false);
            //plt.print1d(autCorr1d[0], "autCorReal");
            //plt.print1d(autCorr1d[1], "autCorImg");
            //plt.print1d(autCorr1d[2], "autCorPow");

            for (ir=0; ir<outSize; ir++) {
              autCorr[it][ir] = autCorr1d[2][ir];
            }

          }
 
          std::vector< std::vector<double> > autCplot(smearedImg[0].size());
          for (uint it=0; it<autCorr.size(); it++) {
            autCplot[it].resize(outSize);
            for (int ilo=0; ilo<autCplot[it].size(); ilo++) {
              if (ilo < 3) {
                autCplot[it][ilo] = 0;
              }
              else {
                autCplot[it][ilo] = autCorr[it][ilo];
              }
              autCplot[it][ilo] = autCorr[it][ilo];
            }
          }
          opts.resize(2);
          vals.resize(2);
          opts[0] = maximum;    vals[0] = "4e-1";
          //opts[1] = minimum;    vals[1] = "-5e2";
          //opts[0] = xLabel;     vals[0] = "R [arb]";
          opts[1] = yLabel;     vals[1] = "Time [ps]";
          //opts[2] = draw;       vals[2] = "CONT4Z";
 
          double rMax = (rMaxRat*NradBins)*(2*PI/(2*11.3));
          cout<<"ANGS: "<<2*PI*(double)((int)(NradBins/2))/(2*11.3)<<"   "<<NradBins*(2*PI/(2*11.3))<<"  "<<rMax<<endl;
          cout<<"shapes: "<<autCplot.size()<<" "<<autCplot[0].size()<<"  "<<legCoeff_map.size() + 1<<endl;
          plt.printRC(autCplot, legDelays, 0, rMax, curDate + "_" + curScan + "_AutCor", opts, vals);


          delete[] legDelays;
          delete[] autDelays;
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
      legCoeff_map[curPosition][i] += imgNorm*(*legCoeffs)[i];
    }
    for (uint ir=0; ir<(*imgSubBkg).size(); ir++) {
      for (uint ic=0; ic<(*imgSubBkg)[ir].size(); ic++) {
        avgImgs_map[curPosition][ir][ic] += imgNorm*(*imgSubBkg)[ir][ic];
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
