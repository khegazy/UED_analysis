#include "alignment.h"
//#include "../simulation/diffractionPattern/simulations.h"
#include <TLegend.h>

using namespace std;


int main(int argc, char* argv[]) {

  if (argc<2) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe 'fileList.txt' !!!"<<endl;
    cerr<<"         Can also run ./analysis.exe 'fileList.txt' 'treeName'"<<endl;
    exit(0);
  }


  ///////////////////////////////////////////////////////////
  /////  Load environment and get the number of events  /////
  ///////////////////////////////////////////////////////////
  
  string fileList(argv[1]);
  string treeName("physics");
  if (argc==3) string treeName(argv[2]);
  analysisClass analysis(fileList, treeName);  // Use this when specifying treeName

  uint64_t Nentries;
  Nentries = analysis.setupEnvironment();  // Alter this function as needed for specific setup

  cout.setf(ios::scientific);
  
  //////////////////////////////////
  /////  Setting up variables  /////
  //////////////////////////////////

  auto iPos = fileList.find("RUN-") + 4;
  auto fPos = fileList.find(".txt");
  std::string runName = fileList.substr(iPos, fPos - iPos);
  cout<<"RunName: "<<runName<<endl;
  alignClass align(runName);
  //parameterClass params(runName);

  /////  Importing variables from command line  /////
  for (int iarg=2; iarg<argc; iarg+=2) {
    if (strcmp(argv[iarg],"-Odir")==0) {
      string str(argv[iarg+1]); 
      align.alignScansOutputDir = str;
    }
    else {
      cerr<<"ERROR!!! Option "<<argv[iarg]<<" does not exist!"<<endl;
      exit(0);
    }
  }

  cout<<"here1"<<endl;
  std::vector<double> reference;
  std::vector<double> referenceCount;

  int FFTsize = align.NradLegBins*2 + align.NautCpadding + 1;
  //int FFTsize = (NdiffInds + NautCpadding)*2 + 1;
  int FTsize = (int)(FFTsize/2) + 1;
  int indHR =  align.holeRat*align.NradLegBins;
  int outSize = FTsize*align.rMaxRat;


  double* qSpace = (double*) fftw_malloc(sizeof(double)*FFTsize);
  fftw_complex* rSpace =
        (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(int)((FFTsize/2) + 1));
  fftw_plan fftB;
  fftB = fftw_plan_dft_r2c_1d(FFTsize, qSpace, rSpace, FFTW_MEASURE);

  cout<<"here2"<<endl;
  /////  Plotting variables  /////
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

  int curScan;
  
  ////////////////////////////////////////
  /////  Compare simulation results  /////
  ////////////////////////////////////////
  
  if (align.compareSims) {
    align.compareSimulations(radicalNames);
  }

  ///////////////////////////
  /////  Aligning Runs  /////
  ///////////////////////////

  if (align.verbose)  {
    std::cout << "Begin to align runs" << endl;
  }

  std::vector<string> runInds;
  std::map< string, std::map< double, double* > >  diffP_arrays;
  std::map< string, int > runShifts;

  ///// Loop through events in the file and saving to maps  /////
  curScan = -1;
  int pInd;
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);

    // Ignore reference images taken before scan
    if (imgIsRef) {
      if (align.verbose) std::cout << "INFO: Adding reference image.\n";

      if (reference.size() < (*azmAvg).size()) {
        reference.resize((*azmAvg).size(), 0);
        referenceCount.resize((*azmAvg).size(), 0);
      }

      std::transform((*azmAvg).begin(), (*azmAvg).end(),
                      reference.begin(), reference.begin(), 
                      std::plus<double>());

      std::for_each(referenceCount.begin(), (referenceCount.begin() + (*azmAvg).size()),
                      [] (double &d) {
                        d++;
                      });
      continue;
    }
    //if (stagePos < (align.t0StagePos - 0.021)*10000) {
    //  continue;
    //}

    cout<<"ORIG: "<<(*azmAvg).size()<<"  "<<(*legCoeffs).size()<<endl;
    ///  Insert entries to respective maps  ///
    align.addEntry(scan, stagePos, azmAvg, legCoeffs, 1);
  }




  ///////////////////////////////////////////////////////////
  /////  Prune data for outliers and sparse time steps  /////
  ///////////////////////////////////////////////////////////


  //align.removeOutliers();

  cout<<"merging scans"<<endl;
  align.mergeScans();

  cout<<"subtracting t0 and norm"<<endl;
  //align.subtractT0andNormalize();

  cout<<"smearing time"<<endl;
  //align.smearTime();


  for (uint ir=0; ir<align.NradAzmBins; ir++) {
    cout<<"azm: "<<ir<<"  "<<align.azimuthalAvg[12][ir]<<endl;
  }


  cout<<"saving"<<endl;
  cout<<"111"<<endl;
  save::saveDat<double>(align.azimuthalAvg, 
      "./plots/data/data_sMsAzmAvgDiff["
      + to_string(align.azimuthalAvg.size()) + ","
      + to_string(align.azimuthalAvg[0].size()) + "].dat");
 
  for (int ilg=0; ilg<align.Nlegendres; ilg++) {
    ///  Saving data  ///
    cout<<"111"<<endl;
    save::saveDat<double>(align.legendres[ilg], "./plots/data/data_sMsL"
        + to_string(ilg) + "Diff["
        + to_string(align.legendres[ilg].size()) + ","
        + to_string(align.legendres[ilg][0].size()) + "].dat");
    cout<<"222"<<endl;
    save::saveDat<double>(align.legendres[ilg][align.legendres[ilg].size()-1], 
        "./plots/data/data_sMsFinalL"
        + to_string(ilg) + "Diff["
        + to_string(align.legendres[ilg][0].size()) + "].dat");


    cout<<"333"<<endl;
    save::saveDat<double>(align.smearedImg[ilg], "./plots/data/data_sMsL"
        + to_string(ilg) + "DiffSmear["
        + to_string(align.smearSTD) + "["
        + to_string(align.smearedImg[ilg].size()) + ","
        + to_string(align.smearedImg[ilg][0].size()) + "].dat");

    cout<<"444"<<endl;
    save::saveDat<double>(align.smearedImg[ilg][align.stagePosInds.size()-1],
        "./plots/data/data_sMsFinalL"
        + to_string(ilg) + "DiffSmear"
        + to_string(align.smearSTD) + "["
        + to_string(align.smearedImg[ilg][0].size()) + "].dat");
  }



  cout<<"start autocorrelation"<<endl;

  ///////////////////////////////////////
  /////  Pair correlation function  /////
  ///////////////////////////////////////
  if (align.verbose) {
    std::cout << "Begin calculating pair correlation functions\n";
  }

  /*
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
  */


  /////////////////////////
  /////  Cleaning up  /////
  /////////////////////////

  if (align.verbose) {
    std::cout << "Cleaning up" << endl;
  }

  fftw_destroy_plan(fftB);
  fftw_free(rSpace);
  fftw_free(qSpace);
  return 1;
}
