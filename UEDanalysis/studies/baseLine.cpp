#include "../analysis.h" 


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

  int start = 100;
  int len = 11;
  vector<int> orderedInds(len);
  vector<double> readout, amplitude, norm;
  vector<double> pressures, pressureDers;
  vector<double> med(len);
  vector<int> scanNum;

  auto iPos = fileList.find("run-");
  if (iPos == string::npos) {
    iPos = fileList.find("RUN-");
  }
  iPos += 4;
  auto fPos = fileList.find(".txt");
  std::string runName = fileList.substr(iPos, 13); //fPos - iPos);
  cout<<"RunName: "<<runName<<endl;
  parameterClass params(runName);

  ///// Loop through events in the file /////
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);

    if (std::find(params.badScans.begin(), params.badScans.end(), scan)
          != params.badScans.end()) continue;

    scanNum.push_back(scan + 0.01*imgNum);
    readout.push_back(readoutNoise);
    norm.push_back(imgNorm);
    pressures.push_back(pressure);
    pressureDers.push_back(pressureDer);

    
    for (int i=start; i<start+len; i++) {
      med[i-start] = (*azmAvg)[i];
    }

    iota(orderedInds.begin(), orderedInds.end(), 0);
    sort(orderedInds.begin(), orderedInds.end(),
        [&med](int i1, int i2)
        {return med[i1] < med[i2];});

    amplitude.push_back(med[orderedInds[len/2]]);

  }

  orderedInds.resize(scanNum.size());

  iota(orderedInds.begin(), orderedInds.end(), 0);
  sort(orderedInds.begin(), orderedInds.end(),
      [&scanNum](int i1, int i2)
      {return scanNum[i1] < scanNum[i2];});

  vector<double> pltReadout(scanNum.size());
  vector<double> pltAmplitude(scanNum.size());
  vector<double> pltNorm(scanNum.size());
  vector<double> pltPressures(scanNum.size());
  vector<double> pltPressureDers(scanNum.size());

  for (int i=0; i<(int)scanNum.size(); i++) {
    pltReadout[i]      = readout[orderedInds[i]];
    pltAmplitude[i]    = amplitude[orderedInds[i]];
    pltNorm[i]         = norm[orderedInds[i]];
    pltPressures[i]    = pressures[orderedInds[i]];
    pltPressureDers[i] = pressureDers[orderedInds[i]];
    if (i>750 && i<1000) {
      cout<<i<<"  "<<scanNum[orderedInds[i]]<<"  "<<pressures[orderedInds[i]]<<"  "<<norm[orderedInds[i]]<<endl;
    }
    if (i>1400 && i<1600) {
      cout<<i<<"  "<<scanNum[orderedInds[i]]<<"  "<<pressures[orderedInds[i]]<<"  "<<norm[orderedInds[i]]<<endl;
    }
    if (i>1750 && i<2200) {
      cout<<i<<"  "<<scanNum[orderedInds[i]]<<"  "<<pressures[orderedInds[i]]<<"  "<<norm[orderedInds[i]]<<endl;
    }
    if (i>2600 && i<2800) {
      cout<<i<<"  "<<scanNum[orderedInds[i]]<<"  "<<pressures[orderedInds[i]]<<"  "<<norm[orderedInds[i]]<<endl;
    }
  }

  plt.print1d(pltReadout, "./plots/readoutNoise_" + runName);
  plt.print1d(pltAmplitude, "./plots/diffAmplitude_" + runName);
  plt.print1d(pltNorm, "./plots/imgNorm_" + runName);
  plt.print1d(pltPressures, "./plots/pressures_" + runName);
  plt.print1d(pltPressureDers, "./plots/pressureDers_" + runName);

  save::saveDat(pltReadout, "./plots/data/readoutNoise_" + runName + ".dat");
  save::saveDat(pltAmplitude, "./plots/data/diffAmplitude_" + runName + ".dat");
  save::saveDat(pltNorm, "./plots/data/imgNorm_" + runName + ".dat");
  save::saveDat(pltPressures, "./plots/data/pressures_" + runName + ".dat");
  save::saveDat(pltPressureDers, "./plots/data/pressureDers_" + runName + ".dat");

  return 1;
}
