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
  vector<int> vCenterC;
  vector<int> vCenterR;
  vector<int> imgID;

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

    imgID.push_back(scan + 0.01*imgNum);
    vCenterR.push_back(centerR);
    vCenterC.push_back(centerC);
  }

  vector<int> orderedInds(imgID.size());
  iota(orderedInds.begin(), orderedInds.end(), 0);
  sort(orderedInds.begin(), orderedInds.end(),
      [&imgID](int i1, int i2)
      {return imgID[i1] < imgID[i2];});

  vector<double> pltCenterR(imgID.size());
  vector<double> pltCenterC(imgID.size());

  for (int i=0; i<(int)imgID.size(); i++) {
    pltCenterR[i] = vCenterR[orderedInds[i]];
    pltCenterC[i] = vCenterC[orderedInds[i]];
  }

  plt.print1d(pltCenterR, "./plots/centerR_" + runName);
  plt.print1d(pltCenterC, "./plots/centerC_" + runName);

  return 1;
}
