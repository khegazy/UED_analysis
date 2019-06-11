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

  auto iPos = fileList.find("run-");
  if (iPos == string::npos) {
    iPos = fileList.find("RUN-");
  }
  iPos += 4;
  auto fPos = fileList.find(".txt");
  std::string runName = fileList.substr(iPos, 13); //fPos - iPos);
  cout<<"RunName: "<<runName<<endl;
  parameterClass params(runName);

  std::vector< std::pair<int, int> > radRanges;
  std::pair<int, int> p{168, 213};
  radRanges.push_back(p);

  std::vector< std::vector<std::string> > pltVals(radRanges.size());
  std::vector< PLOToptions> pltOpts(2);
  pltOpts[0] = maximum; pltOpts[1] = minimum;
  std::vector< std::vector< std::vector<double> > > images(radRanges.size());
  for (int i=0; i<(int)images.size(); i++) {
    images[i].resize((*imgSubBkg).size());
    for (uint ir=0; ir<(*imgSubBkg).size(); ir++) {
      images[i][ir].resize((*imgSubBkg)[ir].size(), 0);
    }
    pltVals[i].resize(pltOpts.size());
    pltVals[i][0] = "90";
    pltVals[i][1] = "50";
  }


  ///// Loop through events in the file /////
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);

    plt.printRC((*imgSubBkg), 
        "./plots/diffPatterns/data-" + runName
        + "_scan-" + to_string(scan)
        + "_stagePos-" + to_string(stagePos));
    save::saveDat<double>((*imgSubBkg),
        "./results/diffPatterns/data-" + runName
        + "_scan-" + to_string(scan)
        + "_stagePos-" + to_string(stagePos)
        + "_bins[" + to_string((*imgSubBkg).size())
        + "," + to_string((*imgSubBkg)[0].size()) + "].dat");

    for (int ir=0; ir<(int)(*imgSubBkg).size(); ir++) {
      for (int ic=0; ic<(int)(*imgSubBkg)[ir].size(); ic++) {
        int rad = std::sqrt(std::pow(ir-centerR, 2) + std::pow(ic-centerC, 2));
        for (uint i=0; i<radRanges.size(); i++) {
          if ((rad > radRanges[i].first) && (rad < radRanges[i].second)) {
            images[i][ir][ic] = (*imgSubBkg)[ir][ic];
          }
        }
      }
    }
        
    for (uint i=0; i<radRanges.size(); i++) {
      plt.printRC(images[i], 
          "./plots/diffPatterns/data-" + runName
          + "_scan-" + to_string(scan)
          + "_stagePos-" + to_string(stagePos)
          + "_range-" + to_string(i),
          pltOpts, pltVals[i]);
    }
  }

  return 1;
}
