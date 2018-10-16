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

  // Get run name
  auto iPos = fileList.find("run-");
  if (iPos == string::npos) {
    iPos = fileList.find("Run-");
    cout<<"here1"<<endl;
  }
  else if (iPos == string::npos) {
    iPos = fileList.find("Run-");
    cout<<"here2"<<endl;
  }
  cout<<"pos: "<<iPos<<endl;
  iPos += 4;
  auto fPos = fileList.find(".txt");
  std::string runName = fileList.substr(iPos, 13); //fPos - iPos);
  cout<<"RunName: "<<runName<<endl;
  parameterClass params(runName);

  std::vector<double> reference(params.NradAzmBins);
  save::importDat<double>(reference,
       "../mergeScans/results/referenceAzm-" + runName +
       "[" + to_string(params.NradAzmBins) + "].dat");
  std::vector<PLOToptions> opts(2);
  std::vector<string> vals(2);
  opts[0] = minimum;  opts[1] = maximum;
  vals[0] = "-5";     vals[1] = "5";

  float Nref = 0;
  std::vector< std::vector<double> > plotScan;
  std::vector<double> refAzmAvg(params.NradAzmBins,0);
  ///// Loop through events in the file /////
  for (uint64_t jevt=0; jevt<Nentries; jevt++) {
    analysis.loadEvent(jevt);

    if (imgIsRef) {

      if (std::find(params.badScans.begin(), params.badScans.end(), scan)
            != params.badScans.end()) {
        continue;
      }
      auto badImgItr = params.badImages.find(scan);
      if (badImgItr != params.badImages.end()) {
        if (std::find(badImgItr->second.begin(), badImgItr->second.end(), stagePos)
              != badImgItr->second.end()) {
          continue;
        }
      }

      for (int iq=0; iq<params.NradAzmBins; iq++) {
        refAzmAvg[iq] = (refAzmAvg[iq]*Nref + (*filtAzmAvg)[iq])/(Nref + 1);
      }
      Nref++;
    }
    else {
      break;
    }
  }


  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);

    std::vector<double> curAzmAvg(params.NradAzmBins);
    plotScan.push_back(curAzmAvg);

    if (std::find(params.badScans.begin(), params.badScans.end(), scan)
          != params.badScans.end()) {
      for (uint i=0; i<(*filtAzmAvg).size(); i++) {
        plotScan[ievt][i] = 0;
      }
      continue;
    }
    auto badImgItr = params.badImages.find(scan);
    if (badImgItr != params.badImages.end()) {
      if (std::find(badImgItr->second.begin(), badImgItr->second.end(), stagePos)
            != badImgItr->second.end()) {
        for (uint i=0; i<(*filtAzmAvg).size(); i++) {
          plotScan[ievt][i] = 0;
        }
        continue;
      }
    }

    for (uint i=0; i<(*filtAzmAvg).size(); i++) {
      if ((*filtAzmAvg)[i] != NANVAL) {
        plotScan[ievt][i] = ((*filtAzmAvg)[i] - refAzmAvg[i])*i;
      }
      else {
        plotScan[ievt][i] = 0;
      }
    }
  }

  plt.printXY(plotScan,
      "/reg/ued/ana/scratch/nitroBenzene/badScanSearch/azmAvgs-"
      + runName + "_scan-" + to_string(scan),
      opts,
      vals);

  return 1;
}
