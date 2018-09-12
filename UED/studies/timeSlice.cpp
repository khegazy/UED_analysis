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

  map<int, map<int, vector<double> > > azmAvgs;

  auto iPos = fileList.find("run-");
  if (iPos == string::npos) {
    iPos = fileList.find("RUN-");
  }
  iPos += 4;
  auto fPos = fileList.find(".txt");
  std::string runName = fileList.substr(iPos, 13); //fPos - iPos);
  cout<<"RunName: "<<runName<<endl;
  parameterClass params(runName);

  // Plot options
  vector<PLOToptions> opts(2);
  vector<string>  vals(2);
  opts[0] = maximum;  vals[0] = "75";
  opts[1] = minimum;  vals[1] = "-75";

  // Get sms norms
  vector<double> atmAzmDiff(params.NradAzmBins, 0.0);
  vector<double> sMsAzmNorm(params.NradAzmBins, 0.0);
  string simFileNameSuffix = 
      "_Qmax-" + to_string(params.maxQazm)
      + "_Ieb-" + to_string(params.Iebeam)
      + "_scrnD-" + to_string(params.screenDist)
      + "_elE-" + to_string(params.elEnergy) 
      + "_Bins[" + to_string(params.NradAzmBins) + "].dat";

  save::importDat<double>(atmAzmDiff, params.simReferenceDir + "/"
              + params.molName + "_atmDiffractionPatternLineOut"
              + simFileNameSuffix);

  for (int ir=0; ir<params.NradAzmBins; ir++) {
    sMsAzmNorm[ir] = (params.maxQazm*(ir+0.5)/params.NradAzmBins)
                        /(atmAzmDiff[ir]);
  }


  ///// Loop through events in the file /////
  int curScan = scan;
  double refCount = 0;
  vector<double> localRef(params.NradAzmBins, 0);
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);

    // Skip bad scans
    if (std::find(params.badScans.begin(), params.badScans.end(), scan)
          != params.badScans.end()) continue;


    if (scan != curScan) {
      std::for_each(localRef.begin(), localRef.end(),
          [&refCount] (double &d) {
            d /= refCount;
          });

      for (auto& stgItr : azmAvgs) {
        for (int i=0; i<params.NradAzmBins; i++) {
          stgItr.second[curScan][i] -= localRef[i];
        }
      }

      std::fill(localRef.begin(), localRef.end(), 0);
      refCount = 0;
      curScan = scan;
    }

    if (imgIsRef) {
      std::transform((*azmAvg).begin(), (*azmAvg).end(),
          localRef.begin(), localRef.begin(),
          std::plus<double>());
      refCount++;
    }

    // Only looking for these time slices
    //if ((stagePos != 1548450) && (stagePos != 1546950)) 
    if ((stagePos != 1548450) && (stagePos != 1549200)) 
      continue;

    azmAvgs[stagePos][scan].resize(params.NradAzmBins);
    for (int i=0; i<params.NradAzmBins; i++) {
      azmAvgs[stagePos][scan][i] = (*azmAvg)[i];
    }
  }

  std::for_each(localRef.begin(), localRef.end(),
      [&refCount] (double &d) {
        d /= refCount;
      });

  for (auto& stgItr : azmAvgs) {
    for (int i=0; i<params.NradAzmBins; i++) {
      stgItr.second[curScan][i] -= localRef[i];
    }
  }

  /////  Subtracting average  /////
  for (auto& stgItr : azmAvgs) {
  /*
    vector<double> mean(params.NradAzmBins, 0);
    for (auto& scnItr : stgItr.second) {
      for (int i=0; i<params.NradAzmBins; i++) {
        mean[i] += scnItr.second[i];
      }
    }

    for (int i=0; i<params.NradAzmBins; i++) {
      mean[i] /= (double)stgItr.second.size();
    }

    for (auto& scnItr : stgItr.second) {
      for (int i=0; i<params.NradAzmBins; i++) {
        scnItr.second[i] -= mean[i];
        //scnItr.second[i] *= sMsAzmNorm[i];
      }
    }
    */

    for (auto& scnItr : stgItr.second) {
      plt.print1d(scnItr.second, "./plots/timeSlice/timeSlice_" 
          + to_string(stgItr.first) 
          + "_" + to_string(scnItr.first),
          opts, vals);

      save::saveDat<double>(scnItr.second, "./plots/data/timeSlice-"
          + to_string(stgItr.first)
          + "_" + to_string(scnItr.first) + ".dat");
    }
  }

  return 1;
}
