#include "../analysis.h" 
#include "/reg/neh/home/khegazy/baseTools/tools/parameters.h"

using namespace std;


struct bkgStruct {

  int stagePos;
  int time;
  double UVcount;
  double bnkTemp;
  double hbTemp;
  std::vector<double> diffIntensities;
};

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
  }
  else if (iPos == string::npos) {
    iPos = fileList.find("Run-");
  }
  iPos += 4;
  auto fPos = fileList.find(".txt");
  std::string runName = fileList.substr(iPos, 13); //fPos - iPos);
  cout<<"RunName: "<<runName<<endl;
  parameterClass params(runName);

  // Plotting variables
  std::vector<PLOToptions> opts(2);
  std::vector<string> vals(2);
  opts[0] = minimum;  opts[1] = maximum;
  vals[0] = "-5";     vals[1] = "5";

  // Analysis variables
  std::map< int, bkgStruct > valMap;
  std::map< int, std::vector< std::vector<double> > > imgMap;
  std::map< int, double > imgMapCount;
  std::vector< std::vector<double> > plotScan;
  double sum, count;
  int startInd, endInd;

  ///// Loop through events in the file /////
  if (params.verbose)
    cout << "INFO: Starting even loop\n";

  int curScan = scan;
  int64_t scanInd = 0;
  double Nrefs;
  bool makingRef = true;
  std::vector<double> reference(params.NradAzmBins, 0);
  std::vector<double> refSubAzmAvg(params.NradAzmBins, 0);
  std::vector<double> pairCorr(params.maxRbins, 0);
  for (int64_t jevt=0; jevt<Nentries; jevt++) {
    analysis.loadEvent(jevt);

    if (std::find(params.badScans.begin(), params.badScans.end(), scan)
          != params.badScans.end()) {
      continue;
    }

    if (scan != curScan) {
      curScan   = scan;
      scanInd  = jevt;
      makingRef = true;
      Nrefs     = 0;
      std::fill(reference.begin(), reference.end(), 0);
    }
    if (makingRef && imgIsRef) {
      for (uint iq=0; iq<params.NradAzmBins; iq++) {
        reference[iq] += (*azmAvg)[iq];
      }
      Nrefs += 1;
      continue;
    }
    if (makingRef && !imgIsRef) {
      for (uint iq=0; iq<params.NradAzmBins; iq++) {
        reference[iq] /= Nrefs;
      }
      makingRef = false;
      jevt = scanInd - 1;
      continue;
    }

    for (int iq=0; iq<params.NradAzmBins; iq++) {
      refSubAzmAvg[iq] = (*azmAvg)[iq] - reference[iq];
    }

    /////  Calculate Pair Correlation  /////
    /*
    save::saveDat<double>(
        refSubAzmAvg, 
        "./results/currentDiffAzmAvg-TEMPORARY["
          + to_string(params.NradAzmBins) + "].dat");
    system(("./../timeDepStudies/pairCorr.exe " + runName 
        + " -Idir ./results/ -Fname currentDiffAzmAvg-TEMPORARY["
        + to_string(params.NradAzmBins) + "].dat -Osuf -TEMPORARY").c_str());
    save::importDat<double>(
        pairCorr,
        "results/data-" + runName 
        + "-TEMPORARY-pairCorrEven["
        + to_string(params.maxRbins) + "].dat")
    */

    valMap[timeStamp].stagePos  = stagePos;
    valMap[timeStamp].UVcount   = UVcounts;
    valMap[timeStamp].bnkTemp   = bunkerTemp;
    valMap[timeStamp].hbTemp    = highBayTemp;
    valMap[timeStamp].time      = timeStamp;
    for (auto const & rng : params.bkgStudyRanges) {
      sum = 0;
      count = 0;
      startInd  = (int)(params.NradAzmBins*rng[0]/params.maxQazm);
      endInd    = (int)(params.NradAzmBins*rng[1]/params.maxQazm);
      for (int ir=startInd; ir<=endInd; ir++) {
        sum += refSubAzmAvg[ir];
        count += 1;
      }
      valMap[timeStamp].diffIntensities.push_back(sum/count);
    }

    if (!imgMap[stagePos].size()) {
      imgMap[stagePos].resize((*imgOrig).size());
      for (uint ir=0; ir<(*imgOrig).size(); ir++) {
        imgMap[stagePos][ir].resize((*imgOrig)[ir].size());
      }
      imgMapCount[stagePos] = 0;
    }

    for (uint ir=0; ir<(*imgOrig).size(); ir++) {
      for (uint ic=0; ic<(*imgOrig)[ir].size(); ic++) {
        imgMap[stagePos][ir][ic] += (*imgOrig)[ir][ic];
        imgMapCount[stagePos] += 1;
      }
    }

    //system("rm results/*TEMPORARY*");
  }

  cout<<"begin saving"<<endl;
  /////  Saving  /////
  std::map<int, std::vector<bkgStruct*> > stgPosMap;

  int ind = 0;
  std::vector<double> UVcountVec(valMap.size());
  std::vector<double> bnkTempVec(valMap.size());
  std::vector<double> hbTempVec(valMap.size());
  std::vector<double> labTimeVec(valMap.size());
  std::vector< std::vector<double> > diffIntVec(params.bkgStudyRanges.size());
  for (uint j=0; j<params.bkgStudyRanges.size(); j++) {
    diffIntVec[j].resize(valMap.size());
  }
  for (auto &bst : valMap) {
    stgPosMap[bst.second.stagePos].push_back(&bst.second);
    UVcountVec[ind] = bst.second.UVcount;
    bnkTempVec[ind] = bst.second.bnkTemp;
    hbTempVec[ind]  = bst.second.hbTemp;
    labTimeVec[ind] = bst.first;
    for (uint j=0; j<params.bkgStudyRanges.size(); j++) {
      diffIntVec[j][ind] = bst.second.diffIntensities[j];
    }
    ind++;
  }
  cout<<"saving time dep"<<endl;
  save::saveDat<double>(labTimeVec, 
      "./results/backgroundStudy/data-" 
        + runName + "-labTime["
        + to_string(UVcountVec.size()) + "].dat");
  save::saveDat<double>(UVcountVec, 
      "./results/backgroundStudy/data-" 
        + runName + "-UVcountOverTime["
        + to_string(UVcountVec.size()) + "].dat");
  save::saveDat<double>(bnkTempVec, 
      "./results/backgroundStudy/data-" 
        + runName + "-bunkerTempOverTime["
        + to_string(bnkTempVec.size()) + "].dat");
  save::saveDat<double>(hbTempVec, 
      "./results/backgroundStudy/data-" 
        + runName + "-highBayTempOverTime["
        + to_string(hbTempVec.size()) + "].dat");
  for (uint j=0; j<params.bkgStudyRanges.size(); j++) {
    save::saveDat<double>(diffIntVec[j], 
        "./results/backgroundStudy/data-" 
          + runName + "-diffractionAvg-"
          + to_string(params.bkgStudyRanges[j][0]) + "-"
          + to_string(params.bkgStudyRanges[j][1]) + "-bins["
          + to_string(diffIntVec[j].size()) + "].dat");
  }

  cout<<"looping over stages"<<endl;
  // Results at each stage position
  for (auto & img : imgMap) {
    for (uint ir=0; ir<img.second.size(); ir++) {
      for (uint ic=0; ic<img.second[ir].size(); ic++) {
        img.second[ir][ic] /= imgMapCount[img.first];
      }
    }

    cout<<"saving image"<<endl;
    save::saveDat<double>(
        img.second, 
        "./results/backgroundStudy/data-" + runName
          + "-stagePos-" + to_string(img.first) 
          + "-averageImage[" + to_string(img.second.size())
          + "," + to_string(img.second[0].size()) + "].dat");
    
    cout<<"filling arrays"<<endl;
    UVcountVec.resize(stgPosMap[img.first].size());
    bnkTempVec.resize(stgPosMap[img.first].size());
    hbTempVec.resize(stgPosMap[img.first].size());
    labTimeVec.resize(stgPosMap[img.first].size());
    for (uint j=0; j<params.bkgStudyRanges.size(); j++) {
      diffIntVec[j].resize(stgPosMap[img.first].size());
    }
    for (uint i=0; i<stgPosMap[img.first].size(); i++) {
      UVcountVec[i] = stgPosMap[img.first][i]->UVcount;
      bnkTempVec[i] = stgPosMap[img.first][i]->bnkTemp;
      hbTempVec[i]  = stgPosMap[img.first][i]->hbTemp;
      labTimeVec[i] = stgPosMap[img.first][i]->time;
      for (uint j=0; j<params.bkgStudyRanges.size(); j++) {
        diffIntVec[j][i] = stgPosMap[img.first][i]->diffIntensities[j];
      }
    }
    cout<<"saving"<<endl;
    save::saveDat<double>(labTimeVec, 
        "./results/backgroundStudy/data-" + runName
          + "-stagePos-" + to_string(img.first)  
          + "-labTime["
          + to_string(UVcountVec.size()) + "].dat");
    save::saveDat<double>(UVcountVec, 
        "./results/backgroundStudy/data-" + runName 
          + "-stagePos-" + to_string(img.first)
          + "-UVcountOverTime[" 
          + to_string(UVcountVec.size()) + "].dat");
    save::saveDat<double>(bnkTempVec, 
        "./results/backgroundStudy/data-" + runName 
          + "-stagePos-" + to_string(img.first)
          + "-bunkerTempOverTime["
          + to_string(bnkTempVec.size()) + "].dat");
    save::saveDat<double>(hbTempVec, 
        "./results/backgroundStudy/data-" + runName 
          + "-stagePos-" + to_string(img.first)
          + "-highBayTempOverTime["
          + to_string(hbTempVec.size()) + "].dat");
    for (uint j=0; j<params.bkgStudyRanges.size(); j++) {
      save::saveDat<double>(diffIntVec[j], 
          "./results/backgroundStudy/data-" 
            + runName + + "-stagePos-"
            + to_string(img.first) + "-diffractionAvg-"
            + to_string(params.bkgStudyRanges[j][0]) + "-"
            + to_string(params.bkgStudyRanges[j][1]) + "-bins[" 
            + to_string(diffIntVec[j].size()) + "].dat");
    }
  }


  cout<<"done"<<endl;

  return 0;
}
