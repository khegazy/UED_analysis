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

  vector<double> refAvg(params.NradAzmBins,0);
  double count = 0;
  // Plot options
  vector<PLOToptions> opts(2);
  vector<string>  vals(2);
  opts[0] = maximum;  vals[0] = "200";
  opts[1] = minimum;  vals[1] = "-200";

  // Get reference
  vector<double> referenceAzm(params.NradAzmBins, 0);
  vector<string> runNames = {"20180627_1551", "20180629_1630", 
                             "20180630_1925", "20180701_0746"};
  for (uint i=0; i<runNames.size(); i++) {
    if (runName.compare(runNames[i]) == 0) continue;

    save::importDat(referenceAzm,
        "../staticDiffraction/results/references-"
        + runNames[i] + ".dat");

    std::transform(referenceAzm.begin(), referenceAzm.end(),
        refAvg.begin(), refAvg.begin(),
        std::plus<double>());
  }

  std::for_each(refAvg.begin(), refAvg.end(),
      [&runNames] (double &d) {
        d /= -1*((double)runNames.size() - 1);
      });

  save::importDat(referenceAzm,
      "../staticDiffraction/results/references-"
      + runName + ".dat");

  std::transform(referenceAzm.begin(), referenceAzm.end(),
      refAvg.begin(), refAvg.begin(),
      std::plus<double>());

  plt.print1d(refAvg, "./plots/referenceCompAllRef_" + runName,
                opts, vals);

  cout<<"filling"<<endl;
  std::fill(refAvg.begin(), refAvg.end(), 0);
  cout<<"filled"<<endl;
 

  // Get sms norms
  vector<double> atmAzmDiff(params.NradAzmBins, 0.0);
  vector<double> sMsAzmNorm(params.NradAzmBins, 0.0);
  string simFileNameSuffix = "_Bins-" + to_string(params.NradAzmBins)
                      + "_Qmax-" + to_string(params.maxQazm)
                      + "_Ieb-" + to_string(params.Iebeam)
                      + "_scrnD-" + to_string(params.screenDist)
                      + "_elE-" + to_string(params.elEnergy) + ".dat";

  save::importDat<double>(atmAzmDiff, params.simReferenceDir + "/"
              + params.molName + "_atmDiffractionPatternLineOut"
              + simFileNameSuffix);

  for (int ir=0; ir<params.NradAzmBins; ir++) {
    sMsAzmNorm[ir] = (params.maxQazm*(ir+0.5)/params.NradAzmBins)
                        /(atmAzmDiff[ir]);
  }


  ///// Loop through events in the file /////
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);

    // Skip bad scans
    if (std::find(params.badScans.begin(), params.badScans.end(), scan)
          != params.badScans.end()) continue;

    // Skip reference images
    if (imgIsRef) continue;

    // Skip post T0 images
    if (stagePos > 1543050) continue;

    std::transform((*azmAvg).begin(), (*azmAvg).end(),
        refAvg.begin(), refAvg.begin(),
        std::plus<double>());

    count++;
  }

  cout<<"normalizing"<<endl;
  std::for_each(refAvg.begin(), refAvg.end(),
           [count] (double &d) {
             d /= -1*count;
           });

  cout<<"adding"<<endl;
  std::transform(referenceAzm.begin(), referenceAzm.end(),
      refAvg.begin(), refAvg.begin(),
      std::plus<double>());

  plt.print1d(refAvg, "./plots/referenceTest_" + runName,
                maximum, "4");


  return 1;
}
