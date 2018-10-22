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
    iPos = fileList.find("RUN-");
  }
  iPos += 4;
  auto fPos = fileList.find(".txt");
  std::string runName = fileList.substr(iPos, 13); //fPos - iPos);
  cout<<"RunName: "<<runName<<endl;
  parameterClass params(runName);

  std::vector<double> reference(params.NradAzmBins);
  save::importDat<double>(reference,
       "../mergeScans/results/referenceAzm-" + runName +
       "[" + to_string(params.NradAzmBins) + "].dat");
  plt.print1d(reference, "./plots/referenceImg_ " + runName + ".dat");

  //std::vector<int> plotDelays = {1543050, 1543150};
  //std::vector<int> plotDelays = {1544550, 1544400, 1544250, 1544100};
  //std::vector<int> plotDelays = {1542350,1542450,1542550,1542650};
  std::vector<int> plotDelays = {1530000, 1530100, 1542350,1542450,1542550,1542650};
  //std::vector<int> plotDelays = {1550700, 1549200};
  std::vector<PLOToptions> opts = {maximum, minimum};
  std::vector<std::string> vals = {"1200", "-200"};

  int count = 0;
  std::vector<double> diff(params.NradAzmBins);
  ///// Loop through events in the file /////
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);

    if (std::find(params.badScans.begin(), params.badScans.end(), scan)
          != params.badScans.end()) continue;

    if (std::find(plotDelays.begin(), plotDelays.end(), stagePos) 
          != plotDelays.end()) {
      plt.printRC((*imgSubBkg), "./plots/image_" 
          + runName + "_scan-" 
          + to_string(scan) + "_pos-" 
          + to_string(stagePos), opts, vals);
      plt.print1d((*filtAzmAvg), "./plots/lineOut_"
          + runName + "_scan-" 
          + to_string(scan) + "_pos-" 
          + to_string(stagePos));
      for (int iq=0; iq<params.NradAzmBins; iq++) {
        if (reference[iq] != 0 && (*filtAzmAvg)[iq] != NANVAL) {
          diff[iq] = (*filtAzmAvg)[iq] - reference[iq];
          if (count) {
            cout<<"COUNT: "<<count<<endl;
            count =0;
          }
        }
        else {
          count++;
          diff[iq] = 0;
        }
      }
      vals[0] = "0.2";
      vals[1] = "-0.2";
      plt.print1d(diff, "./plots/diffLineOut_"
          + runName + "_scan-" 
          + to_string(scan) + "_pos-" 
          + to_string(stagePos), opts, vals);
    }
  }


  return 1;
}
