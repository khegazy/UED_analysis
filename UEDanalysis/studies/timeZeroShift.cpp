#include "../analysis.h" 


using namespace std;


int main(int argc, char* argv[]) {

  if (argc<2) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe 'file.txt' !!!"<<endl;
    cerr<<"         Can also run ./analysis.exe 'file.txt' 'treeName'"<<endl;
    exit(0);
  }

  std::string fileSuffix = "";
  string file(argv[1]);
  for (int iarg=2; iarg<argc; iarg+=2) {
    if (strcmp(argv[iarg], "-fileNameSuffix") == 0) {
      string str(argv[iarg+1]);
      fileSuffix = str;
    }
  }

  auto iPos = file.find("data-");
  iPos += 5;
  auto fPos = file.find("-", iPos);
  std::string runName = file.substr(iPos, 13); //fPos - iPos);
  cout<<"RunName: "<<runName<<endl;
  parameterClass params(runName);

  // Plot options
  vector<PLOToptions> opts(3);
  vector<string>  vals(3);
  opts[0] = maximum;  vals[0] = "10";
  opts[1] = minimum;  vals[1] = "-10";
  opts[2] = xSpan;    vals[2] = "0," + to_string(params.maxQazm);


  // Get Data
  auto sPos = file.rfind("/");
  std::string folder = file.substr(0, sPos);
  std::vector<int> shape  = save::getShape(
                              file.substr(0, sPos), 
                              file.substr(sPos+1, file.length()-sPos-4));
  std::vector< std::vector<double> > tdDiffraction(shape[0]);
  for (int i=0; i<shape[0]; i++) {
    tdDiffraction[i].resize(shape[1]);
  }
  save::importDat<double>(tdDiffraction, file);

  // Get ranges
  double count;
  int startInd, endInd;
  std::vector< std::vector<double> > rangeMeans(params.tZeroQranges.size());
  for (int ir=0; ir<(int)params.tZeroQranges.size(); ir++) {
    startInd  = (int)(params.NradAzmBins
                          *params.tZeroQranges[ir][0]/params.maxQazm);
    endInd    = (int)(params.NradAzmBins
                          *params.tZeroQranges[ir][1]/params.maxQazm);

    rangeMeans[ir].resize(shape[0], 0);
    for (int it=0; it<shape[0]; it++) {
      count = 0;
      rangeMeans[ir][it] = 0;
      for (int iq=startInd; iq<endInd; iq++) {
        rangeMeans[ir][it] += tdDiffraction[it][iq];
        count += 1;
      }
      rangeMeans[ir][it] /= count;
    }
  }

  for (uint ir=0; ir<rangeMeans.size(); ir++) {
    save::saveDat<double>(rangeMeans[ir],
        "./results/data-" + runName 
        + "-Qmean" + to_string(params.tZeroQranges[ir][0])
        + "-" + to_string(params.tZeroQranges[ir][1]) 
        + fileSuffix
        + "-Bins[" + to_string(rangeMeans[ir].size()) + "].dat");
  }

  return 1;
}
