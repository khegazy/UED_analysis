#include "../analysis.h" 
#include <fstream>

using namespace std;


int main(int argc, char* argv[]) {

  if (argc<3) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe runName '/reg/..../scanSearch/sizeXX' !!!\n";
    exit(0);
  }

  // Get command line variables
  std::string fileSuffix = "";
  string folder(argv[2]);
  for (int iarg=2; iarg<argc; iarg+=2) {
    if (strcmp(argv[iarg], "-fileNameSuffix") == 0) {
      string str(argv[iarg+1]);
      fileSuffix = str;
    }
  }

  // Get runName
  std::string runName(argv[1]);

  // Get size of merged scans
  int size;
  auto iPos = folder.find("size") + 4;
  auto fPos = folder.find("/", iPos);
  if (fPos = std::string::npos) {
    fPos = folder.length();
  }
  size = std::atoi(folder.substr(iPos, fPos-iPos).c_str());

  // Get parameters
  parameterClass params(runName);

  // Plot options
  vector<PLOToptions> opts(3);
  vector<string>  vals(3);
  opts[0] = maximum;  vals[0] = "10";
  opts[1] = minimum;  vals[1] = "-10";
  opts[2] = xSpan;    vals[2] = "0," + to_string(params.maxQazm);


  /////////////////////////////////////////////
  /////  Loop through files and evaluate  /////
  /////////////////////////////////////////////
  int index = 1;
  bool foundFile;
  std::string fileID, fileName;
  std::string timeDepFile, referenceFile;
  while (true) {

    // File names
    fileID = "scanLines" + to_string(index)
        + "-" + to_string(index+size-1);
    timeDepFile = 
        "data-" + runName + "-"
        + fileID + "-azmAvgDiff[";
    referenceFile = 
        "data-" + runName + "-"
        + fileID + "-referenceAzm[" 
        + to_string(params.NradAzmBins) + "].dat";

    // Check that files exist
    foundFile = false;
    DIR* dir = opendir(folder.c_str());
    struct dirent* ent;
    while ((ent = readdir(dir)) != NULL) {
      string curFileName(ent->d_name);
      iPos = curFileName.find("scanLines");
      if (iPos == std::string::npos) continue;
      fPos = curFileName.find("-", iPos);
      fPos = curFileName.find("-", fPos+1);
      if (curFileName.substr(iPos, fPos-iPos).compare(fileID) == 0) {
        foundFile = true;
        break;
      }
    }
    if (!foundFile) {
      cout << "Ending loop at index: " << index <<endl;
      break;
    }
    fileID = "-" + fileID;

    // Get Data
    std::vector<int> shape  = save::getShape(folder, timeDepFile);
    timeDepFile +=
        to_string(shape[0]) + "," 
        + to_string(params.NradAzmBins) + "].dat";
    std::vector< std::vector<double> > tdDiffraction(shape[0]);
    for (int i=0; i<shape[0]; i++) {
      tdDiffraction[i].resize(shape[1]);
    }
    save::importDat<double>(tdDiffraction, folder + "/" + timeDepFile);

    
    std::vector<double> reference(params.NradAzmBins, 0);
    save::importDat<double>(reference, folder + "/" + referenceFile); 


    // Get ranges
    double count;
    int startInd, endInd;
    std::vector< std::vector<double> > rangeMeans(params.tZeroQranges.size());
    std::vector<double> rangeRefMeans(params.tZeroQranges.size());
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
      
      count = 0;
      for (int iq=startInd; iq<endInd; iq++) {
        rangeRefMeans[ir] += reference[iq];
        count += 1;
      }
      rangeRefMeans[ir] /= count;

    }

   
    // Tzero ratio
    std::vector<double> ratio(shape[0]);
    for (int it=0; it<shape[0]; it++) {
      ratio[it] = (rangeMeans[params.tZeroRatio[0]][it] 
                  + rangeRefMeans[params.tZeroRatio[0]])/
                  (rangeMeans[params.tZeroRatio[1]][it] 
                  + rangeRefMeans[params.tZeroRatio[1]]);
    }
    save::saveDat<double>(ratio, 
          "./results/data-" + runName 
          + fileSuffix + fileID + "-tZeroRatio-Bins[" 
          + to_string(shape[0]) + "].dat");



    /////  Saving  /////
    for (uint ir=0; ir<rangeMeans.size(); ir++) {
      save::saveDat<double>(rangeMeans[ir],
          "./results/data-" + runName 
          + "-Qmean" + to_string(params.tZeroQranges[ir][0])
          + "-" + to_string(params.tZeroQranges[ir][1]) 
          + fileSuffix + fileID +
          + "-Bins[" + to_string(rangeMeans[ir].size()) + "].dat");
    }

    save::saveDat<double>(rangeRefMeans,
        "./results/data-" + runName 
        + "-referenceQmeans" + fileID + fileSuffix
        + "-Bins[" + to_string(rangeRefMeans.size()) + "].dat");

    index++;
  }
  return 1;
}
