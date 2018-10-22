#include "saving.h"
#include "dirent.h"

using namespace std;

void save::saveDat(std::vector<double> &input, std::string fileName) {
  FILE* output = fopen(fileName.c_str(), "wb");
  fwrite(&input[0], sizeof(double), input.size(), output);
  fclose(output);

  return;
}


std::vector<int> save::getShape(std::string folder, std::string filePrefix) {

  // Get file name
  std::string fileName;
  DIR* dir = opendir(folder.c_str());
  struct dirent* ent;
  while ((ent = readdir(dir)) != NULL) {
    string curFileName(ent->d_name);
    if (filePrefix.compare(curFileName.substr(0,filePrefix.length())) == 0) {
      fileName = curFileName;
      break;
    }
  }

  cout << "INFO: Getting shape from " + fileName <<endl;

  std::vector<int> inds;
  int iPos = fileName.find("[") + 1;
  int fPos = fileName.find(",", iPos);
  cout<<"fileName: "<<fileName<<endl;
  while (fPos != string::npos) {
    cout<<"  "<<fileName.substr(iPos, fPos-iPos)<<endl;
    inds.push_back(stoi(fileName.substr(iPos, fPos-iPos)));
    iPos = fPos + 1;
    fPos = fileName.find(",", iPos);
  }
  fPos = fileName.find("]", iPos);
  cout<<"  "<<fileName.substr(iPos, fPos-iPos)<<endl;
  inds.push_back(stoi(fileName.substr(iPos, fPos-iPos)));

  return inds;
}






/*
template <typename type>
void save::saveDat(std::vector< std::vector<type> > &input, std::string fileName) {
  FILE* output = fopen(fileName.c_str(), "wb");
  for (uint i=0; i<input.size(); i++) {
    fwrite(&input[i][0], sizeof(type), input[i].size(), output);
  }
  fclose(output);

  return;
}

template <typename type>
void save::saveDat(std::vector< std::vector< std::vector<type> > > &input, 
    std::string fileName) {
  FILE* output = fopen(fileName.c_str(), "wb");
  for (uint i=0; i<input.size(); i++) {
    for (uint ii=0; ii<input[i].size(); ii++) {
      fwrite(&input[i][ii][0], sizeof(type), input[i][ii].size(), output);
    }
  }
  fclose(output);

  return;
}
*/
