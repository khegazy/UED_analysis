#ifndef SAVING_H
#define SAVING_H

//#include <chrono>
//#include <thread>
//#include <H5Cpp.h>
#include "tools.h"


namespace save {

  ////////////////////////////
  /////  Saving vectors  /////
  ////////////////////////////

  void saveDat(std::vector<double> &input, std::string fileName);

  template <typename type>
  void saveDat(type* input, int size, std::string fileName);

  template <typename type>
  void saveDat(std::vector<type> &input, std::string fileName);

  template <typename type>
  void saveDat(std::vector< std::vector<type> > &input, std::string fileName);

  template <typename type>
  void saveDat(std::vector< std::vector< std::vector<type> > > &input, 
      std::string fileName);


  /////////////////////////////
  /////  Importing files  /////
  /////////////////////////////

  template <typename type>
  void importDat(
      std::vector<type> &import,
      std::string fileName,
      uint start_ind=0);

  template <typename type>
  void importDat(
      std::vector< std::vector<type> > &import,
      std::string fileName,
      uint start_ind=0);

  template <typename type>
  void importDat(
      std::vector< std::vector< std::vector<type> > > &import,
      std::string fileName,
      uint start_ind=0);


  ////////////////////
  /////  Extras  /////
  ////////////////////

  std::vector<int> getShape(std::string folder, std::string filePrefix);

}


template <typename type>
void save::saveDat(type* input, int size, std::string fileName) {
  cout << "INFO: Saving " + fileName + " ... ";
  FILE* output = fopen(fileName.c_str(), "wb");
  fwrite(input, sizeof(type), size, output);
  fclose(output);
  cout << "Done!" << endl;

  return;
}


template <typename type>
void save::saveDat(std::vector<type> &input, std::string fileName) {
  cout << "INFO: Saving " + fileName + " ... ";
  FILE* output = fopen(fileName.c_str(), "wb");
  fwrite(&input[0], sizeof(type), input.size(), output);
  fclose(output);
  cout << "Done!" << endl;

  return;
}


template <typename type>
void save::saveDat(std::vector< std::vector<type> > &input, std::string fileName) {
  cout << "INFO: Saving " + fileName + " ... ";
  FILE* output = fopen(fileName.c_str(), "wb");
  //cout<<"111"<<endl;
  for (uint i=0; i<input.size(); i++) {
    //cout<<"\ti "<<i<<endl;
    fwrite(&input[i][0], sizeof(type), input[i].size(), output);
  }
  //cout<<"222"<<endl;
  fclose(output);
  cout << "Done!" << endl;

  return;
}


template <typename type>
void save::saveDat(std::vector< std::vector< std::vector<type> > > &input,
    std::string fileName) {
  cout << "INFO: Saving " + fileName + " ... ";
  FILE* output = fopen(fileName.c_str(), "wb");
  for (uint i=0; i<input.size(); i++) {
    for (uint ii=0; ii<input[i].size(); ii++) {
      fwrite(&input[i][ii][0], sizeof(type), input[i][ii].size(), output);
    }
  }
  fclose(output);
  cout << "Done!" << endl;

  return;
}


template <typename type>
void save::importDat(
    std::vector<type> &import,
    std::string fileName,
    uint start_ind) {

  // Check that vector is not empty
  if (!import.size()) {
    std::cerr << "ERROR: Cannot fill vector of size 0!!!\n";
    exit(1);
  }

  FILE* input = fopen(fileName.c_str(), "rb");
  if (input == NULL) {
    std::cerr << "ERROR: Cannot open file " + fileName << 
        " with error: " << strerror(errno) << "!!!\n";
    exit(1);
  }

  fread(&import[start_ind], sizeof(type), import.size(), input);

  fclose(input);
}


template <typename type>
void save::importDat(
    std::vector< std::vector<type> > &import,
    std::string fileName,
    uint start_ind) {

  // Check that vector is not empty
  if (!import.size()) {
    std::cerr << "ERROR: Cannot fill vector of size 0!!!\n";
    exit(1);
  }

  FILE* input = fopen(fileName.c_str(), "rb");
  if (input == NULL) {
    std::cerr << "ERROR: Cannot open file " + fileName << "!!!\n";
    exit(1);
  }

  if (start_ind > import.size()) {
    std::cerr << "ERROR: start_ind > import.size()!!!\n";
    exit(1);
  }

  uint Nrows = (uint)((int)import.size() - (int)start_ind);
  uint Ncols = import[0].size();
  std::vector<type> inpVec(Nrows*Ncols);
  fread(&inpVec[0], sizeof(type), inpVec.size(), input);

  fclose(input);

  for (int ir=0; ir<Nrows; ir++) {
    for (int ic=0; ic<Ncols; ic++) {
      import[start_ind+ir][ic] = inpVec[ir*Ncols + ic];
    }
  }
}


template <typename type>
void save::importDat(
    std::vector< std::vector< std::vector<type> > > &import,
    std::string fileName,
    uint start_ind) {

  // Check that vector is not empty
  if (!import.size()) {
    std::cerr << "ERROR: Cannot fill vector of size 0!!!\n";
    exit(1);
  }

  FILE* input = fopen(fileName.c_str(), "rb");
  if (input == NULL) {
    std::cerr << "ERROR: Cannot open file " + fileName << "!!!\n";
    exit(1);
  }
  
  if (start_ind > import.size()) {
    std::cerr << "ERROR: start_ind > import.size()!!!\n";
    exit(1);
  }

  uint Nrows = (uint)((int)import.size() - (int)start_ind);
  uint Ncols = import[0].size();
  uint Ntrds = import[0][0].size();
  std::vector<type> inpVec(Nrows*Ncols*Ntrds);
  fread(&inpVec[0], sizeof(type), inpVec.size(), input);

  fclose(input);

  for (uint ir=0; ir<Nrows; ir++) {
    for (uint ic=0; ic<Ncols; ic++) {
      for (uint it=0; it<Ntrds; it++) {
        import[start_ind+ir][ic][it] = inpVec[ir*Ncols*Ntrds + ic*Ntrds + it];
      }
    }
  }
}


#endif
