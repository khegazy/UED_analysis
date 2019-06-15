#ifndef SAVING_H
#define SAVING_H

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
  void importDat(std::vector<type> &import, std::string fileName);

  template <typename type>
  void importDat(std::vector< std::vector<type> > &import, std::string fileName);


  ////////////////////
  /////  Extras  /////
  ////////////////////

  std::vector<int> getShape(std::string folder, std::string filePrefix);

}


template <typename type>
void save::saveDat(type* input, int size, std::string fileName) {
  FILE* output = fopen(fileName.c_str(), "wb");
  fwrite(input, sizeof(type), size, output);
  fclose(output);

  return;
}


template <typename type>
void save::saveDat(std::vector<type> &input, std::string fileName) {
  FILE* output = fopen(fileName.c_str(), "wb");
  fwrite(&input[0], sizeof(type), input.size(), output);
  fclose(output);

  return;
}


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


template <typename type>
void save::importDat(std::vector<type> &import, std::string fileName) {
  // Check that vector is not empty
  if (!import.size()) {
    std::cerr << "ERROR: Cannot fill vector of size 0!!!\n";
    exit(0);
  }

  FILE* input = fopen(fileName.c_str(), "rb");
  if (input == NULL) {
    std::cerr << "ERROR: Cannot open file " + fileName << "!!!\n";
    exit(0);
  }

  fread(&import[0], sizeof(type), import.size(), input);

  fclose(input);
}


template <typename type>
void save::importDat(std::vector< std::vector<type> > &import, std::string fileName) {
  // Check that vector is not empty
  if (!import.size()) {
    std::cerr << "ERROR: Cannot fill vector of size 0!!!\n";
    exit(0);
  }

  FILE* input = fopen(fileName.c_str(), "rb");
  if (input == NULL) {
    std::cerr << "ERROR: Cannot open file " + fileName << "!!!\n";
    exit(0);
  }

  int Nrows = import.size();
  int Ncols = import[0].size();
  std::vector<type> inpVec(Nrows*Ncols);
  fread(&inpVec[0], sizeof(type), inpVec.size(), input);

  fclose(input);

  for (int ir=0; ir<Nrows; ir++) {
    for (int ic=0; ic<Ncols; ic++) {
      import[ir][ic] = inpVec[ir*Ncols + ic];
    }
  }
}


#endif
