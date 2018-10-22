#include "saveClass.h"

using namespace std;


SAVEclass::SAVEclass() {

  file = NULL;
  diff = NULL;
}


SAVEclass::SAVEclass(string filename, string treename) {

  newFile(filename, treename);
  diff = NULL;
}


SAVEclass::~SAVEclass() {

  if (file) {
    writeFile();
  }
}


void SAVEclass::newFile(string filename, string treename) {

  file = TFile::Open(filename.c_str(),"recreate");
  tree = new TTree(treename.c_str(), treename.c_str());
  diff = NULL;
}


void SAVEclass::writeFile() {

  if (file) {
    file->Write();
    file->Close();
    file = NULL;
  }
}



void SAVEclass::addInfo(DIFFRACTIONclass* diff_in) {

  diff = diff_in;

  tree->Branch("diffPattern", &diff->diffPattern);
  tree->Branch("diffMolPattern", &diff->diffMolPattern);
  tree->Branch("diffAtmPattern", &diff->diffAtmPattern);
}


void SAVEclass::addInfo(vector<string> names, vector<double*> ptrs) {

  for (uint i=0; i<names.size(); i++) {
    tree->Branch(names[i].c_str(), ptrs[i], (names[i]+"/D").c_str());
  }
}


void SAVEclass::fillEvent() {

  tree->Fill();
}



