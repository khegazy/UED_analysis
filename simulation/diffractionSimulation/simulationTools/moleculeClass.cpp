#include "moleculeClass.h"

using namespace std;

MOLECULEclass::MOLECULEclass() {

  initialize(true);
}


MOLECULEclass::MOLECULEclass(bool deleteAtoms_in) {

  initialize(deleteAtoms_in);
}


MOLECULEclass::~MOLECULEclass() {

//cout<<"in mol dest"<<endl;
  if (deleteAtoms) {
//cout<<"deleteing atoms"<<endl;
    for (int iatm=0; iatm<atoms.size(); iatm++) delete atoms[iatm];
  }

//cout<<"done deleteing atoms"<<endl;
}


void MOLECULEclass::initialize(bool deleteAtoms_in) {

  atoms.clear();
  deleteAtoms = deleteAtoms_in;
}


void MOLECULEclass::addAtom(ATOMclass* atom) {

  atoms.push_back(atom);
}
